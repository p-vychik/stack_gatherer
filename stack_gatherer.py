# start a continuous process from CLI
# Every 100 ms check whether all planes for a full stack have been saved
# We determine whether the stack is done, if new images with a higher timepoint/specimen index,
# or if more than 30 sec has passed since last plane
#
import math
import argparse
import sys
import shutil
import json
import csv
import time
import os
import re
import numpy as np
import napari
import multiprocessing
import traceback
import requests
import threading
import pandas as pd
import matplotlib.pyplot as plt
import cv2
import logging
from tifffile import imread, imwrite, memmap
from tqdm import tqdm
from PIL import Image
import asyncio
from watchfiles import awatch
from threading import Thread
from queue import Queue
from collections import OrderedDict
from qtpy.QtWidgets import QCheckBox
from dataclasses import dataclass
from datetime import datetime
from multiprocessing import Manager, Event

import zmq
from PlottingWindow import Ui_MainWindow
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QMainWindow
from scipy.signal import find_peaks

CURRENT_ACQUISITION_LOGGER = None
STOP_FILE_NAME = "STOP_STACK_GATHERING"
ACQUISITION_META_FILE_PATTERN = "AcquisitionMetadata_"
# BORDER_WIDTH defines the size of the border in pixels between
# projections merged in one image for napari visualization
BORDER_WIDTH = 20
PLT_WIDGET_X_AXIS_LENGTH = 40
PROJECTIONS_QUEUE = OrderedDict()
DRAWN_PROJECTIONS_QUEUE = OrderedDict()
JSON_CONFIGS = {}
MERGE_LIGHT_MODES = False
# napari viewer window
VIEWER = None
# separate matplotlib windows for each specimen
PLT_WIDGETS_DICT = {}
TEMP_DIR = ""
SPECIMENS_QUANTITY_LOADED = False
SPECIMENS_QUANTITY = []
PLOTTING_WINDOW_CREATED = False
plotting_windows_timer = QTimer()
active_stacks = {}
currenty_saving_stacks_locks = {}
currenty_saving_stacks_dict_lock = threading.Lock()
HEARTBEAT_INTERVAL_SEC = 5
command_queue_to_microscope = Queue()
new_microscope_command_event = threading.Event()
HEARTBEAT_FROM_MICROSCOPE_TIMEOUT_sec = 1000

LAYERS_INPUT = multiprocessing.Queue()
PIV_OUTPUT = multiprocessing.Queue()


def setup_main_logger(output_dir):
    log_filename = os.path.join(output_dir, f"stack_gatherer_{datetime.now().strftime('%Y-%b-%d-%H%M%S')}.log")
    logging.basicConfig(filename=log_filename,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%d-%b-%y %H:%M:%S',
                        level=logging.INFO)


def setup_acquisition_log(lapse_output_folder: str):
    global CURRENT_ACQUISITION_LOGGER

    if not os.path.exists(lapse_output_folder):
        os.mkdir(lapse_output_folder)
    if CURRENT_ACQUISITION_LOGGER:
        for handler in CURRENT_ACQUISITION_LOGGER.handlers:
            handler.close()
            CURRENT_ACQUISITION_LOGGER.removeHandler(handler)

    CURRENT_ACQUISITION_LOGGER = logging.getLogger(f"acquisition_log")
    log_filename = os.path.join(lapse_output_folder, f"{os.path.basename(lapse_output_folder).strip('.json')}.log")

    file_handler = logging.FileHandler(log_filename)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    file_handler.setFormatter(formatter)

    CURRENT_ACQUISITION_LOGGER.addHandler(file_handler)
    CURRENT_ACQUISITION_LOGGER.setLevel(logging.INFO)


def logging_broadcast(string: str):
    print(string, file=sys.stderr)
    if CURRENT_ACQUISITION_LOGGER:
        CURRENT_ACQUISITION_LOGGER.info(string)
    else:
        logging.getLogger('main_logger').info(string)
    # logging.info(string)


class PlottingWindow(QMainWindow, Ui_MainWindow):
    window_index = 0

    def __init__(self, window_index):
        super().__init__()
        # Initialize GUI
        self.setupUi(self)
        self.window_index = window_index
        self.setWindowTitle(f"QuickPIV average speed, specimen #{self.window_index}")

    def closeEvent(self, event):
        self.close()

    def plot_curve(self, x, y, x_marker, y_marker):
        # Plot a curve:
        self.plotWidget.canvas.axes.clear()
        y_lim_range = (math.floor(min(y) * 0.9), math.ceil(max(y) * 1.1))
        self.plotWidget.canvas.axes.set_ylim(y_lim_range)
        self.plotWidget.canvas.axes.plot(x, y, color='blue')
        if x_marker and y_marker:
            self.plotWidget.canvas.axes.plot(x_marker, y_marker,
                                             color='red', marker='o',
                                             markersize=12, linewidth=0)
            for xm, ym in zip(x_marker, y_marker):
                self.plotWidget.canvas.axes.annotate(f"{xm}",
                                                     xy=(xm, ym),
                                                     xytext=(5, 5),
                                                     textcoords='offset points',
                                                     fontsize=10,
                                                     color='black')
        self.plotWidget.canvas.axes.set_xticks(x)
        self.plotWidget.canvas.axes.set_xticklabels(x, fontsize=10)
        self.plotWidget.canvas.draw()


class PivProcess(multiprocessing.Process):
    def __init__(self, *args, **kwargs):
        multiprocessing.Process.__init__(self, *args, **kwargs)
        self._pconn, self._cconn = multiprocessing.Pipe()
        self._exception = None

    def run(self):
        try:
            multiprocessing.Process.run(self)
            self._cconn.send(None)
        except Exception as e:
            tb = traceback.format_exc()
            self._cconn.send((e, tb))

    @property
    def exception(self):
        if self._pconn.poll():
            self._exception = self._pconn.recv()
            logging_broadcast(self._exception[0])
        return self._exception


class NotImagePlaneFile(Exception):
    pass


@dataclass
class StackSignature:
    total_num_planes: int
    timelapse_id: str
    specimen: int
    time_point: int
    illumination: int
    camera: int
    channel: int

    def get_attr_excluding(self, exclude_fields=None) -> tuple:
        if isinstance(exclude_fields, tuple) or isinstance(exclude_fields, list):
            return tuple([val for attribute, val in self.__dict__.items() if attribute not in exclude_fields])
        if exclude_fields:
            return tuple([val for attribute, val in self.__dict__.items() if attribute != exclude_fields])
        return (self.total_num_planes, self.timelapse_id, self.specimen, self.time_point,
                self.illumination, self.camera, self.channel)

    def __hash__(self) -> hash:
        return hash(self.get_attr_excluding())

    @property
    def signature(self) -> tuple:
        return self.get_attr_excluding()

    @property
    def signature_no_time(self) -> tuple:
        return self.get_attr_excluding(exclude_fields=["time_point"])


@dataclass
class ProjectionsDictWrapper:
    """
    identifier property defines the lapse mode parameters used except the illumination;
    illumination property is a value that identify the light mode used for a plane;
    signature property is a tuple with values representing parameters used for lapse and stored as object
    of the StackSignature class
    """
    projections: dict
    stack_signature_obj: StackSignature

    @property
    def identifier(self) -> tuple:
        return self.stack_signature_obj.get_attr_excluding("illumination")

    @property
    def illumination(self) -> int:
        return self.stack_signature_obj.illumination

    @property
    def signature(self) -> tuple:
        return self.stack_signature_obj.get_attr_excluding()

    @property
    def stack_signature(self) -> StackSignature:
        return self.stack_signature_obj

    @property
    def specimen(self) -> int:
        return self.stack_signature_obj.specimen


class ImageFile:
    """
    File naming for light-sheet image files.
    Example file name: SPC-0001_TP-0001_ILL-0_CAM-0_CH-01_PL-0001-outOf-0150_blaBla.tif or .bmp
    :param file_path: full path to image
    :type file_path: str
    """
    timelapse_id = ""
    time_point = 0
    specimen = 0
    illumination = 0
    camera = 0
    channel = 0
    plane = 0
    total_num_planes = 0
    additional_info = ""
    extension = ""
    path_to_image_dir = ""

    def __init__(self, file_path):
        self.path_to_image_dir = os.path.dirname(file_path)
        file_name = os.path.basename(file_path)
        split_by = r"timelapseID-|SPC-|TP-|ILL-|CAM-|CH-|PL-|outOf-|\."
        name_parts = re.split(split_by, file_name)

        if len(name_parts) == 10:
            try:
                for i, name_part in enumerate(name_parts):

                    if i == 0:
                        self.dataset_name = name_part.strip("-_")
                    elif i == 1:
                        self.timelapse_id = name_part.strip("-_")
                    elif i == 2:
                        self.specimen = int(name_part.strip("-_"))
                    elif i == 3:
                        self.time_point = int(name_part.strip("-_"))
                    elif i == 4:
                        self.illumination = int(name_part.strip("-_"))
                    elif i == 5:
                        self.camera = int(name_part.strip("-_"))
                    elif i == 6:
                        self.channel = int(name_part.strip("-_"))
                    elif i == 7:
                        self.plane = int(name_part.strip("-_"))
                    elif i == 8:
                        name_and_info = name_part.strip("-_").split("_", 1)
                        if len(name_and_info) == 1:
                            self.total_num_planes = int(name_and_info[0].strip("-_"))
                        else:
                            num_planes, info = name_and_info
                            self.total_num_planes = int(num_planes.strip("-_"))
                            self.additional_info = info
                    elif i == 9:
                        self.extension = name_part
            except ValueError:
                raise NotImagePlaneFile(
                    "This is not a valid filename for a single plane image file!"
                )
        else:
            raise NotImagePlaneFile(
                "Image file name is improperly formatted! Check documentation inside the script. "
                "Expected 10 parts after splitting by %s" % split_by)

        self.extension.lower()

    def get_name(self):
        additional_info = self.additional_info
        dataset_name = self.dataset_name
        if additional_info != "":
            additional_info = "_" + additional_info
        if dataset_name != "":
            dataset_name = dataset_name + "_"
        return (f"{dataset_name}timelapseID-{self.timelapse_id:}_SPC-{self.specimen:04}"
                f"_TP-{self.time_point:04}_ILL-{self.illumination}"
                f"_CAM-{self.camera}_CH-{self.channel:02}"
                f"_PL-{self.plane:04}-outOf-{self.total_num_planes:04}{additional_info}.{self.extension}"
                )

    def get_name_without_extension(self):
        return os.path.splitext(self.get_name())[0]

    def get_stack_name(self):
        additional_info = self.additional_info
        dataset_name = self.dataset_name
        if additional_info != "":
            additional_info = "_" + additional_info
        if dataset_name != "":
            dataset_name = dataset_name + "_"
        return (f"{dataset_name}timelapseID-{self.timelapse_id:}_SPC-{self.specimen:04}"
                f"_TP-{self.time_point:04}_ILL-{self.illumination}"
                f"_CAM-{self.camera}_CH-{self.channel:02}"
                f"_PL-(ZS)-outOf-{self.total_num_planes:04}{additional_info}.{self.extension}"
                )

    def get_stack_path(self):
        return os.path.join(self.path_to_image_dir, self.get_stack_name())

    def get_file_path(self):
        return os.path.join(self.path_to_image_dir, self.get_name())

    def get_stack_signature(self):
        return StackSignature(self.total_num_planes,
                              self.timelapse_id,
                              self.specimen,
                              self.time_point,
                              self.illumination,
                              self.camera,
                              self.channel)


def read_image(image_path):
    if image_path.endswith((".tif", ".tiff")):
        time.sleep(0.1)
        try:
            image = imread(image_path)
            return image
        except Exception as error:
            logging_broadcast(error)
    if image_path.endswith(".bmp"):
        try:
            image = Image.open(image_path)
            if image.mode == "RGB":
                width, height = image.size
                with open(image_path, 'rb') as img_bytes:
                    data = img_bytes.read()
                    header_len = len(data) - height * width * 2
                    no_header_data = data[header_len:]
                    image = Image.frombytes('I;16', (width, height), no_header_data, 'raw')
                return np.asarray(image, dtype="uint16")
            return np.array(image)
        except Exception as error:
            logging_broadcast(error)
    return False


def file_name_merged_illumination_based_on_signature(stack_signature):
    return (f"timelapseID-{stack_signature.timelapse_id}_"
            f"SPC-{stack_signature.specimen}_"
            f"TP-{stack_signature.time_point}_"
            f"ILL-MERGED_CAM-{stack_signature.camera}_"
            f"CH-{stack_signature.channel}_"
            f"Z_MAX_projection")


def fix_image_drift(ref_img, img):
    # ORB detector
    orb = cv2.ORB_create()

    # detect keypoints and descriptors in images
    kp_ref, des_ref = orb.detectAndCompute(ref_img, None)
    kp_img, des_img = orb.detectAndCompute(img, None)
    bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)
    matches = bf.match(des_img, des_ref)
    matches = sorted(matches, key=lambda x: x.distance)

    # extract matches
    points_ref = np.zeros((len(matches), 2), dtype=np.float32)
    points_img = np.zeros((len(matches), 2), dtype=np.float32)

    for i, match in enumerate(matches):
        points_ref[i, :] = kp_ref[match.trainIdx].pt
        points_img[i, :] = kp_img[match.queryIdx].pt

    # find homography
    h, mask = cv2.findHomography(points_img, points_ref, cv2.RANSAC)

    # use homography to warp image
    height, width = ref_img.shape
    aligned_img = cv2.warpPerspective(img, h, (width, height))

    return aligned_img


def plane_to_projection(plane: np.ndarray, output_dictionary: dict):
    """
    Function gets as an input a plane as a numpy array and a dictionary that keys define the projections to be made.
    As the function output, a modified dictionary is returned that values are numpy arrays with transformed matrices.
    for the projection specified by the dictionary key.
    Final result is acquired only after all planes are processed sequentially.

    Define the dimension to collapse for appropriate projection:
    Y projection - combining max values in the rows for each plane
        in the new matrix of size (number of planes X number of rows in a plane)
    X projection - combining max values in the columns (gravity axis) for each plane
        in the new matrix of size (number of planes X number of columns in a plane)
    Z projection - new matrix of the same size as source,
        each value represents a max value among all values along the same depth-axis
    """

    axes_codes = {"Y": 1, "X": 0, "Z": 0}
    axes = output_dictionary.keys()
    for axis in axes:
        array = output_dictionary[axis]
        if axis == "Y" or axis == "X":
            if isinstance(array, np.ndarray):
                array = np.vstack((array, plane.max(axis=axes_codes[axis])))
            else:
                array = plane.max(axis=axes_codes[axis])
        elif axis == "Z":
            if isinstance(array, np.ndarray):
                array = np.stack((array, plane))
                array = array.max(axis=axes_codes[axis])
            else:
                array = plane
        output_dictionary[axis] = array
    return output_dictionary


def collect_files_to_one_stack_get_axial_projections(stack_signature: StackSignature,
                                                     file_list,
                                                     output_file_path,
                                                     shared_dict,
                                                     unprocessed_stack_signatures,
                                                     axes=None,
                                                     anisotropy_factor=1,
                                                     output_dir=None
                                                     ):
    if os.path.exists(output_file_path):
        if stack_signature.signature not in unprocessed_stack_signatures:
            logging_broadcast(f"For image stack {stack_signature.signature} the output file {output_file_path} "
                              f"already exists, processing skipped to prevent corruption "
                              f"of existing files in the folder")
            unprocessed_stack_signatures.add(stack_signature.signature)
        return False
    sample_image = read_image(file_list[0])
    if isinstance(sample_image, bool):
        return False
    shape = (len(file_list), sample_image.shape[0], sample_image.shape[1])
    dtype = sample_image.dtype
    # create an empty OME-TIFF file
    imwrite(output_file_path, shape=shape, dtype=dtype, metadata={'axes': 'ZYX'})
    # memory map numpy array to data in OME-TIFF file
    zyx_stack = memmap(output_file_path)
    logging_broadcast(f"Writing stack to {output_file_path}")
    # prepare input about required projections in the dictionary,
    # the values of the appropriate keys would be the projection matrices
    projections = {}
    projections_files_path = {}
    if axes:
        projections = {axis.upper(): None for axis in set(axes) if axis in "zZxXyY"}
    if projections:
        try:
            file_name = os.path.basename(file_list[0]).replace('.bmp', '.tif')
            for axis in projections.keys():
                projection_folder_path = os.path.join(output_dir, f"{axis}_projections")
                if not os.path.exists(projection_folder_path):
                    os.makedirs(projection_folder_path)
                projections_files_path[axis] = os.path.join(projection_folder_path, file_name)
        except OSError as err:
            logging_broadcast(err)
            return False
    # write data to memory-mapped array
    with tqdm(total=len(file_list), desc="Saving plane") as pbar:
        ref_img = sample_image
        for z in range(shape[0]):
            if z == 0:
                zyx_stack[z] = sample_image
                if shape[0] == 1:
                    # stack which consists of one plane means that on input are Z projections,
                    # clear dictionary to skip processing
                    projections.clear()
                    projections["Z"] = read_image(file_list[z])
                    zyx_stack.flush()
                    pbar.update(1)
                    break
                zyx_stack.flush()
                pbar.update(1)
                continue
            zyx_stack[z] = read_image(file_list[z])
            # doesn't work yet!
            # use BFMatcher to correct shift which may be present
            # aligned_img = fix_image_drift(ref_img, zyx_stack[z])
            # zyx_stack[z] = aligned_img
            # ref_img = aligned_img
            projections = plane_to_projection(zyx_stack[z], projections)
            zyx_stack.flush()
            pbar.update(1)
    if projections:
        # correct anisotropy for X and Y projections by resizing the array with user specified factor
        if "Y" in projections:
            projection_y = projections["Y"]
            projection_y = projection_y.transpose()

            img = Image.fromarray(projection_y)
            projections["Y"] = np.array(img.resize(size=(projections["Y"].shape[0] * anisotropy_factor,
                                                         projections["Y"].shape[1])))
        if "X" in projections:
            projection_x = projections["X"]
            img = Image.fromarray(projection_x)
            projections["X"] = np.array(img.resize(size=(projections["X"].shape[1],
                                                         projections["X"].shape[0] * anisotropy_factor))
                                        )
        if "Z" in projections:
            # push to the queue for PIV calculations
            wrapped_z_projection = ProjectionsDictWrapper({"Z": projections["Z"]}, stack_signature)
            # here we should put the projection with clear identification of species and other parameters
            if stack_signature.specimen in shared_dict:
                shared_dict[stack_signature.specimen].put(wrapped_z_projection)
            else:
                logging_broadcast(f"Specimen signature not found in shared queue, signature: {stack_signature}")
        for axis in projections.keys():
            try:
                imwrite(projections_files_path[axis], projections[axis])
            except Exception as err:
                logging_broadcast(err)

        # we want to store projections derived from planes with different illumination modes to
        # have an option to merge them if the user checked the QCheckBox() in napari viewer GUI interface
        wrapped_projections = ProjectionsDictWrapper(projections, stack_signature)
        if wrapped_projections.signature not in PROJECTIONS_QUEUE:
            PROJECTIONS_QUEUE[wrapped_projections.identifier] = [wrapped_projections]
        else:
            PROJECTIONS_QUEUE[wrapped_projections.identifier].append(wrapped_projections)

    for z in range(shape[0]):
        try:
            os.remove(file_list[z])
        except PermissionError as e:
            logging_broadcast(f"Error: {e}")
    return True


def add_file_to_active_stacks(image_file: ImageFile):
    stack_signature = image_file.get_stack_signature()
    if stack_signature.signature not in active_stacks:
        active_stacks[stack_signature.signature] = {}
        logging_broadcast(f"Adding stack {stack_signature.signature} to active queue.")
    if image_file.plane not in active_stacks[stack_signature.signature]:
        active_stacks[stack_signature.signature][image_file.plane] = image_file
    return stack_signature


def check_stack_and_collect_if_ready(stack_signature: StackSignature,
                                     output_dir,
                                     shared_dict,
                                     unprocessed_stack_signatures,
                                     axes=None,
                                     factor=None
                                     ):
    if len(active_stacks[stack_signature.signature]) < stack_signature.total_num_planes:
        return
    # We have to ensure that two events firing at the same time don't start saving the same stack twice
    with currenty_saving_stacks_dict_lock:
        if stack_signature not in currenty_saving_stacks_locks:
            currenty_saving_stacks_locks[stack_signature.signature] = True
        else:
            return
    file_list = []
    for i, _ in enumerate(active_stacks[stack_signature.signature]):
        # We have to access by index since we can't guarantee that files were added to dict in order of planes
        file_list.append(active_stacks[stack_signature.signature][i].get_file_path())
    sample_file_obj = ImageFile(file_list[0])
    sample_file_obj.extension = "tif"
    stack_path = os.path.join(output_dir, sample_file_obj.get_stack_name())
    if collect_files_to_one_stack_get_axial_projections(stack_signature,
                                                        file_list,
                                                        stack_path,
                                                        shared_dict,
                                                        unprocessed_stack_signatures,
                                                        axes=axes,
                                                        anisotropy_factor=factor,
                                                        output_dir=output_dir
                                                        ):

        global LAYERS_INPUT
        images_dict = get_projections_dict_from_queue()
        if images_dict:
            LAYERS_INPUT.put(images_dict)
        del active_stacks[stack_signature.signature]
        del currenty_saving_stacks_locks[stack_signature.signature]
    else:
        del currenty_saving_stacks_locks[stack_signature.signature]



def load_lapse_parameters_json(file_path: str,
                               destination_folder: str):
    try:
        j_file = open(file_path)
        lapse_parameters = json.load(j_file)
        j_file.close()
    except PermissionError:
        logging_broadcast(f"{file_path} permission error, check if file is already opened")
        return False
    except Exception as e:
        logging_broadcast(
            f"{file_path} failed to parse JSON", e)

    if lapse_parameters:
        timestamp = os.path.basename(file_path).strip(".json").split("_")[1]
        lapse_parameters["output_folder"] = destination_folder
        setup_signature = []
        for spec_entry in lapse_parameters["specimens"]:
            specimen_number = spec_entry["userDefinedIndex"]
            total_num_planes = spec_entry["number_of_planes"]
            for channel in spec_entry["channels"]:
                channel_num = channel["userDefinedIndex"]
                active_illum = channel["lightsheetsEnabled"]
                active_illum = [idx for idx, value in enumerate(active_illum) if value]
                for illum in active_illum:
                    for camera_num, enabled in enumerate(channel["camerasEnabled"]):
                        if enabled:
                            setup_signature.append((total_num_planes,
                                                    timestamp,
                                                    specimen_number,
                                                    illum,
                                                    camera_num,
                                                    channel_num))
        return setup_signature, lapse_parameters
    return False


def load_file_from_input_folder(file_path,
                                output_dir,
                                shared_dict,
                                json_dict,
                                manager,
                                unprocessed_stack_signatures,
                                factor=1,
                                axes=""
                                ):
    global SPECIMENS_QUANTITY
    global SPECIMENS_QUANTITY_LOADED
    global PLOTTING_WINDOW_CREATED

    if os.path.isdir(file_path):
        return
    # load json parameters
    if os.path.basename(file_path).startswith(ACQUISITION_META_FILE_PATTERN) and file_path.endswith(".json"):
        destination_folder = os.path.join(output_dir, os.path.basename(file_path).strip(".json"))
        setup_acquisition_log(destination_folder)
        json_parsing_result = load_lapse_parameters_json(file_path, destination_folder)
        if not json_parsing_result:
            return
        list_of_stack_signatures, metadata_json = json_parsing_result
        json_dict.update(dict.fromkeys(
            list_of_stack_signatures, metadata_json))
        # get specimens quantity
        try:
            SPECIMENS_QUANTITY = [specimen_dict['userDefinedIndex'] for specimen_dict in metadata_json['specimens']]
            SPECIMENS_QUANTITY_LOADED = True
        except Exception as err:
            logging_broadcast(err)
        if len(SPECIMENS_QUANTITY):
            shared_dict.clear()
            for i in SPECIMENS_QUANTITY:
                shared_dict[i] = manager.Queue()
        PLOTTING_WINDOW_CREATED = False
        try:
            if not os.path.exists(os.path.join(destination_folder, os.path.basename(file_path))):
                shutil.move(file_path, destination_folder)
        except Exception as error:
            logging_broadcast(error)
    if not file_path.upper().endswith((".TIF", ".BMP", ".TIFF")):
        return
    # Call the function when a new file is created
    try:
        file = ImageFile(file_path)
    except NotImagePlaneFile:
        # image file is incorrect, move it to the temp_dir
        try:
            shutil.move(file_path, TEMP_DIR)
        except OSError:
            pass
        return
    stack_signature = add_file_to_active_stacks(file)
    # create separate folder for the output based on metadata filename
    try:
        output_dir = json_dict[stack_signature.signature_no_time]["output_folder"]
        if os.path.exists(output_dir):
            check_stack_and_collect_if_ready(stack_signature,
                                             output_dir,
                                             shared_dict,
                                             unprocessed_stack_signatures,
                                             axes,
                                             factor
                                             )
        else:
            logging_broadcast(f"{output_dir} doesn't exist, check if the output folder is correct")
    except KeyError:
        if stack_signature.signature not in unprocessed_stack_signatures:
            unprocessed_stack_signatures.add(stack_signature.signature)
            logging_broadcast(f"No lapse configuration for the plane {stack_signature.signature}, "
                              f"check if AcquisitionMeta file is loaded")


def update_napari_viewer_layer():
    global LAYERS_INPUT
    try:
        if LAYERS_INPUT.qsize() > 0:
            data_input = LAYERS_INPUT.get()
            for axes_names, image in data_input.items():
                if axes_names not in VIEWER.layers:
                    VIEWER.add_image(image, name=axes_names)
                    VIEWER.add_shapes(name=f"Subsample {axes_names} for quickPIV", shape_type='rectangle')
                    VIEWER.layers[f"Subsample {axes_names} for quickPIV"].mouse_drag_callbacks.append(on_rectangle_selection)
                else:
                    layer_image = VIEWER.layers[axes_names]
                    VIEWER.layers[axes_names].data = image
                    if image.dtype == np.dtype('uint16') and layer_image.contrast_limits[-1] <= 255:
                        VIEWER.layers[axes_names].reset_contrast_limits()
                    if image.dtype == np.dtype('uint8') and layer_image.contrast_limits[-1] > 255:
                        VIEWER.layers[axes_names].reset_contrast_limits()
    except Exception as err:
        logging_broadcast(err)



def find_config_files_locations(output_folder: str,
                                input_folder: str,
                                json_file_name: str,
                                unprocessed_lapse_ids: set) -> dict:
    """
    Returns the locations of time-lapse metadata files in JSON format based on corresponding timestamped names
    of image planes from the batch for all files in the output_folder and input_folder.
    Alternatively, if the metadata file name is specified with the json_file_name argument,
    proceed with a direct check for its existence in the output_folder and input_folder

    Args:
    - output_folder (str): Path to the output folder to search for the JSON files.
    - input_folder (str): Path to the input folder to search for the JSON files.
    - json_file_name (str, optional): Name of the JSON file to find. If provided,
      search will be restricted to finding this specific file.
    - unprocessed_lapse_ids (set): Set of lapse ids which were not processed due
    to the absence of correct config file, used for preventing abundant warning messages

    Returns:
    - dict: A dictionary mapping lapse IDs to the full paths of the JSON files found.
            If json_file_name is provided, returns a dictionary with a single entry.
    """
    config_files_locations = {}

    # If json_file_name is provided, directly check its presence
    if json_file_name:
        if os.path.exists(os.path.join(output_folder, json_file_name)):
            config_files_locations[json_file_name] = os.path.join(output_folder, json_file_name)
        elif os.path.exists(os.path.join(output_folder, os.path.splitext(json_file_name)[0], json_file_name)):
            config_files_locations[json_file_name] = os.path.join(output_folder,
                                                                  os.path.splitext(json_file_name)[0],
                                                                  json_file_name)
        elif os.path.exists(os.path.join(input_folder, json_file_name)):
            config_files_locations[json_file_name] = os.path.join(input_folder, json_file_name)
        if not len(config_files_locations):
            logging_broadcast(f"{json_file_name} was not found, images associated with this lapse won't be processed")
        return config_files_locations

    # If json_file_name is not provided, continue with lapse ID-based search
    content = os.listdir(input_folder)
    if len(content) == 0:
        return config_files_locations

    files = [f for f in content if os.path.isfile(os.path.join(input_folder, f))]
    lapse_ids = set()

    for file_name in files:
        if file_name.upper().endswith("JSON"):
            continue
        try:
            lapse_id = file_name.split("_SPC")[0].split("timelapseID-")[1]
            lapse_ids.add(lapse_id)
        except IndexError:
            logging_broadcast(f"problem with {file_name} tame-lapse id processing, it's not"
                              f"flanked by 'timelapseID-' and '_SPC'")

    for lapse_id in lapse_ids:
        if os.path.exists(
                os.path.join(output_folder,
                             f"{ACQUISITION_META_FILE_PATTERN}{lapse_id}",
                             f"{ACQUISITION_META_FILE_PATTERN}{lapse_id}.json")):
            config_files_locations[lapse_id] = os.path.join(output_folder,
                                                            f"{ACQUISITION_META_FILE_PATTERN}{lapse_id}",
                                                            f"{ACQUISITION_META_FILE_PATTERN}{lapse_id}.json")
        elif os.path.exists(os.path.join(output_folder,
                                         f"{ACQUISITION_META_FILE_PATTERN}{lapse_id}.json")):
            config_files_locations[lapse_id] = os.path.join(output_folder,
                                                            f"{ACQUISITION_META_FILE_PATTERN}{lapse_id}.json")
        elif os.path.exists(os.path.join(input_folder,
                                         f"{ACQUISITION_META_FILE_PATTERN}{lapse_id}.json")):
            config_files_locations[lapse_id] = os.path.join(input_folder,
                                                            f"{ACQUISITION_META_FILE_PATTERN}{lapse_id}.json")
        # all known locations are checked for config file without success, display the warning once
        elif lapse_id not in unprocessed_lapse_ids:
            unprocessed_lapse_ids.add(lapse_id)
            logging_broadcast(
                f"Neither the output nor input directories contain the timelapse configuration file "
                f"for lapse id {lapse_id}. Processing of planes won't be initiated.")

    return config_files_locations


async def read_input_files(input_folder,
                           output_dir,
                           shared_queues_of_z_projections,
                           json_config,
                           manager,
                           factor,
                           axes,
                           exit_gracefully: threading.Event,
                           acquisition_metadata_to_process,
                           process_z_projections
                           ):
    try:
        global SPECIMENS_QUANTITY
        global SPECIMENS_QUANTITY_LOADED
        awatch_input_folder = awatch(input_folder)
        unprocessed_lapse_ids = set()
        unprocessed_stack_signatures = set()
        while True:
            try:
                if exit_gracefully.is_set():
                    break

                if process_z_projections:
                    specimen_indices = set()
                    for f in os.listdir(input_folder):
                        if not os.path.isfile(os.path.join(input_folder, f)):
                            continue
                        match = re.search(r"(?<=SPC-)\d+", f)
                        if match:
                            # Convert the matched string to an integer
                            specimen_indices.add(int(match.group()))
                    if specimen_indices:
                        SPECIMENS_QUANTITY.extend(sorted(specimen_indices))
                        SPECIMENS_QUANTITY_LOADED = True

                config_files = find_config_files_locations(output_dir,
                                                           input_folder,
                                                           acquisition_metadata_to_process,
                                                           unprocessed_lapse_ids)
                acquisition_metadata_to_process = ""
                if config_files:
                    for _, path in config_files.items():
                        try:
                            load_file_from_input_folder(path,
                                                        output_dir,
                                                        shared_queues_of_z_projections,
                                                        json_config,
                                                        manager,
                                                        unprocessed_stack_signatures,
                                                        factor,
                                                        axes
                                                        )
                        except Exception as e:
                            logging_broadcast(f"Error processing file {path}: {e}")
                    files = []
                    for f in os.listdir(input_folder):
                        if not os.path.isfile(os.path.join(input_folder, f)):
                            continue
                        if f.endswith(".json"):
                            continue
                        files.append(os.path.join(input_folder, f))
                    files.sort()
                    for path in files:
                        try:
                            load_file_from_input_folder(path,
                                                        output_dir,
                                                        shared_queues_of_z_projections,
                                                        json_config,
                                                        manager,
                                                        unprocessed_stack_signatures,
                                                        factor,
                                                        axes
                                                        )
                        except Exception as e:
                            logging_broadcast(f"Error processing file {path}: {e}")


                changes = await asyncio.wait_for(awatch_input_folder.__anext__(), timeout=5.0)
                paths = []
                for change in changes:
                    event, file_path = change
                    if event.value == 1:  # Assuming '1' means 'file created' or similar
                        paths.append(file_path)
                if paths:
                    paths.sort()
                    for path in paths:
                        try:
                            load_file_from_input_folder(path,
                                                        output_dir,
                                                        shared_queues_of_z_projections,
                                                        json_config,
                                                        manager,
                                                        unprocessed_stack_signatures,
                                                        factor,
                                                        axes
                                                        )
                        except Exception as e:
                            logging_broadcast(f"Error processing file {path}: {e}")
            except asyncio.TimeoutError:
                if exit_gracefully.is_set():
                    break
            except StopAsyncIteration:
                awatch_input_folder = awatch(input_folder)
                continue
    except Exception as err:
        exit_gracefully.set()
        logging_broadcast(err)


def watchfiles_thread(*args):
    asyncio.run(read_input_files(*args))
    logging_broadcast("watchfiles finished")


def run_the_loop(kwargs, exit_gracefully: threading.Event):
    global TEMP_DIR
    global SPECIMENS_QUANTITY_LOADED
    global SPECIMENS_QUANTITY
    global CROP_COORDINATES

    input_dir = kwargs.get("input")
    output_dir = kwargs.get("output")
    axes = kwargs.get("axes", "Z")
    factor = kwargs.get("factor_anisotropy", None)
    TEMP_DIR = kwargs.get("temp_dir", None)
    bot_config_path = kwargs.get("bot_config", "")
    process_z_projections = kwargs.get("process_z_projections", False)
    acquisition_metadata_to_process = kwargs.get("restart", "")
    token, chat_id = "", ""
    quickPIV_parameters = kwargs.get("quickPIV_parameters", "")
    multiquickPIV_jl_path = kwargs.get("multiquickpivjl", "")

    piv_processes = []
    manager = Manager()
    json_config = manager.dict()
    shared_queues_of_z_projections = manager.dict()
    cropping_coordinates_napari_selection = manager.dict()
    stop_file = os.path.join(input_dir, STOP_FILE_NAME)

    # Start the input directory observer
    logging_broadcast(f"Watching {input_dir} for images, and saving stacks to {output_dir}")
    thread = Thread(target=watchfiles_thread, args=(input_dir,
                                                    output_dir,
                                                    shared_queues_of_z_projections,
                                                    json_config,
                                                    manager,
                                                    factor,
                                                    axes,
                                                    exit_gracefully,
                                                    acquisition_metadata_to_process,
                                                    process_z_projections),
                    daemon=False
                    )
    thread.start()
    if multiquickPIV_jl_path:
        quickpiv_run_file = [
            'using PythonCall',
            f"""include("{multiquickPIV_jl_path}")""",
            'function fn(a1, a2)',
            f"   {quickPIV_parameters}",
            '   Us = []',
            '   Vs = []',
            '   pivparams = multi_quickPIV.setPIVParameters(interSize=IA, searchMargin=SM, '
            'step=ST, corr_alg="nsqecc", threshold=TH)',
            '   VF, SN = multi_quickPIV.PIV(a1, a2, pivparams, precision=64)',
            '   push!(Us, VF[1, :, :]);',
            '   push!(Vs, VF[2, :, :])',
            '   avg_speed = [sum(sqrt.(Us[t] .^ 2 + Vs[t] .^ 2)) / length(Us[t]) for t in 1:length(Us)]',
            '   return avg_speed',
            'end']
        quick_piv_runner = os.path.join(output_dir, "avg_speed_quickPIV.jl")
        try:
            if os.path.exists(quick_piv_runner):
                os.remove(quick_piv_runner)
            with open(f"{quick_piv_runner}", 'w') as f:
                f.write(('\n').join(quickpiv_run_file))
        except Exception as err:
            logging_broadcast(err)

    while True:
        migration_detected = manager.Queue()
        avg_speed_data = manager.list()
        stop_process = Event()

        if multiquickPIV_jl_path:
            while not SPECIMENS_QUANTITY_LOADED and not exit_gracefully.is_set():
                time.sleep(1)

            if process_z_projections:
                for i in SPECIMENS_QUANTITY:
                    shared_queues_of_z_projections[i] = manager.Queue()

            if bot_config_path:
                with open(bot_config_path) as fin:
                    try:
                        line = fin.readline()
                        token, chat_id = line.strip().split(" ")
                    except Exception as err:
                        logging_broadcast(err)

            for i, _ in shared_queues_of_z_projections.items():
                piv_process = PivProcess(
                    target=run_piv_process, args=(shared_queues_of_z_projections,
                                                  PIV_OUTPUT,
                                                  avg_speed_data,
                                                  stop_process,
                                                  quick_piv_runner,
                                                  migration_detected,
                                                  i,
                                                  process_z_projections,
                                                  output_dir
                                                  ),
                    name='piv_run',
                )
                piv_processes.append(piv_process)
                piv_process.start()

            # mode to process only max projections
            if process_z_projections:
                logging_broadcast(f"Considering directory as Z projections source: {input_dir}, calculate"
                                  f" average speed and save data and plots to {output_dir}")
                z_projections = os.listdir(input_dir)
                z_projections.sort()
                for projection in z_projections:
                    file_path = os.path.join(input_dir, projection)
                    img_data = read_image(file_path)
                    img_metadata = ImageFile(file_path)
                    if not isinstance(img_data, bool) and not isinstance(img_metadata, bool):
                        stack_signature = img_metadata.get_stack_signature()
                        wrapped_z_projection = ProjectionsDictWrapper({"Z": img_data}, stack_signature)
                        # push projection to the queue for PIV calculations
                        if stack_signature.specimen in shared_queues_of_z_projections:
                            while True:
                                if shared_queues_of_z_projections[stack_signature.specimen].qsize() > 50:
                                    time.sleep(1)
                                else:
                                    break
                            shared_queues_of_z_projections[stack_signature.specimen].put(wrapped_z_projection)
                logging_broadcast(f"Finished with uploading projection files for quickPIV")

        try:

            check_next = time.time() + 60
            speed_list_length = len(avg_speed_data)
            while not os.path.exists(stop_file) and SPECIMENS_QUANTITY_LOADED and not exit_gracefully.is_set():
                if migration_detected.qsize() > 0:
                    migration_event = migration_detected.get()
                    message = f"Detected on {migration_event}!"
                    if token and chat_id:
                        url = (f"https://api.telegram.org/bot{token}/"
                               f"sendMessage?chat_id={chat_id}&text={message}&disable_web_page_preview=true")
                        response = requests.get(url)
                        if response.status_code != 200:
                            logging_broadcast(response.text)
                            logging_broadcast(url)
                    logging_broadcast(message)
                if process_z_projections:
                    # check the data length change approximately every minute
                    # if the avg_speed_data doesn't change - save the data and finish the pipeline
                    if check_next - time.time() <= 0:
                        if speed_list_length == len(avg_speed_data):
                            logging_broadcast(f"No changes in quickPIV output queue, "
                                              f"stop pipeline and save the data")
                            break
                        check_next = time.time() + 60
                        speed_list_length = len(avg_speed_data)
                # Sleep to keep the script running
                time.sleep(1)
            if process_z_projections:
                stop_process.set()
        except KeyboardInterrupt:
            # Gracefully stop the observer if the script is interrupted
            stop_process.set()

        # terminate running quickPIV processes
        if piv_processes:
            for process in piv_processes:
                process.terminate()
                time.sleep(0.1)
                process.join()
            time.sleep(0.5)
            piv_processes = []

        # save PIV data to csv
        if avg_speed_data:
            sorted_speed_data = sorted(avg_speed_data, key=lambda d: d['specimen'])
            current_time = datetime.now()
            with open(os.path.join(output_dir, f"quickPIV_data_{current_time.strftime('%H_%M_%S')}.csv"),
                      'w', newline='') as f_out:
                try:
                    # sort results by specimen index
                    w = csv.DictWriter(f_out, sorted_speed_data[0].keys())
                    w.writeheader()
                    w.writerows(sorted_speed_data)
                    logging_broadcast(f"csv data saved to {output_dir}")
                except IOError:
                    logging_broadcast(f"Attempt to save csv data to {output_dir} failed")
            if process_z_projections:
                df = pd.DataFrame(sorted_speed_data)
                specimen_indices = df["specimen"].unique()
                for specimen_index in specimen_indices:
                    df_subsample = df[df["specimen"] == specimen_index]
                    x, y = df_subsample["time_point"], df_subsample["avg_speed"]
                    x = x.tolist()
                    y = y.tolist()
                    plt.rcParams['figure.figsize'] = [40, 20]
                    fig, ax = plt.subplots()
                    ax.plot(x, y)
                    y = np.asarray(y)
                    peaks, _ = find_peaks(y, prominence=1.5)
                    ax.plot(peaks, y[peaks], color='red', marker='o', markersize=12, linewidth=0)
                    for xm, ym in zip(peaks, y[peaks]):
                        ax.annotate(f"{xm}",
                                    xy=(xm, ym),
                                    xytext=(5, 5),
                                    textcoords='offset points',
                                    fontsize=10,
                                    color='black')
                    dec_x = []
                    for i, num in enumerate(x, start=1):
                        if i % 100 == 0:
                            dec_x.append(num)
                    ax.xaxis.set_ticks(dec_x)
                    ax.set_xticklabels(dec_x, fontsize=12, rotation=0)
                    ax.set(xlabel='frame', ylabel='Avg. speed')
                    plt.title(f"Specimen {specimen_index}")
                    try:
                        plt.savefig(os.path.join(output_dir,
                                                 f"Specimen_{specimen_index}_avg_speed_"
                                                 f"{current_time.strftime('%H_%M_%S')}.png"))
                    except Exception as err:
                        logging_broadcast(err)
        if os.path.exists(stop_file) or process_z_projections or exit_gracefully.is_set():
            exit_gracefully.set()
            break
    thread.join()


def draw_napari_layer(projections_dict):
    # define the padding between pictures
    # zero arrays as vertical and horizontal borders
    # projections_dict = wrapped_projections_dict.projections
    image = None
    image_dict = {}
    dtype_def = "uint8"
    axes_names = ''.join([axis_key for axis_key, _ in projections_dict.items()])
    for _, projection in projections_dict.items():
        dtype_def = projection.dtype
        break
    if len(projections_dict) == 3:
        v_border_array = np.zeros((projections_dict["Z"].shape[0], BORDER_WIDTH), dtype=dtype_def)
        # zero array for horizontal border
        # rows number equals to the BORDER_WIDTH value, columns to width of concatenation of Z, Y and border array
        h_border_array = np.zeros((BORDER_WIDTH,
                                   projections_dict["Z"].shape[1] + BORDER_WIDTH + projections_dict["Y"].shape[1]),
                                  dtype=dtype_def
                                  )
        # extend Z projection with border and Y projection arrays
        z_y = np.hstack((projections_dict["Z"], v_border_array, projections_dict["Y"]))
        # merge Z_Y with horizontal border array
        z_y = np.vstack((z_y, h_border_array))
        x = np.hstack((projections_dict["X"],
                       np.zeros((projections_dict["X"].shape[0], projections_dict["Y"].shape[1] + BORDER_WIDTH),
                                dtype=dtype_def
                                ))
                      )
        image = np.vstack((z_y, x))
    elif len(projections_dict) == 2:
        # place the largest projection in center and arrange the second at the appropriate side
        # if only X and Y - then arrange them in a perpendicular way
        if "Z" in projections_dict and "X" in projections_dict:
            h_border_array = np.zeros((BORDER_WIDTH, projections_dict["Z"].shape[1]), dtype=dtype_def
                                      )
            image = np.vstack((projections_dict["Z"], h_border_array, projections_dict["X"])
                              )
        elif "Z" in projections_dict and "Y" in projections_dict:
            v_border_array = np.zeros((projections_dict["Z"].shape[0], BORDER_WIDTH), dtype=dtype_def
                                      )
            image = np.hstack((projections_dict["Z"], v_border_array, projections_dict["Y"])
                              )
        else:
            # only X and Y projections, arrange them perpendicular
            dummy_array = np.zeros((projections_dict["Y"].shape[0], projections_dict["X"].shape[1]),
                                   dtype=dtype_def
                                   )
            dummy_y = np.hstack((dummy_array, projections_dict["Y"]))
            x_extended = np.hstack((projections_dict["X"],
                                    np.zeros((projections_dict["X"].shape[0],
                                              dummy_y.shape[1] - projections_dict["X"].shape[1]), dtype=dtype_def
                                             ))
                                   )
            image = np.vstack((dummy_y, x_extended))
    elif len(projections_dict) == 1:
        image = projections_dict[axes_names]
    image_dict[axes_names] = image
    return image_dict


def merge_multiple_projections(wrapped_dict_list: tuple):
    # as an input a list of wrapped dictionaries is provided - check class ProjectionsDictWrapper
    # we have to merge appropriate projections according to the maximum intensity
    merged_projections = {}
    for wrapped_dict in wrapped_dict_list[-1]:
        for axis, plane in wrapped_dict.projections.items():
            if axis not in merged_projections:
                merged_projections[axis] = plane
            else:
                array = np.stack((merged_projections[axis], plane))
                # get max projection
                merged_projections[axis] = array.max(axis=0)
    return merged_projections


def get_projections_dict_from_queue():
    global DRAWN_PROJECTIONS_QUEUE
    if len(PROJECTIONS_QUEUE) > 0:
        image_layer_dict = {}
        # store last 4 projections drawn in napari layer for merging channels purpose
        if len(DRAWN_PROJECTIONS_QUEUE) > 4:
            DRAWN_PROJECTIONS_QUEUE.popitem(last=False)
        identifier, projections_dict_list = PROJECTIONS_QUEUE.popitem(last=False)
        if MERGE_LIGHT_MODES:
            # should check the number of illuminations
            if len(projections_dict_list) == 2:
                merged_projections_dict = merge_multiple_projections((identifier, projections_dict_list))
                image_layer_dict = draw_napari_layer(merged_projections_dict)
            elif identifier in DRAWN_PROJECTIONS_QUEUE and projections_dict_list:
                projections_dict_list.extend(DRAWN_PROJECTIONS_QUEUE[identifier])
                merged_projections_dict = merge_multiple_projections((identifier, projections_dict_list))
                image_layer_dict = draw_napari_layer(merged_projections_dict)
                del DRAWN_PROJECTIONS_QUEUE[identifier]
            elif projections_dict_list:
                for wrapped_dict in projections_dict_list:
                    channel_layer = draw_napari_layer(wrapped_dict.projections)
                    for key, val in channel_layer.items():
                        image_layer_dict[(f"Specimen_{wrapped_dict.specimen}_{key}"
                                          f"_ILL_{wrapped_dict.illumination}")] = val
                DRAWN_PROJECTIONS_QUEUE[identifier] = list(projections_dict_list)
        else:
            if projections_dict_list:
                # iterate through the list, assign illumination index to the axis, save all to the new ,
                # data remain unmerged - apply merge and push to napari
                for wrapped_dict in projections_dict_list:
                    channel_layer = draw_napari_layer(wrapped_dict.projections)
                    for key, val in channel_layer.items():
                        image_layer_dict[(f"Specimen_{wrapped_dict.specimen}_{key}"
                                          f"_ILL_{wrapped_dict.illumination}")] = val
                DRAWN_PROJECTIONS_QUEUE[identifier] = list(projections_dict_list)
            else:
                logging_broadcast("Empty projections dictionary")
        if image_layer_dict:
            return image_layer_dict


def bx_trigger():
    global MERGE_LIGHT_MODES
    MERGE_LIGHT_MODES = not MERGE_LIGHT_MODES


def extract_and_display_selection(min_coords, max_coords):
    try:
        for layer in VIEWER.layers:
            if isinstance(layer, napari.layers.Image):
                image_data = layer.data
                selected_region = image_data[
                    min_coords[0]:max_coords[0],
                    min_coords[1]:max_coords[1]
                ]
                # VIEWER.add_image(selected_region, name=f'Selection from {layer.name}')
                imwrite(f"F:/cropping_data/{layer.name}_{time.time()}.tiff", selected_region)
                break
    except Exception as err:
        logging_broadcast(f"extract_and_display_{err}")


def on_rectangle_selection(layer, event):
    shapes_layer = layer
    # Check if the current tool is the rectangle tool
    logging_broadcast(f"mouse release fired {shapes_layer.name}")
    logging_broadcast("release mode fired")
    if len(shapes_layer.data) > 0 and event.type == "mouse_press" and shapes_layer.mode == "add_rectangle":
        rectangle = shapes_layer.data.pop()  # Get the last drawn rectangle
        shapes_layer.data = []
        shapes_layer.data.append(rectangle)
        logging_broadcast(f"{rectangle}")
        min_coords = np.min(rectangle, axis=0).astype(int)
        max_coords = np.max(rectangle, axis=0).astype(int)
        extract_and_display_selection(min_coords, max_coords)


def make_napari_viewer():
    napari_viewer = napari.Viewer()
    bx = QCheckBox('Merge illumination channels')
    bx.setChecked(MERGE_LIGHT_MODES)
    bx.stateChanged.connect(bx_trigger)
    napari_viewer.window.add_dock_widget(bx)
    # Add shapes layer for rectangle selection
    # shapes_layer = napari_viewer.add_shapes(name='Rectangle', shape_type='rectangle')
    # shapes_layer.mode = 'add_rectangle'
    # shapes_layer.events.data.connect(on_rectangle_selection)
    return napari_viewer


def run_piv_process(shared_dict_queue: dict,
                    queue_out: multiprocessing.Queue,
                    avg_speed_data: list,
                    stop_process: Event,
                    piv_path: str,
                    migration_detected: multiprocessing.Queue,
                    queue_number: int,
                    save_merged_illumination: bool,
                    output_dir: str
                    ):
    from juliacall import Main as jl
    jl.include(piv_path)
    piv_projection_queue = Queue()
    projections_to_process = Queue()
    migration_event_frame = set()
    queue_in = shared_dict_queue[queue_number]
    logging_broadcast(f"quickPIV process started, PID: {os.getpid()}")
    t = 0
    x = []
    y = []
    while not stop_process.is_set():
        plane_matrix = queue_in.get()
        if plane_matrix:
            projections_to_process.put(plane_matrix)
        if projections_to_process.qsize() >= 2:
            wrapped_z_p_1 = projections_to_process.get()
            # to calculate avg speed in consecutive way we store
            # the last projection from every pairwise comparison
            wrapped_z_p_2 = projections_to_process.queue[0]
            if wrapped_z_p_1.identifier == wrapped_z_p_2.identifier:
                merged_projections = merge_multiple_projections((wrapped_z_p_1.identifier,
                                                                 (wrapped_z_p_1,
                                                                  wrapped_z_p_2))
                                                                )
                merged_wrapped_proj = ProjectionsDictWrapper(merged_projections,
                                                             wrapped_z_p_2.signature)
                if save_merged_illumination:
                    try:
                        file_name = file_name_merged_illumination_based_on_signature(wrapped_z_p_2.stack_signature)
                        imwrite(os.path.join(output_dir, f"{file_name}.tif"), merged_wrapped_proj.projections["Z"])
                    except Exception as err:
                        migration_detected.put(f"Failed to save the merged Z-projection, {err}")
                piv_projection_queue.put(merged_wrapped_proj)
                projections_to_process.get()
            else:
                piv_projection_queue.put(wrapped_z_p_1)

            if piv_projection_queue.qsize() < 2:
                continue
            m_1 = piv_projection_queue.get().projections["Z"]
            # to compare projections in a consecutive way,
            # we have to leave the last projection in the piv_projection_queue
            m_2 = piv_projection_queue.queue[0].projections["Z"]
            if m_1.shape == m_2.shape:
                try:
                    avg_speed = jl.fn(m_1, m_2)
                    avg_speed = round(avg_speed[-1], 3)
                except Exception as error:
                    raise error
                # current_time = now.strftime("%H:%M:%S")
                t += 1
                avg_speed_data.append({"time_point": t, "specimen": queue_number, "avg_speed": avg_speed})
                try:
                    y.append(avg_speed)
                    x.append(t)
                except Exception as error:
                    raise error
                if len(x) > PLT_WIDGET_X_AXIS_LENGTH:
                    x = x[-PLT_WIDGET_X_AXIS_LENGTH::]
                    y = y[-PLT_WIDGET_X_AXIS_LENGTH::]
                try:
                    peaks, _ = find_peaks(np.asarray(y), prominence=1.5)
                except Exception as error:
                    raise error
                x_marker, y_marker = [], []
                if len(peaks):
                    try:
                        peaks = peaks.tolist()
                        x_marker = [x[i] for i in peaks]
                        y_marker = [y[i] for i in peaks]
                        if len(x) > PLT_WIDGET_X_AXIS_LENGTH:
                            global_peaks_indices = [x + t - PLT_WIDGET_X_AXIS_LENGTH + 1 for x in x_marker]
                        else:
                            global_peaks_indices = x_marker
                        for frame in global_peaks_indices:
                            if frame not in migration_event_frame:
                                migration_event_frame.add(frame)
                                migration_detected.put(f"{frame}, specimen {queue_number}")
                    except Exception as error:
                        raise error
                queue_out.put((queue_number, x, y, x_marker, y_marker))
            else:
                raise ValueError("Projections should have the same size for QuickPIV input")


def update_avg_speed_plot_windows():
    global PIV_OUTPUT
    if PIV_OUTPUT.qsize() > 0:
        index, x, y, x_marker, y_marker = PIV_OUTPUT.get()
        try:
            PLT_WIDGETS_DICT[index].plot_curve(x, y, x_marker, y_marker)
        except KeyError:
            logging_broadcast(f"window with index {index} doesn't exist")


def update_plotting_windows(exit_gracefully: threading.Event,
                            multiquickpivjl: str):
    global PLOTTING_WINDOW_CREATED
    global PLT_WIDGETS_DICT
    global plotting_windows_timer

    if exit_gracefully.is_set():
        if len(PLT_WIDGETS_DICT) > 0:
            for _, window in PLT_WIDGETS_DICT.items():
                window.close()
            plotting_windows_timer.stop()

    if PLOTTING_WINDOW_CREATED or not multiquickpivjl:
        return
    if len(PLT_WIDGETS_DICT) > 0:
        # delete window from previous lapse
        PLT_WIDGETS_DICT = {}
    for i in SPECIMENS_QUANTITY:
        PLT_WIDGETS_DICT[i] = PlottingWindow(i)
        PLT_WIDGETS_DICT[i].show()
    PLOTTING_WINDOW_CREATED = True


def parse_message_from_microscope(message, exit_gracefully: threading.Event):
    time.sleep(0.01)
    try:
        if message.get("type") == "exit":
            logging_broadcast("Received terminate command from microscope.")
            exit_gracefully.set()
        elif message.get("type") == "heartbeat":
            # logging_broadcast("Received heartbeat from the microscope.")
            pass
        else:
            logging_broadcast(f"Received message from the microscope: {str(message)}\nBut there is no action for it.")
    except:
        pass


def heartbeat_and_command_handler(port, exit_gracefully: threading.Event):
    logging_broadcast(f"heartbeat process started")
    context = zmq.Context()
    socket = context.socket(zmq.PAIR)
    socket.connect(f"tcp://localhost:{port}")
    socket.setsockopt(zmq.LINGER, 0)  # Set linger to 0 to prevent blocking on close

    poller = zmq.Poller()
    poller.register(socket, zmq.POLLIN)

    last_heartbeat_sent = time.time() - HEARTBEAT_INTERVAL_SEC
    last_heartbeat_recieved = time.time()  # initialize value with current time

    try:
        while not exit_gracefully.is_set():
            next_heartbeat_time = last_heartbeat_sent + HEARTBEAT_INTERVAL_SEC
            time_until_heartbeat = max(0, next_heartbeat_time - time.time())
            # Timeout in milliseconds, max 50ms
            timeout = min(50, time_until_heartbeat * 1000)

            events = dict(poller.poll(timeout))
            if socket in events:
                message = socket.recv_json()
                # logging_broadcast("Received message: " + str(message))
                parse_message_from_microscope(message, exit_gracefully)
                if message.get("type") == "heartbeat":
                    last_heartbeat_recieved = time.time()

            if time.time() - last_heartbeat_recieved > HEARTBEAT_FROM_MICROSCOPE_TIMEOUT_sec:
                exit_gracefully.set()

            current_time = time.time()
            if current_time >= next_heartbeat_time:
                # logging_broadcast("Sending heartbeat to the microscope.")
                socket.send_json({"type": "heartbeat", "message": "alive"})
                last_heartbeat_sent = current_time

            while not command_queue_to_microscope.empty():
                command = command_queue_to_microscope.get_nowait()
                socket.send_json({"type": "command", "message": command})

    except Exception as e:
        logging_broadcast(f"Heartbeat handler exception: {e}")
        exit_gracefully.set()
    finally:
        logging_broadcast("disconnecting communication socket")
        socket.close()
        poller.unregister(socket)
        context.term()
        logging_broadcast("sending heartbeat stopped")


def correct_argv(argv):
    corrected_argv = []
    i = 0
    while i < len(argv):
        # Check presence of Windows drive letter followed by colon and quotation mark
        # absence of the double slash after the colon results in incorrect arguments parsing
        # when all arguments after the occurrence of quotation mark are parsed as one.

        # look for presence of substrings E:" or E:\\" with a regex
        if re.match(r'^[A-Za-z]:\"|^[A-Za-z]:\\{1,}\"', argv[i]):
            corrected_argv.append(argv[i][:2])
            concatenated_args = re.split(r'\s+', argv[i])
            argv.extend(concatenated_args[1:])
        else:
            corrected_argv.append(argv[i])
        i += 1
    return corrected_argv


def close_napari_viewer(exit_gracefully: threading.Event):
    if exit_gracefully.is_set():
        VIEWER.close()


def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Watch a directory for new file additions and collect .tif or .bmp "
                                                 "files as stacks to output directory.",
                                     argument_default=argparse.SUPPRESS)
    # Define the input_folder argument
    parser.add_argument('-i', '--input', required=True,
                        help="Input folder path to watch.")

    # Define the output_folder argument
    parser.add_argument('-o', '--output', required=True,
                        help="Output folder path to save stacks.")

    # Define axes to make projections
    parser.add_argument('-a', '--axes', required=False,
                        help="Comma separated axes to project on the plane, i.e. X,Y,Z or X, or X")
    # Define port for communication with MicroscopeController software
    parser.add_argument('-p', '--port', required=False, help="Integer value of the port")

    # Anisotropy factor correction
    parser.add_argument('-f', '--factor_anisotropy', type=int,
                        required=True if '--axes' in sys.argv else False,
                        help="Value is used for correcting projection's anisotropic distortions")
    # Move image files with incorrect name to the user provided directory
    parser.add_argument('--temp_dir', required=False,
                        default=None,
                        help="Directory path to store input images with incorrect file name")
    parser.add_argument('--multiquickpivjl', required=True if '-process_z_projections' in sys.argv else False,
                        default=False,
                        help="Path to the multi_quickPIV.jl from quickPIV repository")
    parser.add_argument('--bot_config', required=False, default=False,
                        help="Path to text file with token and chat_id information (one line, space separator) for"
                             " supplying peaks detection on telegram messenger")
    parser.add_argument('-debug_run', required=False, default=False, action='store_true',
                        help="All content of --output specified folder will be DELETED!")
    parser.add_argument('-process_z_projections',
                        required=False, action='store_true',
                        help="Process input files as Z-projections for average speed calculations with quickPIV!")
    # When the microscope controller restarts the stack_gatherer, this argument defines the configuration file
    # of the unfinished image batch with which the processing will start
    parser.add_argument('-r', '--restart', required=False, default=None,
                        help=" defines the name of the AcquisitionMetadata file for unfinished image batch "
                             "with which the processing will start")
    parser.add_argument('--quickPIV_parameters', required=True if '--multiquickpivjl'
                                                                  in sys.argv else False,
                        help="Define IA,SM,ST,TH parameters (more about the arguments in QuickPIV documentation) "
                             "as the string, without spaces e.g. 'IA,SM,ST,TH=(32,32),(16,16),(16,16),30'")

    # event to stop the script gracefully
    exit_gracefully = threading.Event()
    # Parse the command-line arguments
    corrected_args = correct_argv(sys.argv)
    args = parser.parse_args(corrected_args[1:])

    if args.temp_dir is None:
        args.temp_dir = args.output
    setup_main_logger(args.output)
    logging.getLogger('watchfiles').setLevel(logging.WARNING)
    if args.debug_run:
        # for testing purposes only - empty the output folder before loading files
        for filename in os.listdir(args.output):
            file_path = os.path.join(args.output, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                logging_broadcast(e)
    global VIEWER, PLT_WIDGETS_DICT, plotting_windows_timer

    # Allow other computers to attach to debugpy at this IP address and port.
    # debugpy.listen(('0.0.0.0', 5680))
    # logging_broadcast("Waiting for debugger attach")
    # debugpy.wait_for_client()  # Blocks execution until client is attached

    try:
        heartbeat_thread = None
        heartbeat_thread = threading.Thread(
            target=heartbeat_and_command_handler, args=(args.port, exit_gracefully))
        heartbeat_thread.start()
    except Exception as e:
        logging_broadcast(f"{args},\n {e}")
    thread = Thread(target=run_the_loop, args=(vars(args), exit_gracefully))
    thread.start()
    VIEWER = make_napari_viewer()
    timer1 = QTimer()
    timer1.timeout.connect(update_avg_speed_plot_windows)
    timer1.start(25)
    timer2 = QTimer()
    timer2.timeout.connect(update_napari_viewer_layer)
    timer2.start(50)
    plotting_windows_timer.timeout.connect(lambda: update_plotting_windows(exit_gracefully, args.multiquickpivjl))
    plotting_windows_timer.start(1000)
    timer4 = QTimer()
    timer4.timeout.connect(lambda: close_napari_viewer(exit_gracefully))
    timer4.start(1000)
    napari.run()
    logging_broadcast(f"Napari viewer windows was closed, terminating child processes")
    exit_gracefully.set()
    thread.join()
    if heartbeat_thread is not None:
        heartbeat_thread.join()
    logging_broadcast(f"stack_gatherer stopped")


if __name__ == "__main__":
    main()
