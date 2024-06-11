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

# import debugpy

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
PIVJLPATH = ""
SPECIMENS_QUANTITY_LOADED = False
SPECIMENS_QUANTITY = list()
PLOTTING_WINDOW_CREATED = False
plotting_windows_timer = QTimer()
active_stacks = {}
currenty_saving_stacks_locks = {}
currenty_saving_stacks_dict_lock = threading.Lock()
HEARTBEAT_INTERVAL_SEC = 5
command_queue_to_microscope = Queue()
new_microscope_command_event = threading.Event()


LAYERS_INPUT = multiprocessing.Queue()
PIV_OUTPUT = multiprocessing.Queue()


def debugprint(*args):
    print(*args, file=sys.stderr)
# Define the main window class


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
            debugprint(self._exception[0])
        return self._exception


class NotImagePlaneFile(Exception):
    pass


@dataclass
class StackSignature:
    total_num_planes: int
    timelapse_id: str
    time_point: int
    specimen: int
    illumination: int
    camera: int
    channel: int

    def get_attr_excluding(self, exclude_fields=None) -> tuple:
        if isinstance(exclude_fields, tuple) or isinstance(exclude_fields, list):
            return tuple([val for attribute, val in self.__dict__.items() if attribute not in exclude_fields])
        if exclude_fields:
            return tuple([val for attribute, val in self.__dict__.items() if attribute != exclude_fields])
        return (self.total_num_planes, self.timelapse_id, self.time_point, self.specimen,
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


class ImageFile:
    """
    File naming for light-sheet image files.
    Example file name: TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001-outOf-0150_blaBla.tif or .bmp
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
        split_by = r"timelapseID-|TP-|SPC-|ILL-|CAM-|CH-|PL-|outOf-|\."
        name_parts = re.split(split_by, file_name)

        if len(name_parts) == 10:
            try:
                for i, name_part in enumerate(name_parts):

                    if i == 0:
                        self.dataset_name = name_part.strip("-_")
                    elif i == 1:
                        self.timelapse_id = name_part.strip("-_")
                    elif i == 2:
                        self.time_point = int(name_part.strip("-_"))
                    elif i == 3:
                        self.specimen = int(name_part.strip("-_"))
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
        return (f"{dataset_name}timelapseID-{self.timelapse_id:}_TP-{self.time_point:04}"
                f"_SPC-{self.specimen:04}_ILL-{self.illumination}"
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
        return (f"{dataset_name}timelapseID-{self.timelapse_id:}_TP-{self.time_point:04}"
                f"_SPC-{self.specimen:04}_ILL-{self.illumination}"
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
                              self.time_point,
                              self.specimen,
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
            debugprint(error)
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
            debugprint(error)
    return False


def file_name_merged_illumination_based_on_signature(stack_signature):
    return (f"timelapseID-{stack_signature.timelapse_id}_"
            f"TP-{stack_signature.time_point}_"
            f"SPC-{stack_signature.specimen}_"
            f"ILL-MERGED_CAM-{stack_signature.camera}_"
            f"CH-{stack_signature.channel}_"
            f"Z_MAX_projection")


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
                                                     axes=None,
                                                     anisotropy_factor=1,
                                                     output_dir=None):
    if os.path.exists(output_file_path):
        return
    sample_image = read_image(file_list[0])
    if isinstance(sample_image, bool):
        return
    shape = (len(file_list), sample_image.shape[0], sample_image.shape[1])
    dtype = sample_image.dtype
    # create an empty OME-TIFF file
    imwrite(output_file_path, shape=shape, dtype=dtype, metadata={'axes': 'ZYX'})
    # memory map numpy array to data in OME-TIFF file
    zyx_stack = memmap(output_file_path)
    debugprint(f"Writing stack to {output_file_path}")
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
            debugprint(err)
    # write data to memory-mapped array
    with tqdm(total=len(file_list), desc="Saving plane") as pbar:
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
                debugprint(f"Specimen signature not found in shared queue, signature: {stack_signature}")
                # debugprint(shared_dict)
        for axis in projections.keys():
            try:
                imwrite(projections_files_path[axis], projections[axis])
            except Exception as err:
                debugprint(err)

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
            debugprint(f"Error: {e}")


def add_file_to_active_stacks(image_file: ImageFile):
    stack_signature = image_file.get_stack_signature()
    if stack_signature.signature not in active_stacks:
        active_stacks[stack_signature.signature] = {}
        debugprint(f"Adding stack {stack_signature.signature} to active queue.")
    if image_file.plane not in active_stacks[stack_signature.signature]:
        active_stacks[stack_signature.signature][image_file.plane] = image_file
    return stack_signature


def check_stack_and_collect_if_ready(stack_signature: StackSignature, output_dir, shared_dict, axes=None, factor=None):
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
    collect_files_to_one_stack_get_axial_projections(stack_signature,
                                                     file_list,
                                                     stack_path,
                                                     shared_dict,
                                                     axes=axes,
                                                     anisotropy_factor=factor,
                                                     output_dir=output_dir
                                                     )
    global LAYERS_INPUT
    images_dict = get_projections_dict_from_queue()
    if images_dict:
        LAYERS_INPUT.put(images_dict)
    del active_stacks[stack_signature.signature]
    del currenty_saving_stacks_locks[stack_signature.signature]


def load_lapse_parameters_json(file_path: str,
                               destination_folder: str):
    try:
        j_file = open(file_path)
        lapse_parameters = json.load(j_file)
        j_file.close()
    except PermissionError:
        debugprint(f"{file_path} permission error, check if file is already opened")
        return False
    except Exception as e:
        debugprint(
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


def load_file_from_input_folder(file_path, output_dir, shared_dict, json_dict, manager, factor=1, axes=""):
    global SPECIMENS_QUANTITY
    global SPECIMENS_QUANTITY_LOADED
    global PLOTTING_WINDOW_CREATED

    if os.path.isdir(file_path):
        return
    # load json parameters
    if os.path.basename(file_path).startswith(ACQUISITION_META_FILE_PATTERN) and file_path.endswith(".json"):
        destination_folder = os.path.join(output_dir, os.path.basename(file_path).strip(".json"))
        json_parsing_result = load_lapse_parameters_json(
            file_path, destination_folder)
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
            debugprint(err)
        if SPECIMENS_QUANTITY:
            shared_dict.clear()
            for i in SPECIMENS_QUANTITY:
                shared_dict[i] = manager.Queue()
        PLOTTING_WINDOW_CREATED = False
        try:
            if not os.path.exists(destination_folder):
                os.mkdir(destination_folder)
                shutil.move(file_path, destination_folder)
        except Exception as error:
            debugprint(error)
    else:
        if not file_path.endswith((".tif", ".bmp", ".tiff")):
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
                                                 axes,
                                                 factor)
        except KeyError:
            debugprint("No lapse configuration for the plane, check if AcquisitionMeta file is loaded")


def update_napari_viewer_layer():
    global LAYERS_INPUT
    if LAYERS_INPUT.qsize() > 0:
        data_input = LAYERS_INPUT.get()
        for axes_names, image in data_input.items():
            if axes_names not in VIEWER.layers:
                VIEWER.add_image(image, name=axes_names)
            else:
                layer_image = VIEWER.layers[axes_names]
                VIEWER.layers[axes_names].data = image
                if image.dtype == np.dtype('uint16') and layer_image.contrast_limits[-1] <= 255:
                    VIEWER.layers[axes_names].reset_contrast_limits()
                if image.dtype == np.dtype('uint8') and layer_image.contrast_limits[-1] > 255:
                    VIEWER.layers[axes_names].reset_contrast_limits()


async def read_input_files(input_folder,
                           output_dir,
                           shared_queues_of_z_projections,
                           json_config,
                           manager,
                           factor,
                           axes,
                           exit_gracefully: threading.Event):
    try:
        awatch_input_folder = awatch(input_folder)
        unfinished_run_check = False
        while True:
            try:
                if exit_gracefully.is_set():
                    break
                changes = await asyncio.wait_for(awatch_input_folder.__anext__(), timeout=1.0)
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
                                                        factor,
                                                        axes)
                        except Exception as e:
                            print(f"Error processing file {path}: {e}")
            except asyncio.TimeoutError:
                # check the presence of unfinished processing
                if unfinished_run_check:
                    continue
                content = os.listdir(input_folder)
                if not len(content):
                    # proceed with normal run
                    unfinished_run_check = True
                    continue
                unfinished_run_check = True
                files = [f for f in content if os.path.isfile(os.path.join(input_folder, f))]
                # check json configuration in output_dir
                # first get the timelapse id from file batch
                lapse_ids = set()
                json_path = ""
                for f in files:
                    if f.upper().endswith(".TIF") or f.upper().endswith(".BMP"):
                        f_name_part = f.split("_TP")[0]
                        f_name_part = f_name_part.split("timelapseID-")[1]
                        lapse_ids.add(f_name_part)
                if len(lapse_ids) == 0:
                    debugprint(f"incorrect name pattern for existing files in {input_folder}")
                    continue
                if len(lapse_ids) == 1:
                    # default case, find and load the json config and proceed with planes
                    # first check the output folder for the config file
                    id = lapse_ids.pop()
                    # check the output directory for existence of a config file in a distinct folder
                    if os.path.exists(os.path.join(output_dir, f"AcquisitionMetadata_{id}",
                                                   f"AcquisitionMetadata_{id}.json")):
                        json_path = os.path.join(output_dir, f"AcquisitionMetadata_{id}",
                                                   f"AcquisitionMetadata_{id}.json")
                    # check the content of the output directory
                    elif os.path.exists(os.path.join(output_dir, f"AcquisitionMetadata_{id}.json")):
                        json_path = os.path.join(output_dir, f"AcquisitionMetadata_{id}.json")
                    # check the input folder for the config file presence
                    elif os.path.exists(os.path.join(input_folder, f"AcquisitionMetadata_{id}.json")):
                        json_path = os.path.join(input_folder, f"AcquisitionMetadata_{id}.json")
                    else:
                        debugprint(f"output and input directories do not contain the timelapse configuration file,"
                                   f" planes from previous run won't be processed")
                        continue
                else:
                    debugprint(f"Input folder {input_folder} contains image planes from multiple run, "
                               f"resume the processing is not supported")
                    continue
                files = [os.path.join(input_folder, f) for f in files]
                files.sort()
                files.insert(0, json_path)
                for path in files:
                    try:
                        load_file_from_input_folder(path,
                                                    output_dir,
                                                    shared_queues_of_z_projections,
                                                    json_config,
                                                    manager,
                                                    factor,
                                                    axes)
                    except Exception as e:
                        debugprint(f"Error processing file {path}: {e}")
                # Timeout is expected; just check the stop event again
                continue
            except StopAsyncIteration:
                awatch_input_folder = awatch(input_folder)
                continue
    except Exception as err:
        exit_gracefully.set()
        debugprint(err)


def watchfiles_thread(*args):
    asyncio.run(read_input_files(*args))
    debugprint("watchfiles finished")


def run_the_loop(kwargs, exit_gracefully: threading.Event):
    global PIVJLPATH
    global TEMP_DIR
    global SPECIMENS_QUANTITY_LOADED
    global SPECIMENS_QUANTITY

    input_dir = kwargs.get("input")
    output_dir = kwargs.get("output")
    axes = kwargs.get("axes", "Z")
    factor = kwargs.get("factor_anisotropy", None)
    TEMP_DIR = kwargs.get("temp_dir", None)
    PIVJLPATH = kwargs.get("pivjl", None)
    bot_config_path = kwargs.get("bot_config", "")
    specimen_quantity = kwargs.get("specimen_quantity", 0)
    process_z_projections = kwargs.get("process_z_projections", False)
    token, chat_id = "", ""

    piv_processes = []
    manager = Manager()
    json_config = manager.dict()
    shared_queues_of_z_projections = manager.dict()

    # Start the input directory observer
    debugprint(f"Watching {input_dir} for images, and saving stacks to {output_dir}")
    thread = Thread(target=watchfiles_thread, args=(input_dir,
                                                    output_dir,
                                                    shared_queues_of_z_projections,
                                                    json_config,
                                                    manager,
                                                    factor,
                                                    axes,
                                                    exit_gracefully),
                    daemon=False
                    )
    thread.start()

    while True:
        migration_detected = manager.Queue()
        avg_speed_data = manager.list()
        stop_process = Event()

        if PIVJLPATH:
            while not SPECIMENS_QUANTITY_LOADED and not process_z_projections and not exit_gracefully.is_set():
                time.sleep(1)

            if specimen_quantity:
                # update global variable to create plotting windows
                SPECIMENS_QUANTITY = (i + 1 for i in range(0, specimen_quantity))
                for i in SPECIMENS_QUANTITY:
                    shared_queues_of_z_projections[i] = manager.Queue()

            if len(SPECIMENS_QUANTITY):
                SPECIMENS_QUANTITY_LOADED = False

            if bot_config_path:
                with open(bot_config_path) as fin:
                    try:
                        line = fin.readline()
                        token, chat_id = line.strip().split(" ")
                    except Exception as err:
                        debugprint(err)

            for i, _ in shared_queues_of_z_projections.items():
                piv_process = PivProcess(
                    target=run_piv_process, args=(shared_queues_of_z_projections,
                                                  PIV_OUTPUT,
                                                  avg_speed_data,
                                                  stop_process,
                                                  PIVJLPATH,
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
                debugprint(f"Considering directory as Z projections source: {input_dir}, calculate"
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
                debugprint(f"Finished with uploading projection files for quickPIV")

        try:
            stop_file = os.path.join(input_dir, STOP_FILE_NAME)
            counter = 0
            speed_list_len = []
            while not os.path.exists(stop_file) and not SPECIMENS_QUANTITY_LOADED and not exit_gracefully.is_set():
                if migration_detected.qsize() > 0:
                    migration_event = migration_detected.get()
                    message = f"Detected on {migration_event}!"
                    if token and chat_id:
                        url = (f"https://api.telegram.org/bot{token}/"
                               f"sendMessage?chat_id={chat_id}&text={message}&disable_web_page_preview=true")
                        response = requests.get(url)
                        if response.status_code != 200:
                            debugprint(response.text)
                            debugprint(url)
                    debugprint(message)
                if process_z_projections:
                    # counter is used to check the list size approximately every minute
                    # if the avg_speed_data doesn't change - save the data and finish the pipeline
                    counter += 1
                    if counter % 30 == 0:
                        speed_list_len.append(len(avg_speed_data))
                        if len(speed_list_len) > 2:
                            speed_list_len.pop(0)
                            if speed_list_len[-1] == speed_list_len[-2]:
                                debugprint(f"No changes in quickPIV output queue, stop pipeline and save the data")
                                break
                # Sleep to keep the script running
                time.sleep(1)
            if PIVJLPATH:
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
                    debugprint(f"csv data saved to {output_dir}")
                except IOError:
                    debugprint(f"Attempt to save csv data to {output_dir} failed")
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
                        debugprint(err)
        if os.path.exists(stop_file) or process_z_projections or exit_gracefully.is_set():
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
                        image_layer_dict[f"{key}_ILL_{wrapped_dict.illumination}"] = val
                DRAWN_PROJECTIONS_QUEUE[identifier] = list(projections_dict_list)
        else:
            if projections_dict_list:
                # iterate through the list, assign illumination index to the axis, save all to the new ,
                # data remain unmerged - apply merge and push to napari
                for wrapped_dict in projections_dict_list:
                    channel_layer = draw_napari_layer(wrapped_dict.projections)
                    for key, val in channel_layer.items():
                        image_layer_dict[f"{key}_ILL_{wrapped_dict.illumination}"] = val
                DRAWN_PROJECTIONS_QUEUE[identifier] = list(projections_dict_list)
            else:
                debugprint("Empty projections dictionary")
        if image_layer_dict:
            return image_layer_dict


def bx_trigger():
    global MERGE_LIGHT_MODES
    MERGE_LIGHT_MODES = not MERGE_LIGHT_MODES


def make_napari_viewer():
    napari_viewer = napari.Viewer()
    bx = QCheckBox('Merge illumination channels')
    bx.setChecked(MERGE_LIGHT_MODES)
    bx.stateChanged.connect(bx_trigger)
    napari_viewer.window.add_dock_widget(bx)
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
    debugprint(f"quickPIV process started, PID: {os.getpid()}")
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
                x.append(t)
                avg_speed_data.append({"time_point": t, "specimen": queue_number + 1, "avg_speed": avg_speed})
                try:
                    y.append(avg_speed)
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
                                migration_detected.put(f"{frame}, specimen {queue_number + 1}")
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
            debugprint(f"window with index {index} doesn't exist")


def update_plotting_windows(exit_gracefully: threading.Event):
    global PLOTTING_WINDOW_CREATED
    global PLT_WIDGETS_DICT
    global plotting_windows_timer

    if exit_gracefully.is_set():
        if len(PLT_WIDGETS_DICT) > 0:
            for _, window in PLT_WIDGETS_DICT.items():
                window.close()
            return
            plotting_windows_timer.stop()

    if PLOTTING_WINDOW_CREATED:
        return
    if PIVJLPATH:
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
            debugprint("Recieved terminate command from microscope.")
            exit_gracefully.set()
        elif message.get("type") == "heartbeat":
            debugprint("Received heartbeat from the microscope.")
    except:
        pass


def heartbeat_and_command_handler(port, exit_gracefully: threading.Event):
    debugprint(f"heartbeat process started")
    context = zmq.Context()
    socket = context.socket(zmq.PAIR)
    socket.connect(f"tcp://localhost:{port}")
    socket.setsockopt(zmq.LINGER, 0)  # Set linger to 0 to prevent blocking on close

    poller = zmq.Poller()
    poller.register(socket, zmq.POLLIN)

    last_heartbeat_sent = time.time() - HEARTBEAT_INTERVAL_SEC
    last_heartbeat_recieved = time.time()  # initialize value with current time
    incomming_heartbeats_timeout = 60  # timeout in seconds

    try:
        while not exit_gracefully.is_set():
            next_heartbeat_time = last_heartbeat_sent + HEARTBEAT_INTERVAL_SEC
            time_until_heartbeat = max(0, next_heartbeat_time - time.time())
            # Timeout in milliseconds, max 50ms
            timeout = min(50, time_until_heartbeat * 1000)

            events = dict(poller.poll(timeout))
            if socket in events:
                message = socket.recv_json()
                debugprint("Received message: " + str(message))
                parse_message_from_microscope(message, exit_gracefully)
                if message.get("type") == "heartbeat":
                    last_heartbeat_recieved = time.time()

            if time.time() - last_heartbeat_recieved > incomming_heartbeats_timeout:
                exit_gracefully.set()

            current_time = time.time()
            if current_time >= next_heartbeat_time:
                debugprint("Sending heartbeat to the microscope.")
                socket.send_json({"type": "heartbeat", "message": "alive"})
                last_heartbeat_sent = current_time

            while not command_queue_to_microscope.empty():
                command = command_queue_to_microscope.get_nowait()
                socket.send_json({"type": "command", "message": command})

    except Exception as e:
        debugprint(f"Heartbeat handler exception: {e}")
        exit_gracefully.set()
    finally:
        debugprint("disconnecting communication socket")
        socket.close()
        poller.unregister(socket)
        context.term()
        debugprint("sending heartbeat stopped")


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
    # Specimen quantity
    parser.add_argument('--specimen_quantity', type=int,
                        required=True if '-process_z_projections' in sys.argv else False,
                        help="Value for number of specimens in the stage, defines number of plotting windows "
                             "if json config is not available")
    # Move image files with incorrect name to the user provided directory
    parser.add_argument('--temp_dir', required=False,
                        default=None,
                        help="Directory path to store input images with incorrect file name")
    parser.add_argument('--pivjl', required=True if '-process_z_projections' in sys.argv else False,
                        default=False,
                        help="Path to auxiliary julia script for average migration speed calculation with quickPIV")
    parser.add_argument('--bot_config', required=False, default=False,
                        help="Path to text file with token and chat_id information (one line, space separator) for"
                             " supplying peaks detection on telegram messenger")
    parser.add_argument('-debug_run', required=False, default=False, action='store_true',
                        help="All content of --output specified folder will be DELETED!")
    parser.add_argument('-process_z_projections',
                        required=True if '--specimen_quantity' in sys.argv else False, action='store_true',
                        help="Process input files as Z-projections for average speed calculations with quickPIV!")

    # event to stop the script gracefully
    exit_gracefully = threading.Event()
    # Parse the command-line arguments
    corrected_args = correct_argv(sys.argv)
    args = parser.parse_args(corrected_args[1:])
    if args.temp_dir is None:
        args.temp_dir = args.output
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
                debugprint(e)
    global VIEWER, PLT_WIDGETS_DICT, plotting_windows_timer

    # Allow other computers to attach to debugpy at this IP address and port.
    # debugpy.listen(('0.0.0.0', 5680))
    # debugprint("Waiting for debugger attach")
    # debugpy.wait_for_client()  # Blocks execution until client is attached

    try:
        heartbeat_thread = threading.Thread(
            target=heartbeat_and_command_handler, args=(args.port, exit_gracefully))
        heartbeat_thread.start()
    except Exception as e:
        debugprint(f"{args},\n {e}")
    thread = Thread(target=run_the_loop, args=(vars(args), exit_gracefully))
    thread.start()
    VIEWER = make_napari_viewer()
    timer1 = QTimer()
    timer1.timeout.connect(update_avg_speed_plot_windows)
    timer1.start(25)
    timer2 = QTimer()
    timer2.timeout.connect(update_napari_viewer_layer)
    timer2.start(50)
    plotting_windows_timer.timeout.connect(lambda: update_plotting_windows(exit_gracefully))
    plotting_windows_timer.start(1000)
    timer4 = QTimer()
    timer4.timeout.connect(lambda: close_napari_viewer(exit_gracefully))
    timer4.start(1000)
    debugprint(f"stack_gatherer started, PID: {os.getpid()}")
    napari.run()
    debugprint(f"Napari viewer windows was closed, terminating child processes")
    exit_gracefully.set()
    thread.join()
    heartbeat_thread.join()
    debugprint(f"stack_gatherer stopped")


if __name__ == "__main__":
    main()
