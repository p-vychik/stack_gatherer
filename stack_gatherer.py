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
from numpy.typing import NDArray
import napari
import multiprocessing
import requests
import threading
import pandas as pd
import matplotlib.pyplot as plt
import logging
from tifffile import imread, imwrite, memmap
from tqdm import tqdm
from PIL import Image
import asyncio
from watchfiles import awatch
from threading import Thread
import queue
from queue import Queue
from collections import OrderedDict
from qtpy.QtWidgets import QCheckBox
from dataclasses import dataclass, field
from typing import List
from typing import Union
from datetime import datetime
from multiprocessing import Manager
import zmq
from PlottingWindow import Ui_MainWindow
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QMainWindow
from scipy.signal import find_peaks
# drift correction import
from dexp.processing.registration.translation_nd import register_translation_nd
from scipy.ndimage import shift as scipy_shift

CURRENT_ACQUISITION_LOGGER = None
STOP_FILE_NAME = "STOP_STACK_GATHERING"
ACQUISITION_META_FILE_PATTERN = "AcquisitionMetadata_"
# BORDER_WIDTH defines the size of the border in pixels between
# projections merged in one image for napari visualization
BORDER_WIDTH = 20
PLT_WIDGET_X_AXIS_LENGTH = 40
MERGE_LIGHT_MODES = False
# separate matplotlib windows for each specimen
PLT_WIDGETS_DICT = {}
plotting_windows_timer = QTimer()
HEARTBEAT_INTERVAL_SEC = 5
command_queue_to_microscope = Queue()
new_microscope_command_event = threading.Event()
HEARTBEAT_FROM_MICROSCOPE_TIMEOUT_sec = 10
ACQUISITION_METADATA_FILE_TO_PROCESS = ""


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
        try:
            # Plot a curve:
            self.plotWidget.canvas.axes.clear()
            y_lim_range = (math.floor(min(y) * 0.9), math.ceil(max(y) * 1.1))
            self.plotWidget.canvas.axes.set_ylim(y_lim_range)
            if len(x) != len(y):
                # HAVE TO FIX DIFFERENT SIZE OF LISTS!!!
                x = x[1::]
            self.plotWidget.canvas.axes.plot(x, y, color='blue')
            if x_marker and y_marker:
                self.plotWidget.canvas.axes.plot(x_marker, y_marker,
                                                 color='red', marker='o',
                                                 markersize=12, linewidth=0)
            self.plotWidget.canvas.axes.set_xticks(x)
            self.plotWidget.canvas.axes.set_xticklabels(x, fontsize=10)
            self.plotWidget.canvas.draw()
        except Exception as err:
            logging_broadcast(f"plot_curve exception: {err}")


class TimelapseSharedState:
    def __init__(self, manager, args):
        self.args = args
        self.crop_mask_coordinates = manager.dict()
        self.json_config = manager.dict()
        self.maxprojection_queues_dict = manager.dict() if args.pivjl else None  # we use it only for quickPIV
        self.migration_detected = manager.Queue()
        self.avg_speed_data = manager.list()
        self.active_stacks = {}
        self.currently_saving_stacks_dict = manager.dict()
        self.currently_saving_stacks_lock = multiprocessing.Lock()
        self.exit_gracefully = multiprocessing.Event()
        self.napari_layers_to_show_queue = manager.Queue()
        self.quick_piv_output_queue = manager.Queue()
        self.projections_to_napari_dict = OrderedDict()
        self.shown_napari_layers_dict = OrderedDict()
        self.specimen_indices_loaded = False
        self.specimens_indices_list = manager.list()
        self.plotting_window_created = False
        self.temp_dir = ""
        self.pivjl_script_path = ""
        self.restart_timelapse_json = self.args.restart if self.args.restart else ""

    def update_piv_queues(self, manager):
        if len(self.maxprojection_queues_dict) == 0:
            for i in self.specimens_indices_list:
                self.maxprojection_queues_dict[i] = manager.Queue()
                self.crop_mask_coordinates[i] = manager.list()
            logging_broadcast(f"quickPIV queue was updated")

    def reset_shared_queues(self, specimen_indices, manager):
        self.specimens_indices_list.extend(sorted(specimen_indices))
        for index in specimen_indices:
            self.maxprojection_queues_dict[index] = manager.Queue()

    def print_attributes_content_size(self):
        # function to test class attributes size in case of excessive memory use
        for attr_name, attr_value in self.__dict__.items():
            # Check if the attribute is a manager.Queue (using the general Queue type)
            if isinstance(attr_value, multiprocessing.queues.Queue):
                print(f"Attribute: {attr_name}, Type: Queue (proxy), Size: Cannot determine size")

            # Check if the attribute is a manager.list (using ListProxy for managed lists)
            elif isinstance(attr_value, list) or hasattr(attr_value, '__len__'):  # manager.list proxies have len()
                print(f"Attribute: {attr_name}, Type: List (proxy), Size: {len(attr_value)}")

            # Special handling for maxprojection_queues_dict
            elif attr_name == "maxprojection_queues_dict" and isinstance(attr_value, dict):
                print(f"Attribute: {attr_name}, Type: Dict of Queues (proxy), Number of Queues: {len(attr_value)}")
                for key, queue in attr_value.items():
                    if isinstance(queue, multiprocessing.queues.Queue):
                        print(f"  Queue for Key: {key}, Size: Cannot determine size")


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

    @property
    def time_point(self) -> int:
        return self.stack_signature_obj.time_point


@dataclass
class JsonConfigFile:
    _processed: bool
    _path: str
    _lapse_id: str

    def set_processed(self):
        self._processed = True

    @property
    def path(self):
        return self._path

    @property
    def is_processed(self):
        return self._processed

    @property
    def get_lapse_id(self):
        return self._lapse_id


@dataclass
class PivTimeSeries:
    x: List[float] = field(default_factory=list)
    y: List[float] = field(default_factory=list)
    t: int = 0
    peak: List[bool] = field(default_factory=list)


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


def read_image(image_path: str) -> Union[NDArray, bool]:
    if image_path.upper().endswith((".TIF", ".TIFF")):
        time.sleep(0.1)
        try:
            image = imread(image_path)
            return image
        except Exception as error:
            logging_broadcast(f"read_image {error}")
    if image_path.upper().endswith(".BMP"):
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
            logging_broadcast(f"read_image {error}")
    return False


def file_name_merged_illumination_based_on_signature(stack_signature):
    return (f"timelapseID-{stack_signature.timelapse_id}_"
            f"SPC-{stack_signature.specimen}_"
            f"TP-{stack_signature.time_point}_"
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


def collect_files_to_one_stack_get_axial_projections(shared_state: TimelapseSharedState,
                                                     stack_signature: StackSignature,
                                                     file_list: list,
                                                     output_file_path: str,
                                                     axes=None,
                                                     anisotropy_factor=1,
                                                     output_dir=None):
    sample_image = read_image(file_list[0])
    if isinstance(sample_image, bool):
        return
    shape = (len(file_list), sample_image.shape[0], sample_image.shape[1])
    dtype = sample_image.dtype
    # create an empty OME-TIFF file
    if not os.path.exists(output_file_path):
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
        except Exception as err:
            logging_broadcast(f" collect_files_to_one_stack_get_axial_projections {err}")
    # write data to memory-mapped array
    with tqdm(total=len(file_list), desc="Saving plane") as pbar:
        for z in range(shape[0]):
            if z == 0:
                # check if projection's file exists in output folder otherwise it may cause problem with broadcast data
                zyx_stack[z] = sample_image
                if shape[0] == 1:
                    # stack which consists of one plane means that on input are Z projections,
                    # clear dictionary to skip processing
                    projections.clear()
                    projections["Z"] = read_image(file_list[z])
                    pbar.update(1)
                    break
                pbar.update(1)
                continue
            zyx_stack[z] = read_image(file_list[z])
            projections = plane_to_projection(zyx_stack[z], projections)
            # projections = plane_to_projection(zyx_stack[z], projections)
            pbar.update(1)
        zyx_stack.flush()
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
            if shared_state.maxprojection_queues_dict is not None:
                if stack_signature.specimen in shared_state.maxprojection_queues_dict:
                    shared_state.maxprojection_queues_dict[stack_signature.specimen].put(wrapped_z_projection)
                else:
                    logging_broadcast(f"Specimen signature not found in shared queue, signature: {stack_signature}")
        for axis in projections.keys():
            try:
                if not os.path.exists(projections_files_path[axis]):
                    imwrite(projections_files_path[axis], projections[axis])
            except Exception as err:
                logging_broadcast(f"Can't save the projection file, {err}")

        # we want to store projections derived from planes with different illumination modes to
        # have an option to merge them if the user checked the QCheckBox() in napari viewer GUI interface
        wrapped_projections = ProjectionsDictWrapper(projections, stack_signature)
        if wrapped_projections.signature not in shared_state.projections_to_napari_dict:
            shared_state.projections_to_napari_dict[wrapped_projections.identifier] = [wrapped_projections]
        else:
            shared_state.projections_to_napari_dict[wrapped_projections.identifier].append(wrapped_projections)

    for z in range(shape[0]):
        try:
            os.remove(file_list[z])
        except PermissionError as err:
            logging_broadcast(f"Can't remove the plane file: {err}")


def add_file_to_active_stacks(shared_state: TimelapseSharedState,
                              image_file: ImageFile):
    stack_signature = image_file.get_stack_signature()
    if stack_signature.signature not in shared_state.active_stacks:
        shared_state.active_stacks[stack_signature.signature] = {}
        logging_broadcast(f"Adding stack {stack_signature.signature} to active queue.")
    if image_file.plane not in shared_state.active_stacks[stack_signature.signature]:
        shared_state.active_stacks[stack_signature.signature][image_file.plane] = image_file
        return stack_signature
    else:
        return False


def check_stack_and_collect_if_ready(shared_state: TimelapseSharedState,
                                     stack_signature: StackSignature,
                                     output_dir,
                                     axes=None,
                                     factor=None):
    if len(shared_state.active_stacks[stack_signature.signature]) < stack_signature.total_num_planes:
        return
    # We have to ensure that two events firing at the same time don't start saving the same stack twice
    with shared_state.currently_saving_stacks_lock:
        if stack_signature not in shared_state.currently_saving_stacks_dict:
            shared_state.currently_saving_stacks_dict[stack_signature.signature] = True
        else:
            return
        file_list = []
        for i, _ in enumerate(shared_state.active_stacks[stack_signature.signature]):
            # We have to access by index since we can't guarantee that files were added to dict in order of planes
            file_list.append(shared_state.active_stacks[stack_signature.signature][i].get_file_path())
        sample_file_obj = ImageFile(file_list[0])
        sample_file_obj.extension = "tif"
        stack_path = os.path.join(output_dir, sample_file_obj.get_stack_name())
        collect_files_to_one_stack_get_axial_projections(shared_state,
                                                         stack_signature,
                                                         file_list,
                                                         stack_path,
                                                         axes=axes,
                                                         anisotropy_factor=factor,
                                                         output_dir=output_dir
                                                         )

        images_dict = get_projections_dict_from_queue(shared_state)
        if images_dict:
            shared_state.napari_layers_to_show_queue.put(images_dict)
        del shared_state.active_stacks[stack_signature.signature]
        del shared_state.currently_saving_stacks_dict[stack_signature.signature]


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
        logging_broadcast(f"{file_path} failed to parse JSON {e}")
        return False

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


def load_file_from_input_folder(shared_state: TimelapseSharedState,
                                manager,
                                file_path,
                                output_dir,
                                factor=1,
                                axes=""):
    if os.path.isdir(file_path):
        return
    # load json parameters
    if os.path.basename(file_path).startswith(ACQUISITION_META_FILE_PATTERN) and file_path.endswith(".json"):
        logging_broadcast(f"try to parse {file_path}")
        destination_folder = os.path.join(output_dir, os.path.basename(file_path).strip(".json"))
        setup_acquisition_log(destination_folder)
        logging_broadcast(f"new logger is set, prepare to parse JSON {file_path}")
        json_parsing_result = load_lapse_parameters_json(file_path, destination_folder)
        if not json_parsing_result:
            logging_broadcast(f"{file_path} parsing as JSON file failed")
            return
        list_of_stack_signatures, metadata_json = json_parsing_result
        logging_broadcast(f"update lapse config from {file_path}")
        shared_state.json_config.update(dict.fromkeys(list_of_stack_signatures, metadata_json))
        # get specimens quantity
        try:
            shared_state.specimens_indices_list = [specimen_dict['userDefinedIndex'] for specimen_dict in
                                                   metadata_json['specimens']]
            shared_state.specimen_indices_loaded = True
        except Exception as err:
            logging_broadcast(f"Failed to set specimen quantity {err}")
        if len(shared_state.specimens_indices_list) and shared_state.maxprojection_queues_dict is not None:
            shared_state.update_piv_queues(manager)

        shared_state.plotting_window_created = False

        try:
            if not os.path.exists(os.path.join(destination_folder, os.path.basename(file_path))):
                try:
                    shutil.move(file_path, destination_folder)
                except Exception as err:
                    logging_broadcast(f"Attempt to move {file_path} reulted in {err}")
            else:
                logging_broadcast(f"file {os.path.basename(file_path)} exists in"
                                  f"  {os.path.join(destination_folder, os.path.basename(file_path))}, "
                                  f"so it would be removed")
                try:
                    os.remove(file_path)
                except Exception as err:
                    logging_broadcast(f"{err}")
        except Exception as err:
            logging_broadcast(f"load_file_from_input_folder {err}")
    if file_path.upper().endswith((".TIFF", ".BMP", ".TIF")):
        # Call the function when a new file is created
        try:
            file = ImageFile(file_path)
        except NotImagePlaneFile:
            # image file is incorrect, move it to the temp_dir
            try:
                shutil.move(file_path, shared_state.temp_dir)
            except Exception as err:
                logging_broadcast(f"Moving {file_path} to {shared_state.temp_dir} resulted in {err}")
            return
        stack_signature = add_file_to_active_stacks(shared_state, file)
        # create separate folder for the output based on metadata filename
        try:
            if not isinstance(stack_signature, bool):
                output_dir = shared_state.json_config[stack_signature.signature_no_time]["output_folder"]
                if os.path.exists(output_dir):
                    check_stack_and_collect_if_ready(shared_state,
                                                     stack_signature,
                                                     output_dir,
                                                     axes,
                                                     factor)
            else:
                return
        except Exception as err:
            logging_broadcast(f"load_file_from_input_folder was unsuccessful, {err}")


def extract_and_display_selection(shape_layer, min_coords, max_coords, shared_state: TimelapseSharedState):
    try:
        print(shape_layer.name)
        image_layer_name = shape_layer.name.split("Crop mask ")[1]
        # image_layer_name = re.sub(r"TP_\d+", "", image_layer_name)
        specimen = int(image_layer_name.split("SPC_")[1].split("_")[0])
        if specimen in shared_state.crop_mask_coordinates:
            shared_state.crop_mask_coordinates[specimen] = (min_coords, max_coords)
        else:
            logging_broadcast("Can't save cropping mask into the shared dictionary, key index is not present!")
    except Exception as err:
        logging_broadcast(f"Exception extract_and_display: {err}")


def on_rectangle_selection(shapes_layer, event, shared_state: TimelapseSharedState):
    # Check if the current tool is the rectangle tool
    # Start with the initial mouse press event
    if event.type == "mouse_press" and shapes_layer.mode == "add_rectangle":
        # Clear any existing rectangles to ensure only one is visible
        shapes_layer.data = []
    # Start with the initial press
    yield
    while event.type == "mouse_move":
        # Handling the drag/move event
        yield
    if len(shapes_layer.data) > 0 and event.type == "mouse_release" and shapes_layer.mode == "add_rectangle":
        # Get the last drawn rectangle
        rectangle = shapes_layer.data[-1]
        min_coords = np.min(rectangle, axis=0).astype(int)
        max_coords = np.max(rectangle, axis=0).astype(int)
        extract_and_display_selection(shapes_layer, min_coords, max_coords, shared_state)


def update_napari_viewer_layer(shared_state: TimelapseSharedState, napari_viewer) -> None:
    if shared_state.napari_layers_to_show_queue.qsize() > 0:
        data_input = shared_state.napari_layers_to_show_queue.get()
        for axes_names, image in data_input.items():
            match = re.search(r'TP_(\d+)', axes_names)
            if match:
                time_point = int(match.group(1))
                drawn_layer_key = re.sub(r'TP_\d+', f'TP_{time_point - 1}', axes_names)

                if drawn_layer_key not in napari_viewer.layers:
                    # Add a new layer if it doesn't exist
                    layer_image = napari_viewer.add_image(image, name=axes_names)
                    mask_name = re.sub(r'TP_\d+_', "", axes_names)
                    shape_layer = napari_viewer.add_shapes(name=f"Crop mask {mask_name}", shape_type='rectangle')
                    # track changes on shape_layer
                    shape_layer.mouse_drag_callbacks.append(lambda layer, event:
                                                            on_rectangle_selection(layer,
                                                                                   event,
                                                                                   shared_state)
                                                            )

                else:
                    # Retrieve the existing layer
                    layer_image = napari_viewer.layers[drawn_layer_key]

                    # Check if there is an existing selection
                    existing_selection = None
                    if hasattr(layer_image, 'selected_data') and len(layer_image.selected_data) > 0:
                        existing_selection = layer_image.selected_data[0]

                    # Update the layer data and rename it
                    layer_image.data = image
                    layer_image.name = axes_names

                    # Restore the previous selection if it exists
                    if existing_selection:
                        layer_image.selected_data = [existing_selection]
                    # layer_image.mouse_drag_callbacks.append(callback_wrapper)
                # Adjust contrast limits if necessary
                if image.dtype == np.dtype('uint16') and layer_image.contrast_limits[-1] <= 255:
                    napari_viewer.layers[axes_names].reset_contrast_limits()
                if image.dtype == np.dtype('uint8') and layer_image.contrast_limits[-1] > 255:
                    napari_viewer.layers[axes_names].reset_contrast_limits()
        visible_specimens = set()
        for layer in napari_viewer.layers:
            specimen_index = parse_specimen_index(layer.name)
            if layer in napari_viewer.layers.selection or layer.visible:
                visible_specimens.add(specimen_index)
            if specimen_index not in visible_specimens:
                layer.visible = False
        if not len(napari_viewer.layers.selection) and len(napari_viewer.layers) >= 2:
            napari_viewer.layers[-1].visible = True
            napari_viewer.layers[-2].visible = True


def get_lapse_id(file_name):
    try:
        base_name = os.path.basename(file_name)
        if base_name.upper().endswith('.JSON'):
            return base_name.split('_')[1].strip('.json')
        elif base_name.upper().endswith(('.BMP', '.TIFF', '.TIF')):
            return base_name.split("_SPC")[0].split("timelapseID-")[1]
        else:
            logging_broadcast(f"get_lapse_id failed to parse {file_name} for timelapseID")
    except IndexError:
        logging_broadcast(f"IndexError: get_lapse_id failed to parse {file_name} for timelapseID")
        return ""


def find_config_files_locations(config_files_locations: dict,
                                output_folder: str,
                                input_folder: str,
                                json_file_name=None):
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

    Returns:
    - dict: A dictionary mapping lapse IDs to the full paths of the JSON files found.
            If json_file_name is provided, returns a dictionary with a single entry.
    """

    def collect_jsons(file_name):

        base_name = os.path.splitext(os.path.basename(file_name))[0]
        lapse_id = get_lapse_id(file_name)
        if not file_name.endswith('.json'):
            # we deal with image file and should construct the JSON file name based on pattern
            base_name = f"{ACQUISITION_META_FILE_PATTERN}{lapse_id}"
            file_name = f"{ACQUISITION_META_FILE_PATTERN}{lapse_id}.json"

        logging_broadcast(f"start collecting jsons")
        if lapse_id:
            if os.path.exists(os.path.join(output_folder, file_name)):
                config_files_locations[lapse_id] = JsonConfigFile(False,
                                                                  os.path.join(output_folder, file_name),
                                                                  lapse_id
                                                                  )
                logging_broadcast(f"found {config_files_locations[lapse_id].path}")
            elif base_name and os.path.exists(os.path.join(output_folder, base_name, file_name)):
                config_files_locations[lapse_id] = JsonConfigFile(False,
                                                                  os.path.join(output_folder, base_name, file_name),
                                                                  lapse_id
                                                                  )
                logging_broadcast(f"found {config_files_locations[lapse_id].path}")
            elif os.path.exists(os.path.join(input_folder, file_name)):
                config_files_locations[lapse_id] = JsonConfigFile(False,
                                                                  os.path.join(input_folder, file_name),
                                                                  lapse_id
                                                                  )
                logging_broadcast(f"found {config_files_locations[lapse_id].path}")
            else:
                logging_broadcast(f"checked all location unsuccessfully: {os.path.join(output_folder, file_name)} "
                                  f"{os.path.join(output_folder, base_name, file_name)} "
                                  f"{os.path.join(input_folder, file_name)}")

    # If json_file_name is provided, directly check its presence
    if json_file_name:
        collect_jsons(json_file_name)
        return

    # If json_file_name is not provided, collect lapse ID from planes present in the folder and continue
    # with search in input and output folders
    content = list(os.listdir(input_folder))
    if len(content) == 0:
        logging_broadcast(f"input directory doesn't contain images or config files")
        return

    lapse_ids = set()
    file_paths_with_unique_lapse_ids = []

    for file_name in content:
        lapse_id = get_lapse_id(file_name)
        if lapse_id:
            if lapse_id in lapse_ids:
                continue
            lapse_ids.add(lapse_id)
            file_paths_with_unique_lapse_ids.append(file_name)

    if file_paths_with_unique_lapse_ids:
        for file_name in file_paths_with_unique_lapse_ids:
            collect_jsons(file_name)


async def read_input_files(shared_state,
                           manager,
                           input_folder,
                           output_dir,
                           factor,
                           axes,
                           input_files_queue):
    try:
        # Create an asynchronous generator from awatch
        async_gen = awatch(input_folder)
        while True:
            try:
                # Set a short timeout for the next value from the async generator
                changes = await asyncio.wait_for(async_gen.__anext__(), timeout=1)
                # Process changes if any
                paths = []
                for change in changes:
                    if shared_state.exit_gracefully.is_set():
                        break
                    event, file_path = change
                    if event.value == 1:  # File created
                        paths.append(file_path)

                if paths:
                    paths.sort()
                    for path in paths:
                        try:
                            input_files_queue.put(path)
                        except Exception as err:
                            logging_broadcast(f"Error processing file {path}: {err}")
            except asyncio.TimeoutError:
                # During timeout, check if the exit Event is set
                if shared_state.exit_gracefully.is_set():
                    break
            except StopAsyncIteration:
                async_gen = awatch(input_folder)
                continue

    except Exception as err:
        shared_state.exit_gracefully.set()
        logging_broadcast(f"read_input_files {err}")


def watchfiles_thread(*args):
    asyncio.run(read_input_files(*args))
    logging_broadcast("watchfiles finished")


def processing_input_files_thread(shared_state: TimelapseSharedState,
                                  file_path_queue,
                                  manager,
                                  input_dir,
                                  output_dir,
                                  factor,
                                  axes):

    next_check = time.time()
    content = list(os.scandir(input_dir))
    if shared_state.restart_timelapse_json:
        json_paths_dict = {}
        find_config_files_locations(json_paths_dict,
                                    output_dir, input_dir,
                                    shared_state.restart_timelapse_json)
        if len(json_paths_dict):
            for _, json_file in json_paths_dict.items():
                load_file_from_input_folder(shared_state, manager, json_file.path, output_dir, factor, axes)

    while not shared_state.exit_gracefully.is_set():
        try:
            if time.time() - next_check <= 0:
                content = list(os.scandir(input_dir))
                next_check = time.time() + 30
            if content:
                file_paths = [os.path.join(input_dir, f) for f in content]
                for path in file_paths:
                    load_file_from_input_folder(shared_state, manager, path, output_dir, factor, axes)
                content = []
            # Get file path from queue (block until an item is available)
            path = file_path_queue.get(timeout=1)
            load_file_from_input_folder(shared_state, manager, path, output_dir, factor, axes)
        except queue.Empty:
            continue
        except Exception as err:
            logging_broadcast(f"Error processing file from queue: {err}")


def run_the_loop(kwargs, shared_state: TimelapseSharedState, manager: Manager):
    input_dir = kwargs.get("input")
    output_dir = kwargs.get("output")
    axes = kwargs.get("axes", "Z")
    factor = kwargs.get("factor_anisotropy", None)
    shared_state.temp_dir = kwargs.get("temp_dir", None)
    shared_state.pivjl_script_path = kwargs.get("pivjl", "")
    bot_config_path = kwargs.get("bot_config", "")
    process_z_projections = kwargs.get("process_z_projections", False)
    token, chat_id = "", ""
    fix_drift = kwargs.get("correct_drift", False)
    crop_margin = kwargs.get("margin_crop", 100)

    # Start the input directory observer
    logging_broadcast(f"Watching {input_dir} for images, and saving stacks to {output_dir}")
    # queue for input_data
    plane_files_queue = Queue()

    thread = Thread(target=watchfiles_thread, args=(shared_state,
                                                    manager,
                                                    input_dir,
                                                    output_dir,
                                                    factor,
                                                    axes,
                                                    plane_files_queue),
                    daemon=True
                    )
    thread2 = Thread(target=processing_input_files_thread,
                     args=(shared_state, plane_files_queue, manager, input_dir, output_dir, factor, axes),
                     daemon=True
                     )
    if not process_z_projections:
        thread.start()
        thread2.start()

    while True:
        # RESET avg_speed_data and migration_detected for shared_state object!

        if shared_state.pivjl_script_path:
            while (not shared_state.specimen_indices_loaded
                   and not process_z_projections
                   and not shared_state.exit_gracefully.is_set()):
                time.sleep(1)
                if shared_state.exit_gracefully.is_set():
                    break
            if process_z_projections:
                specimen_indices = set()
                for f in os.listdir(input_dir):
                    if not f.upper().endswith((".TIF", ".BMP")):
                        continue
                    match = re.search(r"(?<=SPC-)\d+", f)
                    if match:
                        specimen_indices.add(int(match.group()))
                if specimen_indices:
                    shared_state.reset_shared_queues(specimen_indices, manager)

            if len(shared_state.specimens_indices_list):
                shared_state.specimen_indices_loaded = False

            if bot_config_path:
                with open(bot_config_path) as fin:
                    try:
                        line = fin.readline()
                        token, chat_id = line.strip().split(" ")
                    except Exception as err:
                        logging_broadcast(err)
            piv_process = multiprocessing.Process(
                target=run_piv_process, args=(shared_state,
                                              process_z_projections,
                                              output_dir,
                                              fix_drift,
                                              crop_margin
                                              ),
                name='piv_run', daemon=True)
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
                    if isinstance(img_data, bool):
                        continue
                    img_metadata = ImageFile(file_path)
                    if not isinstance(img_metadata, bool):
                        stack_signature = img_metadata.get_stack_signature()
                        wrapped_z_projection = ProjectionsDictWrapper({"Z": img_data}, stack_signature)
                        # push projection to the queue for PIV calculations
                        if stack_signature.specimen in shared_state.maxprojection_queues_dict:
                            while True:
                                if shared_state.maxprojection_queues_dict[stack_signature.specimen].qsize() > 10:
                                    time.sleep(1)
                                else:
                                    break
                            if wrapped_z_projection.signature not in shared_state.projections_to_napari_dict:
                                shared_state.projections_to_napari_dict[wrapped_z_projection.identifier] = \
                                    [wrapped_z_projection]
                            else:
                                shared_state.projections_to_napari_dict[wrapped_z_projection.identifier].append(
                                    wrapped_z_projection)
                            images_dict = get_projections_dict_from_queue(shared_state)
                            if images_dict:
                                shared_state.napari_layers_to_show_queue.put(images_dict)
                            shared_state.maxprojection_queues_dict[stack_signature.specimen].put(wrapped_z_projection)
                        else:
                            logging_broadcast(f"Specimen index {stack_signature.specimen} was not found, check "
                                              f"if file name format follows the expected name pattern")
                logging_broadcast(f"Finished with uploading projection files for quickPIV")

        try:
            stop_file = os.path.join(input_dir, STOP_FILE_NAME)
            check_next = time.time() + 60
            speed_list_length = len(shared_state.avg_speed_data)
            while (not os.path.exists(stop_file)
                   and not shared_state.specimen_indices_loaded
                   and not shared_state.exit_gracefully.is_set()):
                if shared_state.migration_detected.qsize() > 0:
                    migration_event = shared_state.migration_detected.get()
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
                    # counter is used to check the list size approximately every minute
                    # if the avg_speed_data doesn't change - save the data and finish the pipeline
                    if check_next - time.time() <= 0:
                        if speed_list_length == len(shared_state.avg_speed_data):
                            logging_broadcast(f"No changes in quickPIV output queue, "
                                              f"stop pipeline and save the data")
                            break
                        check_next = time.time() + 60
                        speed_list_length = len(shared_state.avg_speed_data)
                # for tests only
                # if check_next - time.time() <= 0:
                #     shared_state.print_attributes_content_size()
                #     check_next = time.time() + 60
                # Sleep to keep the script running
                time.sleep(1)
            if shared_state.pivjl_script_path:
                shared_state.exit_gracefully.set()
        except KeyboardInterrupt:
            # Gracefully stop the observer if the script is interrupted
            shared_state.exit_gracefully.set()

        # save PIV data to csv
        if shared_state.avg_speed_data:
            sorted_speed_data = sorted(shared_state.avg_speed_data, key=lambda d: d['specimen'])
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
                        logging_broadcast(f"run_the_loop: {err}")
        if (os.path.exists(stop_file)
                or process_z_projections
                or shared_state.exit_gracefully.is_set()):
            break
    if shared_state.pivjl_script_path:
        piv_process.join()
    if not process_z_projections:
        thread2.join()
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


def get_projections_dict_from_queue(shared_state: TimelapseSharedState):
    if len(shared_state.projections_to_napari_dict) > 0:
        image_layer_dict = {}
        # store last 4 projections drawn in napari layer for merging channels purpose
        if len(shared_state.shown_napari_layers_dict) > 4:
            shared_state.shown_napari_layers_dict.popitem(last=False)
        identifier, projections_dict_list = shared_state.projections_to_napari_dict.popitem(last=False)
        if MERGE_LIGHT_MODES:
            # should check the number of illuminations
            if len(projections_dict_list) == 2:
                merged_projections_dict = merge_multiple_projections((identifier, projections_dict_list))
                image_layer_dict = draw_napari_layer(merged_projections_dict)
            elif identifier in shared_state.shown_napari_layers_dict and projections_dict_list:
                projections_dict_list.extend(shared_state.shown_napari_layers_dict[identifier])
                merged_projections_dict = merge_multiple_projections((identifier, projections_dict_list))
                image_layer_dict = draw_napari_layer(merged_projections_dict)
                del shared_state.shown_napari_layers_dict[identifier]
            elif projections_dict_list:
                for wrapped_dict in projections_dict_list:
                    channel_layer = draw_napari_layer(wrapped_dict.projections)
                    for key, val in channel_layer.items():
                        image_layer_dict[(f"SPC_{wrapped_dict.specimen}_TP_{wrapped_dict.time_point}_{key}"
                                          f"_ILL_{wrapped_dict.illumination}")] = val
                shared_state.shown_napari_layers_dict[identifier] = list(projections_dict_list)
        else:
            if projections_dict_list:
                # iterate through the list, assign illumination index to the axis, save all to the new ,
                # data remain unmerged - apply merge and push to napari
                for wrapped_dict in projections_dict_list:
                    channel_layer = draw_napari_layer(wrapped_dict.projections)
                    for key, val in channel_layer.items():
                        image_layer_dict[(f"SPC_{wrapped_dict.specimen}_TP_{wrapped_dict.time_point}_{key}"
                                          f"_ILL_{wrapped_dict.illumination}")] = val
                # shared_state.shown_napari_layers_dict[identifier] = list(projections_dict_list)
            else:
                logging_broadcast("Empty projections dictionary")
        if image_layer_dict:
            return image_layer_dict


def bx_trigger():
    global MERGE_LIGHT_MODES
    MERGE_LIGHT_MODES = not MERGE_LIGHT_MODES


def parse_specimen_index(name):
    """Extract the specimen index from the event value name"""
    try:
        return name.split("SPC_")[1].split("_")[0]
    except (IndexError, AttributeError) as err:
        print(f"Problem with layer name parsing, def on_selection: {err}")
        return None


def update_layer_visibility(viewer, specimen_index):
    """Update the visibility and selection of layers based on the specimen index"""
    for layer in viewer.layers:
        is_target_layer = f"SPC_{specimen_index}_" in layer.name
        layer.visible = is_target_layer


def on_selection(viewer, event):
    """Handle layer selection based on the event"""
    if not event.value:
        return

    specimen_index = parse_specimen_index(event.value.name)
    if specimen_index:
        update_layer_visibility(viewer, specimen_index)


def make_napari_viewer():
    napari_viewer = napari.Viewer()
    bx = QCheckBox('Merge illumination channels')
    bx.setChecked(MERGE_LIGHT_MODES)
    bx.stateChanged.connect(bx_trigger)
    napari_viewer.window.add_dock_widget(bx)
    napari_viewer.layers.selection.events.active.connect(lambda event: on_selection(napari_viewer, event))
    return napari_viewer


def run_piv_process(shared_state: TimelapseSharedState,
                    save_merged_illumination: bool,
                    output_dir: str,
                    fix_drift: bool,
                    crop_margin: int,
                    ):
    from juliacall import Main as jl
    jl.include(shared_state.pivjl_script_path)
    piv_projection_queue, projections_to_process, ts_dict = {}, {}, {}
    # queue_number corresponds to the specimen index from JSON config file or
    # to the SPC_ value in plane name if we process collected z-projections
    for queue_number, _ in shared_state.maxprojection_queues_dict.items():
        piv_projection_queue[queue_number] = Queue()
        projections_to_process[queue_number] = Queue()
        ts_dict[queue_number] = PivTimeSeries(x=[], y=[], t=0, peak=[])
    migration_event_frame = set()
    logging_broadcast(f"quickPIV process started, PID: {os.getpid()}")
    csv_file = os.path.join(output_dir, f"piv_avg_speed_{datetime.now().strftime('%H_%M_%S')}.csv")
    with open(csv_file, mode='a', newline='') as file:
        writer = csv.writer(file)
        if file.tell() == 0:
            writer.writerow(['time_point', 'specimen', 'avg_speed', 'peak'])
        while not shared_state.exit_gracefully.is_set():
            for queue_number, queue_in in shared_state.maxprojection_queues_dict.items():
                if shared_state.exit_gracefully.is_set():
                    break
                # logging_broadcast(f"piv heartbeat #{queue_number}")
                start_time = time.time()
                try:
                    plane_matrix = shared_state.maxprojection_queues_dict[queue_number].get(timeout=0.05)
                except Exception as err:
                    continue
                if plane_matrix:
                    projections_to_process[queue_number].put(plane_matrix)
                if projections_to_process[queue_number].qsize() >= 2:
                    wrapped_z_p_1 = projections_to_process[queue_number].get()
                    # to calculate avg speed in consecutive way we store
                    # the last projection from every pairwise comparison
                    wrapped_z_p_2 = projections_to_process[queue_number].queue[0]
                    if wrapped_z_p_1.identifier == wrapped_z_p_2.identifier:
                        merged_projections = merge_multiple_projections((wrapped_z_p_1.identifier,
                                                                         (wrapped_z_p_1,
                                                                          wrapped_z_p_2))
                                                                        )
                        merged_wrapped_proj = ProjectionsDictWrapper(merged_projections,
                                                                     wrapped_z_p_2.signature)
                        if save_merged_illumination:
                            try:
                                file_name = file_name_merged_illumination_based_on_signature(
                                    wrapped_z_p_2.stack_signature)
                                imwrite(os.path.join(output_dir, f"{file_name}.tif"),
                                        merged_wrapped_proj.projections["Z"])
                            except Exception as err:
                                shared_state.migration_detected.put(f"Failed to save the merged Z-projection, {err}")
                        piv_projection_queue[queue_number].put(merged_wrapped_proj)
                        projections_to_process[queue_number].get()
                    else:
                        piv_projection_queue[queue_number].put(wrapped_z_p_1)

                    if piv_projection_queue[queue_number].qsize() < 2:
                        continue
                    m_1 = piv_projection_queue[queue_number].get().projections["Z"]
                    # to compare projections in a consecutive way,
                    # we have to leave the last projection in the piv_projection_queue
                    m_2 = piv_projection_queue[queue_number].queue[0].projections["Z"]
                    # crop images according to user selection
                    try:
                        if shared_state.crop_mask_coordinates[queue_number]:
                            min_coords, max_coords = shared_state.crop_mask_coordinates[queue_number]
                            m_1 = m_1[min_coords[0]:max_coords[0], min_coords[1]:max_coords[1]]
                            m_2 = m_2[min_coords[0]:max_coords[0], min_coords[1]:max_coords[1]]
                            logging_broadcast(f"piv run with crop mask for specimen {queue_number},"
                                              f" coordinates min {min_coords[0], min_coords[1]}, "
                                              f"max {max_coords[0], max_coords[1]}")
                    except Exception as err:
                        logging_broadcast(f"crop mask resulted in error: {err}")
                    # correct planes drift:
                    if fix_drift:
                        logging_broadcast(f"start drift correction")
                        m_1_copy = m_1.copy()
                        m_2_copy = m_2.copy()
                        # Check if shapes are appropriate for cropping
                        if (m_1.shape[0] <= 2 * crop_margin) or (m_1.shape[1] <= 2 * crop_margin) or \
                                (m_2.shape[0] <= 2 * crop_margin) or (m_2.shape[1] <= 2 * crop_margin):
                            logging_broadcast(f"Image dimensions are too small for the specified crop margin."
                                              f" image 1 shape {m_1.shape}, image 2 shape {m_2.shape}, selected margin"
                                              f" crop value {crop_margin}")
                        else:
                            try:
                                translation_model = register_translation_nd(m_1, m_2)
                                m_1 = m_1[crop_margin:-crop_margin, crop_margin:-crop_margin]
                                corrected_m_2 = scipy_shift(m_2, translation_model.shift_vector)
                                m_2 = corrected_m_2[crop_margin:-crop_margin, crop_margin:-crop_margin]
                            except Exception as err:
                                logging_broadcast(f"drift correction wasn't successful, {err}")
                                m_1 = m_1_copy
                                m_2 = m_2_copy
                    if m_1.shape == m_2.shape:
                        try:
                            start_piv = time.time()
                            avg_speed = jl.fn(m_1, m_2)
                            logging_broadcast(f"piv run took {time.time() - start_piv}")
                            if isinstance(avg_speed[-1], float) or isinstance(avg_speed[-1], int):
                                avg_speed = round(avg_speed[-1], 3)
                            else:
                                logging_broadcast(f"quickPIV returned abnormal value {avg_speed[-1]} "
                                                  f"for specimen #{queue_number}, "
                                                  f"frame #{ts_dict[queue_number].t + 1}")
                        except Exception as error:
                            raise error
                        ts_dict[queue_number].t += 1
                        ts_dict[queue_number].x.append(ts_dict[queue_number].t)
                        # we add default value which might be changed after find_peaks step
                        ts_dict[queue_number].peak.append(False)
                        try:
                            ts_dict[queue_number].y.append(avg_speed)
                        except Exception as error:
                            logging_broadcast(f"failed to update time-series list: {error}")
                        if len(ts_dict[queue_number].x) > PLT_WIDGET_X_AXIS_LENGTH:
                            ts_dict[queue_number].x = ts_dict[queue_number].x[-PLT_WIDGET_X_AXIS_LENGTH::]
                            ts_dict[queue_number].y = ts_dict[queue_number].y[-PLT_WIDGET_X_AXIS_LENGTH::]
                        try:
                            peaks, _ = find_peaks(np.asarray(ts_dict[queue_number].y), prominence=1.5)
                        except Exception as error:
                            raise error
                        x_marker, y_marker = [], []
                        if len(peaks):
                            try:
                                peaks = peaks.tolist()
                                x_marker = [ts_dict[queue_number].x[i] for i in peaks]
                                y_marker = [ts_dict[queue_number].y[i] for i in peaks]
                                if len(ts_dict[queue_number].x) > PLT_WIDGET_X_AXIS_LENGTH:
                                    global_peaks_indices = [x + ts_dict[queue_number].t - PLT_WIDGET_X_AXIS_LENGTH
                                                            + 1 for x in x_marker]
                                else:
                                    global_peaks_indices = x_marker
                                for frame in global_peaks_indices:
                                    if frame not in migration_event_frame:
                                        migration_event_frame.add(frame)
                            except Exception as error:
                                raise error

                        writer.writerow([ts_dict[queue_number].t,
                                         queue_number,
                                         avg_speed,
                                         False])
                        file.flush()
                        shared_state.quick_piv_output_queue.put(
                            (queue_number, ts_dict[queue_number].x, ts_dict[queue_number].y, x_marker, y_marker))
                        logging_broadcast(f"plane from queue #{queue_number} was processed")
                    else:
                        raise ValueError("Projections should have the same size for QuickPIV input")
                logging_broadcast(f"for specimen #{queue_number} run took {time.time() - start_time}")
            time.sleep(0.1)


def update_avg_speed_plot_windows(shared_state: TimelapseSharedState):
    # global PIV_OUTPUT
    if shared_state.quick_piv_output_queue.qsize() > 0:
        index, x, y, x_marker, y_marker = shared_state.quick_piv_output_queue.get()
        try:
            PLT_WIDGETS_DICT[index].plot_curve(x, y, x_marker, y_marker)
        except KeyError:
            logging_broadcast(f"window with index {index} doesn't exist")


def update_plotting_windows(shared_state: TimelapseSharedState):
    # global PLOTTING_WINDOW_CREATED
    global PLT_WIDGETS_DICT
    global plotting_windows_timer

    if shared_state.exit_gracefully.is_set():
        if len(PLT_WIDGETS_DICT) > 0:
            for _, window in PLT_WIDGETS_DICT.items():
                window.close()
            plotting_windows_timer.stop()

    if shared_state.plotting_window_created:
        return
    if shared_state.pivjl_script_path:
        if len(PLT_WIDGETS_DICT) > 0:
            # delete window from previous lapse
            PLT_WIDGETS_DICT = {}
        if len(shared_state.specimens_indices_list):
            for i in shared_state.specimens_indices_list:
                PLT_WIDGETS_DICT[i] = PlottingWindow(i)
                PLT_WIDGETS_DICT[i].show()
            shared_state.plotting_window_created = True


def parse_message_from_microscope(message, exit_gracefully: threading.Event):
    time.sleep(0.01)
    try:
        if message.get("type") == "exit":
            logging_broadcast("Recieved terminate command from microscope.")
            exit_gracefully.set()
        elif message.get("type") == "heartbeat":
            # logging_broadcast("Received heartbeat from the microscope.")
            pass
        else:
            logging_broadcast(f"Received message from the microscope: {str(message)}\nBut there is no action for it.")
    except:
        pass


def heartbeat_and_command_handler(port, shared_state: TimelapseSharedState):
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
        while not shared_state.exit_gracefully.is_set():
            next_heartbeat_time = last_heartbeat_sent + HEARTBEAT_INTERVAL_SEC
            time_until_heartbeat = max(0, next_heartbeat_time - time.time())
            # Timeout in milliseconds, max 50ms
            timeout = min(50, time_until_heartbeat * 1000)

            events = dict(poller.poll(timeout))
            if socket in events:
                message = socket.recv_json()
                # logging_broadcast("Received message: " + str(message))
                parse_message_from_microscope(message, shared_state.exit_gracefully)
                if message.get("type") == "heartbeat":
                    last_heartbeat_recieved = time.time()

            if time.time() - last_heartbeat_recieved > HEARTBEAT_FROM_MICROSCOPE_TIMEOUT_sec:
                logging_broadcast(f"Heartbeat from microscope missed, restart")
                shared_state.exit_gracefully.set()

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
        shared_state.exit_gracefully.set()
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


def close_napari_viewer(shared_state: TimelapseSharedState, napari_viewer):
    if shared_state.exit_gracefully.is_set():
        napari_viewer.close()


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
    parser.add_argument('--temp_dir', required=False,
                        default=None,
                        help="Directory path to store input images with incorrect file name")
    parser.add_argument('--pivjl', required=True if '-process_z_projections' in sys.argv else False,
                        default=None,
                        help="Path to auxiliary julia script for average migration speed calculation with quickPIV")
    parser.add_argument('--bot_config', required=False, default=False,
                        help="Path to text file with token and chat_id information (one line, space separator) for"
                             " supplying peaks detection on telegram messenger")
    parser.add_argument('-debug_run', required=False, default=False, action='store_true',
                        help="All content of --output specified folder will be DELETED!")
    parser.add_argument('-process_z_projections',
                        required=True if '--specimen_quantity' in sys.argv else False, action='store_true',
                        help="Process input files as Z-projections for average speed calculations with quickPIV!")
    # Drift correction arguments
    parser.add_argument('-correct_drift',
                        required=True if ('--margin_crop' in sys.argv or "-m" in sys.argv) else False,
                        action='store_true',
                        help="Process input files with dexp drift correction")
    parser.add_argument('-m', '--margin_crop', type=int,
                        required=True if '-correct_drift' in sys.argv else False,
                        help="Value is used for cropping image's margin to get rid of zero valued pixels "
                             "which appear after image shift correction")
    # When the microscope controller restarts the stack_gatherer, this argument defines the configuration file
    # of the unfinished image batch with which the processing will start
    parser.add_argument('-r', '--restart', required=False, default=None,
                        help=" defines the name of the AcquisitionMetadata file for unfinished image batch "
                             "with which the processing will start")

    # event to stop the script gracefully
    # exit_gracefully = threading.Event()
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
    global plotting_windows_timer
    if args.restart:
        global ACQUISITION_METADATA_FILE_TO_PROCESS
        ACQUISITION_METADATA_FILE_TO_PROCESS = args.restart

    # Allow other computers to attach to debugpy at this IP address and port.
    # debugpy.listen(('0.0.0.0', 5680))
    # logging_broadcast("Waiting for debugger attach")
    # debugpy.wait_for_client()  # Blocks execution until client is attached
    heartbeat_thread = None
    manager = Manager()
    shared_state = TimelapseSharedState(manager, args)

    try:
        heartbeat_thread = threading.Thread(
            target=heartbeat_and_command_handler, args=(args.port, shared_state), daemon=True)
        heartbeat_thread.start()
    except Exception as e:
        logging_broadcast(f"{args},\n {e}")
    thread = Thread(target=run_the_loop, args=(vars(args), shared_state, manager), daemon=True)
    thread.start()
    napari_viewer = make_napari_viewer()
    timer1 = QTimer()
    timer1.timeout.connect(lambda: update_avg_speed_plot_windows(shared_state))
    timer1.start(25)
    timer2 = QTimer()
    timer2.timeout.connect(lambda: update_napari_viewer_layer(shared_state, napari_viewer))
    timer2.start(10)
    plotting_windows_timer.timeout.connect(lambda: update_plotting_windows(shared_state))
    plotting_windows_timer.start(1000)
    timer4 = QTimer()
    timer4.timeout.connect(lambda: close_napari_viewer(shared_state, napari_viewer))
    timer4.start(1000)
    napari.run()
    logging_broadcast(f"Napari viewer windows was closed, terminating child processes")
    shared_state.exit_gracefully.set()
    thread.join()
    if heartbeat_thread is not None:
        heartbeat_thread.join()
    logging_broadcast(f"stack_gatherer stopped")


if __name__ == "__main__":
    main()
