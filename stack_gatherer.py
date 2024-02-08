# start a continuous process from CLI
# Every 100 ms check whether all planes for a full stack have been saved
# We determine whether the stack is done, if new images with a higher timepoint/specimen index,
# or if more than 30 sec has passed since last plane
#
import math
import argparse
import sys
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer
import time
import os
import re
import numpy as np
from tifffile import imread, imwrite, memmap
from tqdm import tqdm
import threading
from PIL import Image
import napari
from napari.qt.threading import thread_worker
from threading import Thread
from queue import Queue
from collections import OrderedDict
from qtpy.QtWidgets import QCheckBox
from dataclasses import dataclass
import shutil
import json
import csv
# for PIV results visualisation with matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
from datetime import datetime
np.set_printoptions(threshold=np.inf)

STOP_FILE_NAME = "STOP_STACK_GATHERING"
ACQUISITION_META_FILE_PATTERN = "AcquisitionMetadata_"
# BORDER_WIDTH defines the size of the border in pixels between
# projections merged in one image for napari visualization
BORDER_WIDTH = 20
# queue for storing dictionaries with projections generated in collect_files_to_one_stack_get_axial_projections()
MAX_SPEED = 1
PROCESSED_Z_PROJECTIONS = Queue()
PROJECTIONS_QUEUE = OrderedDict()
DRAWN_PROJECTIONS_QUEUE = OrderedDict()
JSON_CONFIGS = {}
MERGE_LIGHT_MODES = False
VIEWER = None
TEMP_DIR = None
PIVJLPATH = ""
plt.rcParams.update({'figure.autolayout': True})
style.use('fivethirtyeight')
FIG = plt.figure()
DYNAMIC_CANVAS = FigureCanvas(Figure(figsize=(5, 3)))
# added due to UserWarning: Tight layout not applied

PIV_PLOT_CANVAS = DYNAMIC_CANVAS.figure.subplots()
PIL_DATA_STORAGE = OrderedDict()
active_stacks = {}
currenty_saving_stacks_locks = {}
currenty_saving_stacks_dict_lock = threading.Lock()


class NotImagePlaneFile(Exception):
    pass


@dataclass
class StackSignature:
    total_num_planes: int
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
        return self.total_num_planes, self.time_point, self.specimen, self.illumination, self.camera, self.channel

    def __hash__(self):
        return hash(self.get_attr_excluding())

    @property
    def signature(self):
        return self.get_attr_excluding()


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


class ImageFile:
    """
    File naming for light-sheet image files.
    Example file name: TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001-outOf-0150_blaBla.tif or .bmp
    :param file_path: full path to image
    :type file_path: str
    """
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
        split_by = r"TP-|SPC-|ILL-|CAM-|CH-|PL-|outOf-|\."
        name_parts = re.split(split_by, file_name)
        
        if len(name_parts) == 9:
            try:
                for i, name_part in enumerate(name_parts):
                    
                    if i == 0:
                        self.dataset_name = name_part.strip("-_")
                    elif i == 1:
                        self.time_point = int(name_part.strip("-_"))
                    elif i == 2:
                        self.specimen = int(name_part.strip("-_"))
                    elif i == 3:
                        self.illumination = int(name_part.strip("-_"))
                    elif i == 4:
                        self.camera = int(name_part.strip("-_"))
                    elif i == 5:
                        self.channel = int(name_part.strip("-_"))
                    elif i == 6:
                        self.plane = int(name_part.strip("-_"))
                    elif i == 7:
                        name_and_info = name_part.strip("-_").split("_", 1)
                        if len(name_and_info) == 1:
                            self.total_num_planes = int(name_and_info[0].strip("-_"))
                        else:
                            num_planes, info = name_and_info
                            self.total_num_planes = int(num_planes.strip("-_"))
                            self.additional_info = info
                    elif i == 8:
                        self.extension = name_part
            except ValueError:
                raise NotImagePlaneFile(
                    "This is not a valid filename for a single plane image file!"
                )
        else:
            raise NotImagePlaneFile(
                "Image file name is improperly formatted! Check documentation inside the script. "
                "Expected 8 parts after splitting by %s" % split_by)

        self.extension.lower()

    def get_name(self):
        additional_info = self.additional_info
        dataset_name = self.dataset_name
        if additional_info != "":
            additional_info = "_" + additional_info
        if dataset_name != "":
            dataset_name = dataset_name + "_"
        return (f"{dataset_name}TP-{self.time_point:04}"
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
        return (f"{dataset_name}TP-{self.time_point:04}"
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
                              self.time_point,
                              self.specimen,
                              self.illumination,
                              self.camera,
                              self.channel)


def read_image(image_path):
    if image_path.endswith((".tif", ".tiff")):
        return imread(image_path)
    if image_path.endswith(".bmp"):
        image = Image.open(image_path)
        return np.array(image)
    return False


def plane_to_projection(plane, output_dictionary):

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


def collect_files_to_one_stack_get_axial_projections(stack_signature,
                                                     file_list,
                                                     output_file_path,
                                                     axes=None,
                                                     anisotropy_factor=1,
                                                     output_dir=None):
    if os.path.exists(output_file_path):
        return
    sample_image = read_image(file_list[0])
    shape = (len(file_list), sample_image.shape[0], sample_image.shape[1])
    dtype = sample_image.dtype
    # create an empty OME-TIFF file
    imwrite(output_file_path, shape=shape, dtype=dtype, metadata={'axes': 'ZYX'})
    # memory map numpy array to data in OME-TIFF file
    zyx_stack = memmap(output_file_path)
    print(f"Writing stack to {output_file_path}")
    # prepare input about required projections in the dictionary,
    # the values of the appropriate keys would be the projection matrices
    projections = {}
    projections_files_path = {}
    if axes:
        if ',' in axes:
            projections = {axis.upper(): None for axis in set(axes) if axis in "zZxXyY"}
        elif axes in "zZxXyY":
            projections = {axes.upper().strip(): None}
        if projections:
            try:
                file_name = os.path.basename(file_list[0]).replace('.bmp', '.tif')
                for axis in projections.keys():
                    projection_folder_path = os.path.join(output_dir, f"{axis}_projections")
                    if not os.path.exists(projection_folder_path):
                        os.makedirs(projection_folder_path)
                    projections_files_path[axis] = os.path.join(projection_folder_path, file_name)
            except OSError as err:
                print(err)
    # write data to memory-mapped array
    with tqdm(total=len(file_list), desc="Saving plane") as pbar:
        for z in range(shape[0]):
            if z == 0:
                zyx_stack[z] = sample_image
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
            PROCESSED_Z_PROJECTIONS.put(wrapped_z_projection)
            # PROCESSED_Z_PROJECTIONS.put(projections["Z"])

        for axis in projections.keys():
            try:
                imwrite(projections_files_path[axis], projections[axis])
            except IOError as err:
                print(err)

        # we want to store projections derived from planes with different illumination modes to
        # have an option to merge them if the user checked the QCheckBox() in napari viewer GUI interface

        wrapped_projections = ProjectionsDictWrapper(projections, stack_signature)
        if wrapped_projections.signature not in PROJECTIONS_QUEUE:
            PROJECTIONS_QUEUE[wrapped_projections.identifier] = [wrapped_projections]
        else:
            PROJECTIONS_QUEUE[wrapped_projections.identifier].append(wrapped_projections)

    for z in range(shape[0]):
        os.remove(file_list[z])


def add_file_to_active_stacks(image_file: ImageFile):

    stack_signature = image_file.get_stack_signature()
    if stack_signature.signature not in active_stacks:
        active_stacks[stack_signature.signature] = {}
        print(f"Adding stack {stack_signature.signature} to active queue.")
    if image_file.plane not in active_stacks[stack_signature.signature]:
        active_stacks[stack_signature.signature][image_file.plane] = image_file
    return stack_signature


def check_stack_and_collect_if_ready(stack_signature, output_dir, axes=None, factor=None):
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
                                                     axes=axes,
                                                     anisotropy_factor=factor,
                                                     output_dir=output_dir)
    del active_stacks[stack_signature.signature]
    del currenty_saving_stacks_locks[stack_signature.signature]


# Define the event handler class
class MyHandler(FileSystemEventHandler):
    output_dir = ""
    factor = 1
    axes = ""

    def __init__(self, output_dir, factor, axes):
        self.output_dir = output_dir
        self.factor = int(factor)
        self.axes = axes

    def on_created(self, event):
        if event.is_directory:
            return
        file_path = event.src_path
        # load json parameters
        if os.path.basename(file_path).startswith(ACQUISITION_META_FILE_PATTERN) and file_path.endswith(".json"):
            destination_folder = os.path.join(self.output_dir, os.path.basename(file_path).strip(".json"))
            lapse_parameters = None
            try:
                j_file = open(file_path)
                lapse_parameters = json.load(j_file)
                j_file.close()
            except PermissionError:
                print(f"{file_path} permission error, check if file is already opened")
            if lapse_parameters:
                lapse_parameters["output_folder"] = destination_folder
                setup_signature = []
                # in the filename specimen's index usually starts from 1
                for specimen_number, spec_entry in enumerate(lapse_parameters["specimens"], start=1):
                    timepoint = lapse_parameters["timelapse"]["timepoints"]
                    total_num_planes = spec_entry["number_of_planes"]
                    for channel_num, channel in enumerate(spec_entry["channels"]):
                        for camera_num, enabled in enumerate(channel["camerasEnabled"]):
                            if enabled:
                                setup_signature.append((total_num_planes,
                                                        timepoint,
                                                        specimen_number,
                                                        channel["illuminationLine"],
                                                        camera_num,
                                                        channel_num))
                JSON_CONFIGS.update(dict.fromkeys(setup_signature, lapse_parameters))
            try:
                if not os.path.exists(destination_folder):
                    os.mkdir(destination_folder)
                    shutil.move(file_path, destination_folder)
            except Exception as error:
                print(error)
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
            # while file.get_stack_signature() not in JSON_CONFIGS:
            #     time.sleep(0.1)
            stack_signature = add_file_to_active_stacks(file)
            # create separate folder for the output based on metadata filename
            try:
                output_dir = JSON_CONFIGS[stack_signature.signature]["output_folder"]
                if os.path.exists(output_dir):
                    check_stack_and_collect_if_ready(stack_signature, output_dir, self.axes, self.factor)
            except KeyError:
                print("No lapse configuration for the plane, check if appropriate AcquisitionMeta file is loaded")


def run_the_loop(kwargs):
    global TEMP_DIR
    input_dir = kwargs.get("input")
    output_dir = kwargs.get("output")
    axes = kwargs.get("axes", "Z")
    factor = kwargs.get("factor_anisotropy", None)
    TEMP_DIR = kwargs.get("temp_dir", None)

    # Create an observer and attach the event handler
    observer = Observer()
    observer.schedule(MyHandler(output_dir=output_dir, factor=factor, axes=axes), path=input_dir, recursive=False)
    # Start the observer
    observer.start()
    print(f"Watching {input_dir} for images, and saving stacks to {output_dir}")

    try:
        stop_file = os.path.join(input_dir, STOP_FILE_NAME)
        while (not os.path.exists(stop_file)):
            time.sleep(1)  # Sleep to keep the script running
    except KeyboardInterrupt:
        # Gracefully stop the observer if the script is interrupted
        observer.stop()
    # save PIV data to csv
    if PIL_DATA_STORAGE:
        with open(os.path.join(output_dir, "quckPIV_data.csv"), 'w') as f_out:
            w = csv.writer(f_out)
            try:
                w.writerows(PIL_DATA_STORAGE.items())
            except IOError:
                print(f"Attempt to save csv data to {output_dir} failed")
    # Wait for the observer to complete
    observer.join()


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


def update_layer(layers_dict):
    for axes_names, image in layers_dict.items():
        if axes_names not in VIEWER.layers:
            VIEWER.add_image(image, name=axes_names)
        else:
            VIEWER.layers[axes_names].data = image


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


@thread_worker(connect={'yielded': update_layer})
def get_projections_dict_from_queue():
    global DRAWN_PROJECTIONS_QUEUE
    while True:
        time.sleep(0.2)
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
                    print("Empty projections dictionary")
            if image_layer_dict:
                yield image_layer_dict


def bx_trigger():
    global MERGE_LIGHT_MODES
    MERGE_LIGHT_MODES = not MERGE_LIGHT_MODES


def make_napari_viewer():
    global VIEWER
    VIEWER = napari.Viewer()
    bx = QCheckBox('Merge illumination channels')
    bx.setChecked(MERGE_LIGHT_MODES)
    bx.stateChanged.connect(bx_trigger)
    VIEWER.window.add_dock_widget(bx)
    if PIVJLPATH:
        PIV_PLOT_CANVAS.plot(0,0)
        VIEWER.window.add_dock_widget(DYNAMIC_CANVAS, area='bottom', name='PIL data')


def plot_pil_data(piv_data):
    if piv_data:
        y_range = math.ceil(MAX_SPEED * 1.5)
        PIV_PLOT_CANVAS.clear()
        PIV_PLOT_CANVAS.set_ylim([0, y_range])
        x, y = zip(*piv_data.items())
        PIV_PLOT_CANVAS.plot(x, y)
        PIV_PLOT_CANVAS.set_xticklabels(x, fontsize=10, rotation=30)
        # for tick in PIV_PLOT_CANVAS.get_xticklabels():
        #     tick.set_rotation(30)
        PIV_PLOT_CANVAS.figure.canvas.draw()
        FIG.tight_layout()


def PIV_worker():
    if PIVJLPATH:
        from juliacall import Main as jl
        jl.include(PIVJLPATH)
        global PIL_DATA_STORAGE
        global MAX_SPEED
        piv_projection_queue = Queue()
        while True:
            time.sleep(0.1)
            if PROCESSED_Z_PROJECTIONS.qsize() >= 2:
                wrapped_z_p_1 = PROCESSED_Z_PROJECTIONS.get()
                # to calculate avg speed in consecutive way we store
                # the last projection from every pairwise comparison
                wrapped_z_p_2 = PROCESSED_Z_PROJECTIONS.queue[0]
                
                if wrapped_z_p_1.identifier == wrapped_z_p_2.identifier:
                    merged_projections = merge_multiple_projections((wrapped_z_p_1.identifier,
                                                                     (wrapped_z_p_1, wrapped_z_p_2)))
                    merged_wrapped_proj = ProjectionsDictWrapper(merged_projections, wrapped_z_p_2.signature)
                    piv_projection_queue.put(merged_wrapped_proj)
                    PROCESSED_Z_PROJECTIONS.get()
                else:
                    piv_projection_queue.put(wrapped_z_p_1)

                if piv_projection_queue.qsize() < 2: 
                    continue

                m_1 = piv_projection_queue.get().projections["Z"]
                m_2 = piv_projection_queue.get().projections["Z"]
                if m_1.shape == m_2.shape:
                    avg_speed = jl.fn(m_1, m_2)
                    avg_speed = round(avg_speed[-1], 3)
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    if avg_speed:
                        MAX_SPEED = max(avg_speed, MAX_SPEED)
                    PIL_DATA_STORAGE[current_time] = avg_speed
                    if len(PIL_DATA_STORAGE) > 15:
                        # PIL_DATA_STORAGE.popitem(last=False)
                        plot_pil_data(PIL_DATA_STORAGE[-15::])
                    else:
                        plot_pil_data(PIL_DATA_STORAGE)
                else:
                    print("Projections should have the same size for QuickPIV input")


def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Watch a directory for new file additions and collect .tif or .bmp "
                                                 "files as stacks to output directory.")
    # Define the input_folder argument
    parser.add_argument('-i', '--input', required=True, help="Input folder path to watch.")
    
    # Define the output_folder argument
    parser.add_argument('-o', '--output', required=True, help="Output folder path to save stacks.")

    # Define axes to make projections
    parser.add_argument('-a', '--axes', required=False, help="Comma separated axes to project "
                                                             "on the plane, i.e. X,Y,Z or X, or X")

    # Anisotropy factor correction
    parser.add_argument('-f', '--factor_anisotropy', required=True if '-a' in sys.argv else False,
                        help="Value is used for correcting projection's anisotropic distortions")
    # Move image files with incorrect name to the user provided directory
    parser.add_argument('-d', '--temp_dir', required=True,
                        help="Directory path to store input images with incorrect file name")
    parser.add_argument('--pivjl', required=False, default=False,
                        help="Path to auxiliary julia script for average migration speed calculation with quickPIV")
    parser.add_argument('-debug_run', required=False, default=False, action='store_true',
                        help="All content of --output specified folder will be DELETED!")
    global PIVJLPATH
    # Parse the command-line arguments
    args = parser.parse_args()
    if args.pivjl:
        if os.path.isfile(args.pivjl):
            PIVJLPATH = args.pivjl
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
                print(e)

    make_napari_viewer()
    if PIVJLPATH:
        piv_worker_thread = Thread(target=PIV_worker)
        piv_worker_thread.start()
    get_projections_dict_from_queue()
    thread = Thread(target=run_the_loop, args=(vars(args), ))
    thread.start()
    napari.run()
    thread.join()


if __name__ == "__main__":
    main()
