# start a continuous process from CLI
# Every 100 ms check whether all planes for a full stack have been saved
# We determine whether the stack is done, if new images with a higher timepoint/specimen index,
# or if more than 30 sec has passed since last plane
# 

import argparse
import sys
from copy import deepcopy
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
# from queue import Queue
from collections import OrderedDict
from qtpy.QtWidgets import QCheckBox
import shutil
# for PIV results visualisation with matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
import random
from datetime import datetime

STOP_FILE_NAME = "STOP_STACK_GATHERING"
# BORDER_WIDTH defines the size of the border in pixels between
# projections merged in one image for napari visualization
BORDER_WIDTH = 20
# queue for storing dictionaries with projections generated in collect_files_to_one_stack_get_axial_projections()
# PROJECTIONS_QUEUE = Queue()
PROJECTIONS_QUEUE = OrderedDict()
MERGE_LIGHT_MODES = False
VIEWER = None
TEMP_DIR = None

plt.rcParams.update({'figure.autolayout': True})
style.use('fivethirtyeight')
FIG = plt.figure()
DYNAMIC_CANVAS = FigureCanvas(Figure(figsize=(5, 3)))
AX1 = DYNAMIC_CANVAS.figure.subplots()
PIL_DATA_STORAGE = OrderedDict()

active_stacks = {}
currenty_saving_stacks_locks = {}
currenty_saving_stacks_dict_lock = threading.Lock()


class NotImagePlaneFile(Exception):
    pass

class ProjectionsDictWrapper:
    projections = {}
    signature = None
    identifier = None
    illumination = None

    def __init__(self, dictionary, stack_signature):
        self.projections = dictionary.copy()
        self.signature = tuple([val for val in stack_signature])
        self.identifier = (stack_signature[0], stack_signature[1], stack_signature[2], stack_signature[4])
        self.illumination = stack_signature[3]


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
        split_by = "TP-|SPC-|ILL-|CAM-|CH-|PL-|outOf-|\."
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
        return (self.total_num_planes, self.time_point, self.specimen, self.illumination, self.camera, self.channel)


def read_image(image_path):
    if image_path.endswith((".tif", ".tiff")):
        return imread(image_path)
    if image_path.endswith(".bmp"):
        image = Image.open(image_path)
        return np.array(image)
    return False


def plane_to_projection(plane, output_dictionary):

    '''
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
    '''

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
        # .fromarray() mode parameter defines the type and depth of the pixel, here 8-bit grayscale is set with mode="L"
        if "Y" in projections:
            projection_y = projections["Y"]
            projection_y = projection_y.transpose()

            img = Image.fromarray(projection_y, mode="L")
            projections["Y"] = np.array(img.resize(size=(projections["Y"].shape[0] * anisotropy_factor,
                                                         projections["Y"].shape[1])))
        if "X" in projections:
            projection_x = projections["X"]
            img = Image.fromarray(projection_x, mode="L")
            projections["X"] = np.array(img.resize(size=(projections["X"].shape[1],
                                                         projections["X"].shape[0] * anisotropy_factor)))
        for axis in projections.keys():
            try:
                imwrite(projections_files_path[axis], projections[axis])
            except IOError as err:
                print(err)

        # we want to store projections derived from planes with different illumination modes to
        # have an option to merge them if the user checked the QCheckBox() in napari viewer interface

        wrapped_projections = ProjectionsDictWrapper(projections, stack_signature)
        if wrapped_projections.identifier not in PROJECTIONS_QUEUE:
            PROJECTIONS_QUEUE[wrapped_projections.identifier] = [wrapped_projections]
        else:
            PROJECTIONS_QUEUE[wrapped_projections.identifier].append(wrapped_projections)


    for z in range(shape[0]):
        os.remove(file_list[z])


def add_file_to_active_stacks(image_file: ImageFile):

    stack_signature = image_file.get_stack_signature()
    if stack_signature not in active_stacks:
        active_stacks[stack_signature] = {}
        print(f"Adding stack {stack_signature} to active queue.")
    if image_file.plane not in active_stacks[stack_signature]:
        active_stacks[stack_signature][image_file.plane] = image_file 
    return stack_signature


def check_stack_and_collect_if_ready(stack_signature, output_dir, axes=None, factor=None):
    if len(active_stacks[stack_signature]) < stack_signature[0]:
        return
    # We have to ensure that two events firing at the same time don't start saving the same stack twice
    with currenty_saving_stacks_dict_lock:
        if stack_signature not in currenty_saving_stacks_locks:
            currenty_saving_stacks_locks[stack_signature] = True
        else:
            return
    file_list = []
    for i, _ in enumerate(active_stacks[stack_signature]):
        # We have to access by index since we can't guarantee that files were added to dict in order of planes
        file_list.append(active_stacks[stack_signature][i].get_file_path())
    sample_file_obj = ImageFile(file_list[0])
    sample_file_obj.extension = "tif"
    stack_path = os.path.join(output_dir, sample_file_obj.get_stack_name())
    collect_files_to_one_stack_get_axial_projections(stack_signature,
                                                    file_list,
                                                    stack_path,
                                                    axes=axes,
                                                    anisotropy_factor=factor,
                                                    output_dir=output_dir
                                                    )
    del active_stacks[stack_signature]
    del currenty_saving_stacks_locks[stack_signature]


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
        check_stack_and_collect_if_ready(stack_signature, self.output_dir, self.axes, self.factor)


def run_the_loop(kwargs):
    input_dir = kwargs.get("input")
    output_dir = kwargs.get("output")
    axes = kwargs.get("axes", None)
    factor = kwargs.get("factor_anisotropy", None)
    global TEMP_DIR
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

    # Wait for the observer to complete
    observer.join()


def draw_napari_layer(projections_dict):
    # define the padding between pictures
    # zero arrays as vertical and horizontal borders
    # projections_dict = wrapped_projections_dict.projections
    image = None
    image_dict = {}
    axes_names = ''.join([axis_key for axis_key, _ in projections_dict.items()])
    if len(projections_dict) == 3:
        v_border_array = np.zeros((projections_dict["Z"].shape[0], BORDER_WIDTH))
        # zero array for horizontal border
        # rows number equals to the BORDER_WIDTH value, columns to width of concatenation of Z, Y and border array)
        h_border_array = np.zeros((BORDER_WIDTH,
                                   projections_dict["Z"].shape[1] + BORDER_WIDTH + projections_dict["Y"].shape[1]))
        # extend Z projection with border and Y projection arrays
        z_y = np.hstack((projections_dict["Z"], v_border_array, projections_dict["Y"]))
        # merge Z_Y with horizontal border array
        z_y = np.vstack((z_y, h_border_array))
        x = np.hstack((projections_dict["X"],
                       np.zeros((projections_dict["X"].shape[0], projections_dict["Y"].shape[1] + BORDER_WIDTH))))
        image = np.vstack((z_y, x))
    elif len(projections_dict) == 2:
        # place largest projection in center and arrange the second at the appropriate side
        # if only X and Y - then arrange them in a perpendicular way
        if "Z" in new_image and "X" in new_image:
            h_border_array = np.zeros((BORDER_WIDTH, projections_dict["Z"].shape[1]))
            image = np.vstack((projections_dict["Z"], h_border_array, projections_dict["X"]))
        elif "Z" in projections_dict and "Y" in projections_dict:
            v_border_array = np.zeros((projections_dict["Z"].shape[0], BORDER_WIDTH))
            image = np.hstack((projections_dict["Z"], v_border_array, projections_dict["Y"]))
        else:
            # only X and Y projections, arrange them perpendicular
            dummy_array = np.zeros((projections_dict["Y"].shape[0], projections_dict["X"].shape[1]))
            dummy_y = np.hstack((dummy_array, projections_dict["Y"]))
            x_extended = np.hstack((projections_dict["X"],
                                    np.zeros((projections_dict["X"].shape[0],
                                              dummy_y.shape[1] - projections_dict["X"].shape[1]))))
            image = np.vstack((dummy_y, x_extended))
    elif len(projections_dict) == 1:
        image = projections_dict[axes_names]
    image_dict[axes_names] = image
    return image_dict


def update_layer(layers_dict):
    for axes_names, image in layers_dict.items():
    #if isinstance(image, np.ndarray):
        if axes_names not in VIEWER.layers:
            VIEWER.add_image(image, name=axes_names)
        else:
            VIEWER.layers[axes_names].data = image

    #else:
        #print("Concatenating projections failed!")


def merge_multiple_projections(wrapped_dict_list: list):
    # as an input a list of wrapped dictionaries is provided - check class ProjectionsDictWrapper
    # we have to merge appropriate projections according to the maximum intensity
    merged_projections = {}
    for wrapped_dict in wrapped_dict_list[-1]:
        for axis, plane in wrapped_dict.projections.items():
            if axis not in merged_projections:
                merged_projections[axis] = plane
            else:
                array = np.stack((merged_projections[axis], plane))
                merged_projections[axis] = array.max(axis=0)
    return merged_projections


@thread_worker(connect={'yielded': update_layer})
def get_projections_dict_from_queue():
    while True:
        call_pill()
        plot_pil_data()
        time.sleep(0.5)
        #if not PROJECTIONS_QUEUE.empty():
        #    yield PROJECTIONS_QUEUE.get()
        if len(PROJECTIONS_QUEUE):
            if MERGE_LIGHT_MODES:
                _, first_entry = next(iter(PROJECTIONS_QUEUE.items()))
                # should check the number of light channels
                if len(first_entry) == 2:
                    projections_dict = merge_multiple_projections(PROJECTIONS_QUEUE.popitem(last=False))
                    layer_image_dict = draw_napari_layer(projections_dict)
                    yield layer_image_dict
            else:
                _, first_entry = next(iter(PROJECTIONS_QUEUE.items()))
                if len(first_entry) > 1:
                    # iterate through the list, assign illumination index to the axis, save all to the new ,
                    # data remain unmerged - apply merge and push to napari
                    _, first_entry_list = PROJECTIONS_QUEUE.popitem(last=False)
                    image_layer_dict = {}
                    for wrapped_dict in first_entry_list:
                        channel_layer = draw_napari_layer(wrapped_dict.projections)
                        for key, val in channel_layer.items():
                            image_layer_dict[f"{key}_ILL_{wrapped_dict.illumination}"] = val
                    yield image_layer_dict
                    # projections_dict = merge_multiple_projections(wrapped_dict)
                    # yield projections_dict
                else:
                    pass
                    # layer_image_dict = draw_napari_layer(first_entry[0].projections)
                    # yield layer_image_dict
                    # yield projections_list[0].projections


def bx_trigger():
    global MERGE_LIGHT_MODES
    MERGE_LIGHT_MODES = not MERGE_LIGHT_MODES


def make_napari_viewer():
    global VIEWER
    VIEWER = napari.Viewer()
    bx = QCheckBox('Merge illumination')
    bx.setChecked(MERGE_LIGHT_MODES)
    bx.stateChanged.connect(bx_trigger)
    VIEWER.window.add_dock_widget(bx)
    AX1.plot()
    VIEWER.window.add_dock_widget(DYNAMIC_CANVAS, area='bottom', name='PIL data')


def plot_pil_data():

    if PIL_DATA_STORAGE:
        AX1.clear()
        x,y = zip(*PIL_DATA_STORAGE.items())
        AX1.plot(x, y)
        for tick in AX1.get_xticklabels():
            tick.set_rotation(45)
        AX1.figure.canvas.draw()


def call_pill():
    global PIL_DATA_STORAGE
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    PIL_DATA_STORAGE[current_time] = random.randint(1,20)
    if len(PIL_DATA_STORAGE) > 40:
        PIL_DATA_STORAGE.popitem(last=False)


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
    parser.add_argument('-d', '--temp_dir', required=False,
                        help="Directory path to store input images with incorrect file name")
    # Parse the command-line arguments
    args = parser.parse_args()
    make_napari_viewer()
    get_projections_dict_from_queue()
    thread = Thread(target=run_the_loop, args=(vars(args), ))
    thread.start()
    napari.run()
    thread.join()


if __name__ == "__main__":
    main()
