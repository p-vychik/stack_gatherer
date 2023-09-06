
# start a continuous process from CLI
# Every 100 ms check whether all planes for a full stack have been saved
# We determine whether the stack is done, if new images with a higher timepoint/specimen index,
# or if more than 30 sec has passed since last plane
# 

import argparse
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer
import time
import os
import re
import numpy
from tifffile import imread, imwrite, memmap
from tqdm import tqdm

active_stacks = {}

class ImageFile:
    """
    File naming for light-sheet image files.
    Example file name: TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001-outOf-0150_blaBla.tif
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
        else:
            raise Exception(
                "Image file name is improperly formatted! Check documentation inside the script. Expected 8 parts after splitting by %s" % split_by)

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
        additional_info = self.additional_info
        if additional_info != "":
            additional_info = "_" + additional_info
        return (f"{self.dataset_name}_TP-{self.time_point:04}"
                f"_SPC-{self.specimen:04}_ILL-{self.illumination}"
                f"_CAM-{self.camera}_CH-{self.channel:02}_PL-{self.plane:04}-outOf-{self.total_num_planes:04}{additional_info}" 
                )
    def get_file_path(self):
        return os.path.join(self.path_to_image_dir, self.get_name())
    def get_stack_signature(self):
        return (self.total_num_planes, self.time_point, self.specimen, self.illumination, self.camera)


STOP_FILE_NAME = "STOP_STACK_GATHERING"

def collect_files_to_one_stack(file_list, output_file_name):
    sample_image = imread(file_list[0])
    shape = (len(file_list), sample_image.shape[0], sample_image.shape[1])
    dtype = sample_image.dtype

    # create an empty OME-TIFF file
    imwrite(output_file_name, shape=shape, dtype=dtype, metadata={'axes': 'ZYX'})

    # memory map numpy array to data in OME-TIFF file
    zyx_stack = memmap(output_file_name)
    print(f"Writing stack to {output_file_name}")
    # write data to memory-mapped array
    with tqdm(total=len(file_list), desc="Saving plane") as pbar:
        for t in range(shape[0]):
            if t == 0:
                zyx_stack[t] = sample_image
                zyx_stack.flush()
                pbar.update(1)
                continue
            zyx_stack[t] = imread(file_list[t])
            zyx_stack.flush()
            pbar.update(1)



def add_file_to_active_stacks(file_path):
    file = ImageFile(file_path)
    stack_signature = file.get_stack_signature()
    if stack_signature not in active_stacks:
        active_stacks[stack_signature] = {}
        print(f"Adding stack {stack_signature} to active queue.")
    if file.plane not in active_stacks[stack_signature]:
        active_stacks[stack_signature][file.plane] = file 
    return stack_signature

def check_if_a_stack_is_ready(stack_signature):
    if len(active_stacks[stack_signature]) < stack_signature[0]:
        return
    file_list = []
    for i, _ in enumerate(active_stacks[stack_signature]):
        # We have to access by index since we can't gurantee that files were added to dict in order of planes
        file_list.append(active_stacks[stack_signature][i])
    

# Define the event handler class


class MyHandler(FileSystemEventHandler):
    def on_created(self, event):
        if event.is_directory:
            return
        file_path =  event.src_path
        # Call the function when a new file is created
        stack_signature = add_file_to_active_stacks(file_path)
        check_if_a_stack_is_ready(stack_signature)

def run_the_loop(input_dir):
    # Create an observer and attach the event handler
    observer = Observer()
    observer.schedule(MyHandler(), path=input_dir, recursive=False)

    # Start the observer
    observer.start()

    try:
        stop_file = os.path.join(input_dir, STOP_FILE_NAME)
        while (not os.path.exists(stop_file)):
            time.sleep(1)  # Sleep to keep the script running
    except KeyboardInterrupt:
        # Gracefully stop the observer if the script is interrupted
        observer.stop()

    # Wait for the observer to complete
    observer.join()

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Watch a directory for new file additions.")
    parser.add_argument("directory", help="The directory to watch for new files.")

    # Parse the command line arguments
    args = parser.parse_args()
    run_the_loop(args.directory)




if __name__ == "__main__":
    main()
