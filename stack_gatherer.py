
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
                        self.total_num_planes = int(name_and_info.strip("-_"))
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


STOP_FILE_NAME = "STOP_STACK_GATHERING"

def add_file_to_active_stacks(file_path):
    file_obj = ImageFile(os.path.basename(file_path))
    stack_signature = ()

def check_if_a_stack_is_ready(file_path):
    # Replace this with your logic to check if a stack is ready
    print(f"Checking if stack is ready for file: {file_path}")

# Define the event handler class


class MyHandler(FileSystemEventHandler):
    def on_created(self, event):
        if event.is_directory:
            return
        file_path =  event.src_path
        # Call the function when a new file is created
        add_file_to_active_stacks(file_path)
        check_if_a_stack_is_ready(file_path)


def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Watch a directory for new file additions.")
    parser.add_argument("directory", help="The directory to watch for new files.")

    # Parse the command line arguments
    args = parser.parse_args()

    # Create an observer and attach the event handler
    observer = Observer()
    observer.schedule(MyHandler(), path=args.directory, recursive=False)

    # Start the observer
    observer.start()

    try:
        stop_file = os.path.join(args.directory, STOP_FILE_NAME)
        while (not os.path.exists(stop_file)):
            time.sleep(1)  # Sleep to keep the script running
    except KeyboardInterrupt:
        # Gracefully stop the observer if the script is interrupted
        observer.stop()

    # Wait for the observer to complete
    observer.join()


if __name__ == "__main__":
    main()
