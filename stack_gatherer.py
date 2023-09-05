
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
    Example file name: TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001_.tif
    :param file_name: full image file name
    :type file_name: str
    """
    time_point = 0
    specimen = 0
    illumination = 0
    camera = 0
    channel = 0
    plane = 0
    additional_info = ""
    extension = ""

    def __init__(self, file_name):
        split_by = "TP-|SPC-|ILL-|CAM-|CH-|PL-|\."
        name_parts = re.split(split_by, file_name)
        
        if len(name_parts) == 8:
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
                    name_and_info = name_part.strip("-_").split("_", 1)
                    if len(name_and_info) == 1:
                        self.plane = int(plane.strip("-_"))
                    else:
                        plane, info = name_and_info
                        self.plane = int(plane.strip("-_"))
                        self.additional_info = info
                elif i == 7:
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
                f"_CAM-{self.camera}_CH-{self.channel:02}_PL-{self.plane:04}{additional_info}.{self.extension}" 
                )

    def get_name_without_extension(self):
        additional_info = self.additional_info
        if additional_info != "":
            additional_info = "_" + additional_info
        return (f"{self.dataset_name}_TP-{self.time_point:04}"
                f"_SPC-{self.specimen:04}_ILL-{self.illumination}"
                f"_CAM-{self.camera}_CH-{self.channel:02}_PL-{self.plane:04}{additional_info}" 
                )


STOP_FILE_NAME = "STOP_STACK_GATHERING"

def check_if_a_stack_is_ready(filename):
    # Replace this with your logic to check if a stack is ready
    print(f"Checking if stack is ready for file: {filename}")

# Define the event handler class


class MyHandler(FileSystemEventHandler):
    def on_created(self, event):
        if event.is_directory:
            return
        # Call the function when a new file is created
        check_if_a_stack_is_ready(event.src_path)


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
