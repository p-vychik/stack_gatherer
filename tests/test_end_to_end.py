
from distutils.dir_util import mkpath
import sys
import time
import pytest
import shutil
import threading

 
# setting path
sys.path.append('../stack_gatherer')
from stack_gatherer import *

import numpy
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer
import time
import os
import re
import numpy
from tifffile import imread, imwrite, memmap
from tqdm import tqdm


@pytest.fixture
def setup_image_file_streaming_to_dir():
    images_dir = "test_images"
    mkpath(images_dir)
    acquisition_dir = "test_acquisition"
    mkpath(acquisition_dir)
    num_planes = 30
    shape = (3000, 4000)
    stack_shape = (num_planes,) + shape
    dtype = 'uint8'
    file_list = []
    for plane in range(num_planes):
        for ch in range(4):
            file_path = os.path.join(images_dir, f"TP-0001_SPC-0001_ILL-0_CAM-0_CH-{ch:02}_PL-{plane:04}-outOf-{num_planes:04}.tif")
            imwrite(file_path, numpy.random.randint(0, 2 ** 7, shape, dtype=dtype))
            file_list.append(file_path)
    yield (file_list, acquisition_dir, stack_shape)
    os.remove(images_dir)

    
def image_feeder(file_paths, output_dir):
    for file in file_paths:
        file_name = os.path.basename(file)
        shutil.move(file, os.path.join(output_dir, file_name))
        time.sleep(0.01)

# def test_end_to_end(setup_image_file_streaming_to_dir):
#     file_list_for_feeding, acquisition_dir, stack_shape = setup_image_file_streaming_to_dir

#     saving_thread = threading.Thread(target=run_the_loop(acquisition_dir))
#     saving_thread.start()
#     feeding_thread = threading.Thread(target=image_feeder(file_list_for_feeding, acquisition_dir))
#     feeding_thread.start()
#     feeding_thread.join()
#     stop_file = os.path.join(acquisition_dir, STOP_FILE_NAME)
#     with open(stop_file, 'w') as empty_file:
#         pass  
#     saving_thread.join()
#     os.remove(stop_file)
#     sample_file = ImageFile(file_list_for_feeding[0])
#     sample_file.path_to_image_dir = acquisition_dir
#     for ch in range(stack_shape[0]):
#         sample_file.channel = ch
#         stack_path = sample_file.get_stack_path()
#         assert stack_shape == imread(stack_path).shape
        

def test_spawning_io_bound_function_on_a_separate_thread():
    def io_function():
        # Simulate I/O operations (e.g., reading/writing files or making network requests)
        for _ in range(5):
            print("I/O operation in progress...")
            time.sleep(1)
        print("I/O operation complete")

    # Create a Thread object and pass the io_function as the target
    io_thread = threading.Thread(target=io_function)

    # Start the thread to execute the function concurrently
    io_thread.start()
    for i in range(5):
        print(f"Hell {i}")
    # Optionally, wait for the thread to complete
    io_thread.join()

    # The main thread can continue executing other tasks while the I/O thread runs
    print("Main thread continues...")
    assert True



