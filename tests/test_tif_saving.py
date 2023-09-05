
import sys
import time

 
# setting path
sys.path.append('../stack_gatherer')
from stack_gatherer import collect_files_to_one_stack

import numpy
from tifffile import imread, imwrite
import os

def test_collecting_files_to_stack():
    num_planes = 50
    shape = (3000, 4000)
    dtype = 'uint8'
    output_file = "output_stack.ome.tif"
    file_list = []
    for i in range(num_planes):
        file_name = f"plane_{i}.tif"
        imwrite(file_name, numpy.random.randint(0, 2 ** 7, shape, dtype=dtype))
        file_list.append(file_name)
    start_time = time.time()
    collect_files_to_one_stack(file_list, output_file)

    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Execution time: {execution_time} seconds")

    assert imread(output_file).shape == (num_planes,) + shape
    for file_path in file_list:
        try:
            os.remove(file_path)
            print(f"Deleted: {file_path}")
        except Exception as e:
            print(f"Error deleting {file_path}: {e}")