
import sys
import os

 
# setting path
sys.path.append('../stack_gatherer')
from stack_gatherer import ImageFile

def test_same_filename_is_returned():
    input = "dataseTPtName_timelapseID-20240214-192141_TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001-outOf-0150_ffdfs.tif"
    assert input == ImageFile(input).get_name()
    input = "dataseTPtName_timelapseID-20240214-192141_TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001-outOf-0150.tif"
    assert input == ImageFile(input).get_name()
    input = "timelapseID-20240214-192141_TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001-outOf-0150.tif"
    assert input == ImageFile(input).get_name()

def test_same_filename_withot_extension_is_returned():
    input = "dataseTPtName_timelapseID-20240214-192141_TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001-outOf-0150_ffdfs.tif"
    input_no_extension = "dataseTPtName_timelapseID-20240214-192141_TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001-outOf-0150_ffdfs"
    assert input_no_extension == ImageFile(input).get_name_without_extension()

def test_stack_name_is_returned_correctly():
    input = "dataseTPtName_timelapseID-20240214-192141_TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001-outOf-0150_ffdfs.tif"
    stack_name = "dataseTPtName_timelapseID-20240214-192141_TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-(ZS)-outOf-0150_ffdfs.tif"
    assert stack_name == ImageFile(input).get_stack_name()
    
def test_params_values():
    input = r"foo/bar/dataseTPtName_timelapseID-20240214-192141_TP-0001_SPC-0101_ILL-3_CAM-2_CH-11_PL-0777-outOf-1000_ffdfs.tif"
    obj = ImageFile(input)
    assert obj.dataset_name == "dataseTPtName"
    assert obj.timelapse_id == "20240214-192141"
    assert obj.time_point == 1
    assert obj.specimen == 101
    assert obj.illumination == 3
    assert obj.camera == 2
    assert obj.channel == 11
    assert obj.plane == 777
    assert obj.total_num_planes == 1000
    assert obj.additional_info == "ffdfs"
    assert obj.extension == "tif"
    assert obj.path_to_image_dir == r"foo/bar"

def test_same_full_path_is_returned():
    file_name = "dataseTPtName_timelapseID-20240214-192141_TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001-outOf-0150_ffdfs.tif"
    input = os.path.join("foo", "bar", file_name)
    assert input == ImageFile(input).get_file_path()

