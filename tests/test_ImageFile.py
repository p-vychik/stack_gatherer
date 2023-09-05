
import sys
 
# setting path
sys.path.append('../stack_gatherer')
from stack_gatherer import ImageFile

def test_same_filename_is_returned():
    input = "dataseTPtName_TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001_ffdfs.tif"
    assert input == ImageFile(input).get_name()

def test_same_filename_withot_extension_is_returned():
    input = "dataseTPtName_TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001_ffdfs.tif"
    input_no_extension = "dataseTPtName_TP-0001_SPC-0001_ILL-0_CAM-0_CH-01_PL-0001_ffdfs"
    assert input_no_extension == ImageFile(input).get_name_without_extension()

def test_params_values():
    input = "dataseTPtName_TP-0001_SPC-0101_ILL-3_CAM-2_CH-11_PL-0777_ffdfs.tif"
    obj = ImageFile(input)
    assert obj.dataset_name == "dataseTPtName"
    assert obj.time_point == 1
    assert obj.specimen == 101
    assert obj.illumination == 3
    assert obj.camera == 2
    assert obj.channel == 11
    assert obj.plane == 777
    assert obj.additional_info == "ffdfs"
    assert obj.extension == "tif"