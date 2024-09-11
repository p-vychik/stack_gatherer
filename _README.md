INSTALL

1. install anaconda for local user
https://www.anaconda.com/

2. create separate environment (tested with python 3.10 and julia 1.10)
conda create -n pivenv python=3.10

3. activate the created environment
conda activate pivenv

4. install julia manager
conda install -c conda-forge juliaup

5. choose julia version
juliaup add 1.10~x64

6. install python packages
pip install watchdog matplotlib juliacall napari PyQt5

7. run julia in the created environment, in REPL mode set python path:

ENV["PYTHON"] = "C:\Users\..\anaconda3\envs\pivenv\python.exe"
ENV["JULIA_CONDAPKG_BACKEND"] = "Current"

8. run package manager by pressing the key "]" and install packages:

add PythonCall
add FFTW

9. clone git repositories:
git clone https://github.com/p-vychik/stack_gatherer
git clone https://github.com/Marc-3d/multi_quickPIV

10. edit path in get_avg_speed_quickpiv.jl (located in stack_gatherer directory):
include("here should be full path to multi_quickPIV/src/multi_quickPIV.jl")

11. run stack_gatherer with appropriate paths
(pivenv) C:>python C:\Users\..\Documents\stack_gatherer\stack_gatherer.py -i C:\Users\..\Documents\test_input -o C:\Users\..\Documents\test_output -a Z,X,y -f 6 --temp_dir C:\Users\..\Documents\temp --pivjl C:\Users\..\Documents\stack_gatherer\get_avg_speed_quickpiv.jl