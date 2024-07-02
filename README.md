
# stack_gatherer.py

`stack_gatherer.py` is an auxiliary Python 3 script designed to process raw image planes derived from a light-sheet fluorescence microscope. It provides functionalities such as:

	•  Reading and processing image files in TIFF and bitmap format.

	•  Live visualization of Z, X, Y projections in Napari viewer and saving the projection images to the user-defined directory.

	•  Calculating cell average migration speed with Julia script [QuickPIV](https://github.com/Marc-3d/multi_quickPIV), live plotting of the speed data in a matplotlib widget window, and reporting peaks presence on the graph with a Telegram message or socket communication.

	•  Interacting with the microscope controller software.

## Installation
To install and set up the stack_gatherer environment, you will need to use the terminal 
or command prompt on your computer:

#### stack_gatherer.py was tested in Python 3.10 environment for Windows 10 only!
    
#### Using Conda:
1. Open your terminal

	•  On Windows, you can search for "Command Prompt" or "Anaconda Prompt" if you have Anaconda installed

	•  On macOS or Linux, you can open the "Terminal" application

1. Create a new Conda environment

	•  Run the following command to create a new environment named stack_gatherer with Python version 3.10
	```bash
	conda create -n stack_gatherer python=3.10
	```
1. Activate the new environment

	•  Before installing packages, activate the new environment using:
	```bash
	conda activate stack_gatherer
	```
1. Install required packages

	•  Ensure you have a requirements.txt file in the current directory. This file contains a list of packages needed for stack_gatherer.py

	•  Install the packages with the following command:
	```bash
	conda install --yes --file requirements.txt
	```
#### Using Pip (without Conda):

If you prefer not to use Conda or do not have it installed, you can use pip, the Python package installer.
1. Ensure you have Python 3.10 installed

	•  You can check your Python version by running in your terminal:
	```
	python --version
	```

	•  If you don't have Python 3.10, download and install it using the guideline for your Operating System.

1. Install packages with pip

	•  Run the following command to install the required packages listed in requirements.txt:

	```bash
	pip install -r requirements.txt
	```

	•  This will install the packages in your current Python environment.

Remember to run these commands in the directory where your requirements.txt file is located, or provide the path to it.
For example, if your requirements.txt is in the Downloads folder, 
you would use --file ~/Downloads/requirements.txt for Conda or -r ~/Downloads/requirements.txt for pip.

## Functionality
##### Watch Directory: 
The script continuously observes the input directory for any new files.

##### Collect and Process Images: 
When new image files (.tif, .tiff, or .bmp) are added, they are collected, processed, and organized into appropriate stacks.

##### Create Projections: 
Different projections are made from the processed image data (X, Y, Z) based on the specified axes.

##### Calculation of Average Migration Speed:
* Note! Additionally, a Julia environment should be installed and configured to use this functionality. Refer to the section "Special Usage" – "Julia Install and Environment Configuration".
* Utilizes the quickPIV auxiliary Julia script to calculate the average migration speed of specimen cells and plot the speed graphs in separate matplotlib widgets windows. 

##### Napari Visualization: 
Opens a Napari viewer to display the projections.

##### Telegram Alerts: 
Sends Telegram messages for specific events, such as the detection of migration peaks in the images.

## Technical Description of the Workflow
##### Watching Input Directory for Image Planes (watchfiles module):

•  Monitors an input folder for a JSON configuration file for the timelapse which defines the number of image planes in the stack, specimen quantity on the stage, illumination modes – the specified parameters are used for processing and visualization of raw image data.

•  Arrival of image files with a correct file extension (.bmp, .tiff, .tif) triggers a check if the file is already added to the image stack with the method **check_stack_and_collect_if_ready()**. The check is based on a file signature derived from the file name mask:
```
f"{dataset_name}timelapseID-{self.timelapse_id}_SPC-{self.specimen}_TP-{self.time_point}_ILL-{self.illumination}_CAM-{self.camera}_CH-{self.channel}_PL-{self.plane}-outOf-{self.total_num_planes}{additional_info}.{self.extension}"
```

•  The signature of the image file is constructed as a tuple of values, extracted from the appropriate fields defined by the file name mask.

##### Collecting and Processing Image Files
•  After an image successfully passes the signature check, the method **check_stack_and_collect_if_ready()** compares the plane number with the expected number of images in the stack defined by the lapse parameter total_num_planes. If the microscope has supplied all planes for processing, the method **collect_files_to_one_stack_get_axial_projections()** initiates the consecutive generation of projections from stack planes.

##### Calculating Average Speed
* Calculation of the cells migration is based on the separate Julia script (quickPIV). The implementation within the `stack_gatherer.py` workflow requires input of two Z-projections of the same specimen from consecutive time points. The source of the projections is a shared_dict object which keeps several shared multiprocessing queues accessible by specimen index, the total number of queues is defined by the number of specimens provided with the timelapse configuration in JSON file.
* Current data supply and execution of Julia code is dependent on the juliacall package. Its main limitation consists of inability to be called from python threads if it was imported on the global scope and not being able to be imported multiple times in the local scope. Both issues define the usage of juliacall in an infinite loop within separate python process.
* After signal defining the end of the timelapse is received (currently the role of that signal plays the creation of empty file named `STOP_STACK_GATHERING` within input directory), average speed values for every specimen and framepoint are saved within output folder as csv file  named `quickPIV_data_{%H_%M_%S}.csv`, where H, M, S – are hour, minutes, and seconds of time file creation

##### Telegram Alerts and Socket Communications
* `stack_gatherer.py` is able to inform the user with Telegram messages for such events, as the migration peaks in the images of specimens. The message format looks like “Detected on {frame_number}, specimen {specimen_number}! 
* To activate notifications in Telegram one has to provide as an argument with a key --bot_config a path to the separate config file containing a one line the bot token value and the chat id, separated with space:
`bot_token chat_id`

* One can get bot token from the developers, to get chat id one has to proceed with the following steps:
  1. find Telegram bot named `stack_gatherer.py` with a search menu
  1. write random message to the bot
  1. execute command in the terminal or PowerShell. JSON response will contain key “chat” and a corresponding dictionary containing keys “id”, “first_name”, “last_name” etc. You need to copy the “id” value for the “chat” key and put it in the file in the above-described format.
```bash
curl https://api.telegram.org/bot{put_bot_token_here_without_curl_brackets}/getUpdates
```


* Socket communication with microscope controller software is supported, one needs to provide a port value of a listener with -p key. General functionality of the socket communication is to inform the controller that script works via sending a “heartbeat” message, and informing about the peak detection


## Usage
* Input Folder (Required), -i, --input:  
Specify the path to the input folder to watch for new image files.
* Output Folder (Required), -o, --output:  
Specify the path to save the generated image stacks and processed data. 
* Axes (Optional), -a, --axes:  
Enumerate without spaces the axes to project (e.g., --axes X,Y,Z) during the processing of the image data.  
* Anisotropy Factor (Required if –axes key was provided), -f, --factor_anisotropy:  
Provide a value for correcting the anisotropic distortions of the projections. 
* Process Z-projections (Optional), -process_z_projections:  
Used only in special cases outside of the timelapse to calculate and plot average migration speed from available z-projection images. Requires –pivjl,  --specimen_quantity  
* Specimen Quantity (Required if  -process_z_projections was set), --specimen_quantity:  
Used only in special cases outside of the timelapse. Specify the number of specimens for processing already collected Z-projections and creating plotting windows.  
* Temporary Directory (Optional), --temp_dir:  
Define a directory to store input images with incorrect filenames before processing.  
* Path to run file for quickPIV call, --pivjl (Required if  -process_z_projections was set):  
Can be used either during timelapse to calculate average migration speed in a live mode, or in -process_z_projections mode  
* Debug Run (Optional), -debug_run:  
Use this flag to clear the output folder before loading new files for testing purposes.  
* Path to Telegram bot configuration file  (Optional), --bot_config:  
The configuration file provides credentials to receive messages informing about average speed peaks in the processed images, can be used both in live timelapses or to process the ready z-projections.  
* Process Z-Projections (Optional):  
Activate this option to process the input files as Z-projections for average speed calculations with quickPIV.

## General Usage Example:
 
```bash
python stack_gatherer.py -i path/to/input -o path/to/output -a X,Y,Z -f 1 --temp_dir /temp  --pivjl -p 5555
```
1. -i (input path) here is a folder where microscope controller saves the raw planes and timelapse configuration file
1. -a key sets projections to display in napari viewer
1. -f key gives an anisotropy factor value - the coefficient to correct distortions of X and Y projections, in future it will be calculated based on the timelapse configuration properties
1. --temp_dir path to the directory where all files not passing the name mask filter will be moved
1. path to the `get_avg_speed_quickpiv.jl` script which runs the QuickPIV package for calculation of the average migration speed of the cells.
1. -p defines the port number to communicate with microscope controller software

## Special Usage
Independent processing of Z-projections is supported. As an output, a CSV file containing the average speed data is created, also matplotlib plots are saved as PNG files displaying the change of average speed for every individual specimen marking detected peaks with a red dot.

##### Julia Install and Environment Configuration
1. Install Anaconda for the local user:Anaconda2. Create a separate environment (tested with Python 3.10 and Julia 1.10):
```bash
conda create -n pivenv python=3.10
```
2. Activate the created environment:
```bash
conda activate pivenv
```
3. Install Julia manager:
```bash
conda install -c conda-forge juliaup
```
4. Choose Julia version:
```bash
juliaup add 1.10~x64
```
5. Install Python packages:
```bash
pip install watchdog matplotlib juliacall napari PyQt5
```
6. Run Julia in the created environment, in REPL mode set the Python path:
```bash
ENV["PYTHON"] = "C:\\Users\\..\\anaconda3\\envs\\pivenv\\python.exe"
ENV["JULIA_CONDAPKG_BACKEND"] = "Current"
```
7. Run package manager by pressing the key ] and install packages:
```bash
add PythonCall add FFTW
```
8. Clone the git repository:
```bash
git clone https://github.com/Marc-3d/multi_quickPIV
```
9. Edit the path in `get_avg_speed_quickpiv.jl` (located in the stack_gatherer directory):
```
include("here should be full path to multi_quickPIV/src/multi_quickPIV.jl")
```
