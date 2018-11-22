changing the file

# preprocessingfNIRS
UCLA SCN lab preprocessing function for fNIRS data, collected with NIRx 

# About 
This is a pretty simple function (more of a functionalized script, really) that automates the SCN lab's current preprocessing pipeline for fNIRS data, collected with NIRx NirScout and NirSport systems. Currently, it only comes in Matlab flavor right now (written by Shannon Burns, 2018, MIT License). Python version is in development. See below for more details on each step of the pipeline, and how to use this code yourself. (Note: only lightly maintained by one person, so don't be surprised to find bugs right now. Please report any you see!)

# Dependencies
preprocessingfNIRS depends on two publicly available Matlab packages - inpaint_nans and Homer2. The first is small and included with the preprocessingfNIRS repo. The second is a larger collection of various fNIRS processing functions and can be downloaded at https://www.nitrc.org/projects/homer2.

# Contents
- preprocessingfNIRS is structured as 4-step pipeline: extracting raw data from the NIRx files; removing bad channels; creating auxillary variables used by the Homer2 .nirs data format; preprocessing the data with our chosen algorithms from Homer2. Please see comments in main .m file to read about specific functionality. The preprocessingfNIRS.m file has almost everything you need written in it as a main function and subfunctions. Exceptions are the dependencies mentioned above. 

- The inpain_nans dependency. A small helper function that interpolates missing values in a timecourse. 

- The "testQCoD.m" helper function. This tells you which channels are being marked as good and bad, and also visualizes the power spectral density of the raw signal in each fNIRS channel. Doing so is our current approach to finding bad channels (good channels have a preponderance of signal in slow frequencies, pure noise has random amounts of signal at every frequency). We're using a modified version of the quartile coefficient of dispersion (https://en.wikipedia.org/wiki/Quartile_coefficient_of_dispersion), abbreviated QCoD, to automatically decide which channels have good or bad signal. Essentially, it sums the frequency amplitudes in the first and third quartiles of the frequency range, and then compares them via (Q1-Q3)/(Q1+Q3). Larger QCoD is cleaner signal. Default threshold is set to 0.1 in preprocessingfNIRS. In testQCoD, change this to <0.1 to test allowing greater noise in the signal, or change to >0.1 for more stringency. Use this test function to check if the currently selected QCoD threshold for bad channel marking is what you want. Can be run on raw data, but does not do any other preprocessing. 

- extractNIRxData. Reads the raw nirs data from .wl1 and wl2 files, as well as meta data from .hdr file. 

- removeBadChannels. Creates a mask for which channels to keep or mark as noise (see above entry for testQCoD for removal method).

- getMiscNirsVars. Creates auxillary variables from data that we don't need if just reading a matrix of preprocessed data, but is needed if saving a .nirs file. 

- fNIRSFilterPipeline. Our current filtering pipeline. Add or remove methods from here to change up exactly how the data is preprocessed (this is where the Homer2 dependency comes in). 

# How to Use
- To download: clone or pull repo to your desired local directory. Then add folder and subfolders to your Matlab path via: 
addpath(genpath('[YOUR DIRECTORY]'));

- preprocessingfNIRS arguments: the function takes 3 arguments - raw_dir (folder path with all single-subject or dyad-level raw data files - path must be relative to current matlab path), dataprefix (string. Prefix of every folder name that should be considered a data folder. E.g., MIN for MIN_101, MIN_102, etc.), and dyads (0 or 1. 1 if hyperscanning, 0 if single subject). No output arguments, but saves a .nirs file of raw data and a .mat file of z-scored and non z-scored oxy, deoxy, and totaloxy matrices (timepoint x channel).

- preprocessingfNIRS requirements: besides dependencies installed, data must be in the following structure: MAIN_DIRECTORY/SUBJECT/SCAN/raw files. Can also be MAIN_DIRECTORY/SUBJECT/raw if only one scan per person was collected. If dyads, structure will be MAIN_DIRECTORY/DYAD/SCAN/SUBJ1_OR_2/raw files, or without the scan level if only one scan per dyad was collected. All NIRx files must be together in each raw file folder. All subjects and scans must start with a study-specific prefix.

- testQCoD arguments: testsubjpath (relative path to one scan folder with all raw nirx files), QCoDthresh (QCoD threshold to test out), and suppressPlot (0 or 1. 0 to let plots happen, 1 to keep them from coming up). Output is the channel mask, a vector of 1s and 0s where 1s correspond to good channels and 0s to bad. 

# Known Difficulties
Since this is still being worked on, here are some known issues that might pop up and how to fix them right now:
 - If you get an "Index exceeds matrix dimensions" error in hmrMotionArtifact for a subject that's not the first file: Check the SD Mask structure in the .hdr of that subject to see if it matches the channel structure of the selected probeInfo file. If the wrong probeInfo file was chosen, this will throw the error. also happens if the wrong montage was selected in recording, Simply copy-paste the correct SD Mask and ChannelDistance list into the .hdr file from a subject's .hdr file that had the correct montage.

- Delete all false start files from the data directory, or will cause script to error out. 

# FAQ Not Covered Above
None yet, so ask me things if you're confused! 
