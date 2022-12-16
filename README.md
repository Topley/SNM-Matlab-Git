# SNM-lab-software
## Description

This repo houses the primary functions and scripts used in the SNM lab. 
Currently, the main branch is based on the Coco and TRD protocols. The goal is to have a main codebase to take from when starting new projects.
As more projects are developing new techniques, we will add to the repo. 


## Structure

The folders are divided into general topics in which the functions were originally developed. 
Some functions are used in multiple analyses, while others like teh startup file is run only once during or after opening MatLab

### startup.m
This file can be run after opening MatLab, but should be run before running any other scripts. 
This file will define some default settings, like making tableau the default color scheme, add the SNM lab code directory to the path, 
and let you pick which subject folder to start in
It is highly suggested you pick the outermost folder with all of the subject data in it as most of the code is designed to run from this location
while also having flexibility in the number of subjects/files being processed

### Basic functions
This folder contains basic functions like producing a structure of ISIs for each spike train. 
These are functions used in almost every script and analysis

### Coco Code
This folder contains functions and scripts that will get, analyze, and plot auxillary and EMG signals from multiple muscles 

### Coherence Analysis
This folder contains functions specifically for our coherence analyses. This includes intramuscular and intermuscular coherence,
comparing CSTs, batch processing, and producing specific figures

### Delta F Analysis
Same as above, but for delta F analyses

### Event Triggered Analysis
Same as previous anlayses, but based on events depending on the project

### Kinematics
These are functions for kinematic analyses. This includes reading c3d and trb files, creating joint rotation matrices, 
calculating euler and axis-angle joint angles, and saving data into .m files for analyses in MatLab

### MU Matching 
This folder is not yet started, but will contain files to process and group motor units that match across contractions

### Plotting functions
self explanatory

### Preprocessing 
This folder has software to preprocess OTB and c3d files for the current setups the lab has operating. 
Also within this folder are plotting scripts to preview files and the file that runs the cleaning software. 

IMPORTANT - the actual cleaning software is not in this repo and will remain private as it needs authorization for someone to access it
