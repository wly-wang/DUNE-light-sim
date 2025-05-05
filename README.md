Step 0: initialise software environment

Step 1: generate samples

Step 2: merge samples

Step 3: fit GH functions

Step 4: validate


Based on https://github.com/Li-Jiaoyang97/DUNE-light-sim

Relevant tutorial https://cdcvs.fnal.gov/redmine/projects/sbn-analysis-group/wiki/Tutorial_3_Semi-Analytic_mode_How_to_generate_the_correction_curves

#-----------------------------------------------------------------#
# The above serves as the basis structure of how the code works,  #
# below will be a more comprehensive guide on how the code works. #
#-----------------------------------------------------------------#

N.B.: Something to notice is that many names are named protoDUNE-XX, or lateral-xx, many of the names DOES NOT represent ANY specific geometry properties being hard-coded in, the names are just never changed for code consistency, each content and piece of output WILL be explained inside this doc.

FOR ANY QUESTIONS REGARDING THE RUNNING OF THE CODE, please don't hesistate to email me: L. William Wang (L.Wang-121@sms.ed.ac.uk) (U. of Edinburgh), or message me on Slack.

#---------------#
# Custom DUNESW #
#---------------#

Custom DUNESW needs to be set up in case of the need for a custom detector geometry .gdml file, (if custom geometry modification is made in dunecore

For more information, refer to "setting up your workspace" at the DUNE LArSoft Tutorial: https://indico.cern.ch/event/1461779/contributions/6319609/

#-------------------#
# Sample Generation #
#-------------------#

Currently, the code is adapted to run on local Edinburgh PPE computer that utilizes 30 cores,the plan is to integrate the current code onto Larger EDDIE cluster and also the DUNE computing grid.

$$$ Local Edinburgh Computer Method

- TWO SAMPLE GENERATION FOLDER: ONE FOR ARGON, ANOTHER FOR XENON-DOPED ARGON
	-This is for the ease of sample generation, and to avoid minimal large scale editing during one's usage of the code.

- The editing of the fhicl file used to run full optical simulation:
	- Line 1-14: Relevant include file, NO NEED for change in the service files, all possible ones are included for both DUNEFD and ProtoDUNE
	- Line 26-27: This needs to be edited accordingly depending on the relevant geometry, i.e. whether the study is conducted on ProtoDUNE or DUNE-FD.
	- Line 32-39: Relevant geometry .gdml settings:
		- Line 35: Enable only when using STANDARD/NON-CUSTOM geometry. e.g.: if using standard dune10kt fd-hd geometry, set as "dune10kt_v6"
		- Line 38-39: Enable only if when using NON-STANDARD/CUSTOM geometry. i.e.: if using STANDARD geometry, AuxDet geometry is automatically set in the service files.
	- Line 80-94: Analyzers:
		- Line 89: MUST SET TO TRUE, this is THE output tree that is used by the whole code structure.
		- Line 90-92: These can be set to false, but set to TRUE if one would like to obtain more information about specific OpDet or detailed info on detected hits from LArAna SimPhotonCounter.
	- Line 101: "stream1: [ out1 ]" -> Will produce artroot files that can be utilized for, e.g. eventdumping
		    "stream1: [  ]" -> Will NOT produce artroot files.

- Fhicl files shall be the SAME for both ARGON & XENON-DOPED ARGON.

- Python script used to Start the generation:
	- ARGON VERSION:
		- Line 15-17: p -> value of energy for argon shall be 9.69. n -> number of photons emitted per photon origin (later would be looped 100 times, that would be the total number of photon emitted per photon origin), i.e. 100,000 means 10M photons per photon origin. 
	- XENON-DOPED ARGON VERSION:
		- Line 15-17: p -> value of energy for xenon-doped argon shall be 7.13, n remains the same, same logic.

	-Common b/w ARGON and XENON-DOPED VERSION:
		- Line 28-35: x,y,z coordinate for the simulated measured photon origins. Generally enough to just simulate a quarter of the detector, like the values currently set, dense closer to edge of drift direction of detector, less dense in the middle, span across the whole drift direction; y-value span half of the height of detector, and z value span avg across half of beam-direction of detector. TOTAL NUMBER OF MEASURED PHOTON ORIGINS = x_vals * y_vals * z_vals. Recommendation, DON'T modify too much of x and y other than make sure to match the edge/border dimension of x and y for the geometry one's investigating. Change only the z-direction coordinate depending on how many points one wants to generation.

		- Line 60: Number of cores to use, 30 max for Edinburgh PPE computing server.
	
	- RUNNING THE GENERATION:
		- ```python multiprocess.py```

#----------------#
# Sample Merging #
#----------------#

After the generation stage, there will be xx number of folder named "process xx", 30 cores means 30 processes running at the same time, so 30 process folders. The goal of the merging step is to merge all of these output root file into 1 master analysis .root file. 

First thing to do before editing: DELETE EXISTING .txt files: 'List_of_files.txt' & 'protodune_optical_mapping.txt' & 'output.root'

- makeChannelMap.sh:
	- This script obtains and creates a list of optical detectors in the geometry and their coordinate, information obtained from the LArSoft output log that is stored at "/generate/processXX/log". 
	- THIS IS THE FIRST THING TO RUN DURING MERGE STEP. 

- makeFileList.sh:
	- This script is used to create a text list of all of the root files that needs to be merged. 
	- THIS NEEDS TO BE RAN NEXT AFTER makeChannelMap.sh
	- If one needs additional points (or more photon per photon origin), instead of regenerating, one can generate extra amount and then just merge onto the existing generated files. Through modifying the makeFileList.sh. 

- src/main.cpp:
	- This is the main script for merging. 
	- Line 24: This should match as the number of photon per measured photon origin value used in the generation stage. Again, total number of photon per measured photon origin value would be '100*genPhotons', i.e. genPhotons=100,000 -> total photon per point=10M. 
	- Line 52: Name of the master output analysis .root file. 
	- Line 121: The value AFTER 100* needs to be changed the same way as Line 24.

- RUNNING OF THE MERGING STEP:
	- ```source makeFileList.sh```
	- ```source makeChannelMap.sh```
	- ```make```
	- ```bin/FastOpticalMerge List_of_files.txt```

#------------------#
# Fit GH Functions #
#------------------#

Two separate versions for the Fit part of the code: 1 for Argon-only and another for Xenon-doped Argon. This part of the code is the main code part that is used to fit the GH correction curve onto our "ideal/pure geometrical" photon detector response to account for the effect of Rayleigh Scattering. 

- makeChannelMap.sh (same as in the merge step of the code):
	- This script obtains and creates a list of optical detectors in the geometry and their coordinate, information obtained from the LArSoft output log that is stored at "/generate/processXX/log". 
	- THIS IS THE FIRST THING TO RUN DURING THE FIT STEP. 

-src/main.cpp:
	- This is the main script for the fit part of the code
	- Line 40: min_number_entries, set for the minimum value of a hit, this is mainly used when the number of generated photon per photon origin is not a large enough. 50-100 is about right, no need to tune unless analysis specific on requiring specific hit values. 
	- Line 47: Self-explanatory, the center of the YZ plane coordinate of the detector. CRUCIAL PARAMETER TO CHANGE, DON'T FORGET!!
	- Line 57: The absorption length, 20m for Argon, while it's 80m for xenon-doped argon. 
	- Line 131: The commented out line in such format is a line of code that can be utilized (if needed) to look at only specific photon detector responses. As long as one understands fully the indexing and where the interested photon detectors are. 
	- Line 285: Leave this to true unless modification made to photon detectors. This code is now adapted/optimized to run for mainly just DUNE, with the photon detector type being always set to X-ARAPUCAs. The same code is also utilized by SBND, for example, but adding in a bunch of conditions confuses the user and is best to strictly just look at DUNE-related geometries. 
	- Line 288: Number of offset angled bin, used to look at the angled space-points coming to the photon detector, i.e. 0 deg means directly seen by the photon detector, hitting straight on, while 90 deg means it's hitting on the side perpendicular to the photon detector. 
	- Line 297-321: This snippet of code is very crucial to this part of the code.
		Line 298-299: These sets the range of the distance d to the photon detector. 
		Line 


