Step 0: initialise software environment

Step 1: generate samples

Step 2: merge samples

Step 3: fit GH functions

Step 4: validate


Based on https://github.com/Li-Jiaoyang97/DUNE-light-sim

Relevant tutorial https://cdcvs.fnal.gov/redmine/projects/sbn-analysis-group/wiki/Tutorial_3_Semi-Analytic_mode_How_to_generate_the_correction_curves

#--------------------------------------------------------------#
The above serves as the basis structure of how the code works,
below will be a more comprehensive guide on how the code works.
#--------------------------------------------------------------#

#---------------#
# Custom DUNESW #
#---------------#

Custom DUNESW needs to be set up in case of the need for a custom detector geometry .gdml file, (if custom geometry modification is made in dunecore

For more information, refer to "setting up your workspace" at the DUNE LArSoft Tutorial: https://indico.cern.ch/event/1461779/contributions/6319609/

#-------------------#
# Sample Generation #
#-------------------#
