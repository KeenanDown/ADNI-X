#!/bin/bash
# Historical SGE example; local paths and notification address were removed.
# Read ../README.md and adapt paths before use; processing may overwrite files.
#$ -pe smp 1
#$ -m bae
#$ -l h_rt=240:0:0
#$ -l h_vmem=64G

# Specify the project-run code
PROJCODE=6Jan25

# Select the folder which is called 'ADNI' inside of the exported files.
ADNI_INPUT_FOLDER=/path/to/scratch/$USER/$PROJCODE/$PROJCODE/ADNI

# Provide the folder requested from Clinica containing all of the additional study data.
STUDY_DATA_FOLDER=/path/to/project/StudyData-$PROJCODE

# Name the output folder.
OUTPUT_FOLDER=/path/to/scratch/$USER/ADNIBIDS-$PROJCODE



# Load anaconda.
module load anaconda3

# Activate the clinica environment.
conda activate clinicaEnv

# Convert the ADNI data to BIDS.
clinica convert adni-to-bids -np ${NSLOTS} $ADNI_INPUT_FOLDER $STUDY_DATA_FOLDER $OUTPUT_FOLDER
