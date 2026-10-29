#!/bin/bash
# Historical SGE example; local paths and notification address were removed.
# Read ../README.md and adapt paths before use; processing may overwrite files.
#$ -pe smp 1
#$ -cwd
#$ -o ~
#$ -j y
#$ -l rocky
#$ -l h_vmem=16G
#$ -l h_rt=240:0:0
#$ -t 1-1154

# This will load each fMRI session and create a MATLAB PhiID tensor (116x116x16)
# and save it as a .mat in the same directory.
# Note that it also does the parcellation, so we *don't* need to run S8a first.
# S8bx saves files in place, whereas S8b makes a new folder.

# Specify the project-run code
PROJCODE=20Nov23

# Specify the subject list. This list will have sub- before each subject.
SUBJECT_LIST=/path/to/project/Datasets/ADNI/20Nov23/S8bx-116-Redo0.txt

# Specify the directory for ELPh.
ELPH_DIRECTORY=/path/to/project/Software/elph

# Specify the number of parcels as either 116 or 232
PARCEL_NUM=116

# Specify the location of the parcellation.
PARCELLATION_LOCATION=/path/to/project/Software/ADNI-PID/Parcellations/Schaefer2018_${PARCEL_NUM}Parcels_7Networks_order_FSLMNI152_2mm.nii

# Specify the parcellation name
PARCELLATION_NAME=ExtSchaefer$PARCEL_NUM

# Specify the subject name by looking up in the SUBJECT_LIST.
SUBJECT=$(sed -n "${SGE_TASK_ID}p" $SUBJECT_LIST)

# Make a new directory for the files with tensors.
OUTPUT=/path/to/project/Datasets/ADNI/20Nov23/S8bOut-$PROJCODE

# Create a file with a list of the BOLD run locations.
find $OUTPUT -type f -name "*$SUBJECT*preproc_bold.nii*" > $OUTPUT/bold_niiLocations_$SUBJECT

# Load MATLAB ready for processing.
module load matlab

# Do the matlab script to each of them.
while read BOLD_FILE; do
	matlab -nodisplay -nosplash -r "create_phiID_tensors_v3('$BOLD_FILE', '$ELPH_DIRECTORY', '$PARCELLATION_LOCATION', '$PARCELLATION_NAME'); quit"
done < $OUTPUT/bold_niiLocations_$SUBJECT

# Destroy the file bold_gzLocations
rm $OUTPUT/bold_niiLocations_$SUBJECT
