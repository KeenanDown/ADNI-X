#!/bin/bash
# Historical SGE example; local paths and notification address were removed.
# Read ../README.md and adapt paths before use; processing may overwrite files.
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=10:0:0
#$ -l centos
#$ -t 1-226 # This needs to be the number of subjects.

# Specify the project-run code
PROJCODE=Nov20

# Specify the fMRIPrepOut location
FMRIOUTLOCATION=/path/to/project/Datasets/ADNI/20Nov23/S4Out-20Nov23

# Specify where to save the list of subjects
SUBJECT_LIST=/path/to/project/Datasets/ADNI/20Nov23/S5Redo1.txt

# If the subject list doesn't exist, then make it.
if [ ! -f "$SUBJECT_LIST" ] && [ "${SGE_TASK_ID}" -eq 1 ]; then
    cd "$FMRIOUTLOCATION"
    for dir in ./*/; do
        basename "$dir";
    done | grep -v 'logs' | grep -v 'sourcedata' > "$SUBJECT_LIST"
fi

# If not the first job, wait until it exists.
if [ ${SGE_TASK_ID} != "1" ]; then
    while [ ! -f "$SUBJECT_LIST" ]; do
        echo "Waiting for $SUBJECT_LIST to be created by the first job..."
        sleep 10 # Waits for 10 seconds before checking again
    done
    echo "$SUBJECT_LIST has been detected. Proceeding with the job..."
fi

# Load anaconda and activate the NiPy environment.
module load anaconda3
conda activate NiPy

# Provide the location of the RegressFilter.py file.
PYTHON_LOCATION=/path/to/project/Software/ADNI-PID/RegressFilter.py

# Select the current subject from that list.
SUBJECT=$(sed -n "${SGE_TASK_ID}p" $SUBJECT_LIST)

# Check if the SUBJECT variable is empty
if [ -z $SUBJECT ]; then
    echo "Error: Unable to find subject at position ${SGE_TASK_ID} in $SUBJECT_LIST"
    exit 1
fi

# Specify the input BIDS.
INPUT=/path/to/project/Datasets/ADNI/20Nov23/S4Out-20Nov23
OUTPUT=/path/to/project/Datasets/ADNI/20Nov23/S5Out-20Nov23

# Copy the entire old directory in only the first job.
if [ ${SGE_TASK_ID} == 1 ]; then
    if [ ! -f "$OUTPUT/CheckCopied-Cleaned.temp" ]; then
        cp -r $INPUT $OUTPUT
        echo Copy Complete! >> $OUTPUT/CheckCopied-Cleaned.temp
    fi
elif [ ${SGE_TASK_ID} != 1 ]; then
    while ! [ -f "$OUTPUT/CheckCopied-Cleaned.temp" ]; do
        sleep 10
    done
fi

# Name a working directory for saving file names. Make it if it doesn't exist.
FILEDIRECTORY=/path/to/scratch/$USER/$PROJCODE-FileNames/$SUBJECT

if ! [ -d "$FILEDIRECTORY" ]
then
	mkdir -p $FILEDIRECTORY
fi

# Create a list of all of the preprocessed bold files inside of the SUBJECT's directory.
find $OUTPUT/$SUBJECT -type f -name "*$SUBJECT*desc-preproc_bold.nii.gz" > $FILEDIRECTORY/niigzLocations.temp

# Create a list of all of the confound locations.
find $OUTPUT/$SUBJECT -type f -name "*$SUBJECT*desc-confounds_timeseries.tsv" > $FILEDIRECTORY/confoundLocations.temp

# Create a list of all of the .json locations.
find $OUTPUT/$SUBJECT -type f -name "*$SUBJECT*desc-preproc_bold.json" > $FILEDIRECTORY/jsonLocations.temp

# Create a file with a list of the new file names.
cat $FILEDIRECTORY/niigzLocations.temp | sed 's/\(.*\)desc-preproc_bold.nii.gz/\1desc-preproc-cleaned_bold.nii.gz/'\
    > $FILEDIRECTORY/cleanedniigzLocations.temp

# Get the total number of BOLD files in niigzLocations
FILE_COUNT=$(wc -l < "$FILEDIRECTORY/niigzLocations.temp")

# Loop over each file in the list.
for i in $(seq 1 "$FILE_COUNT"); do

    # Specify the input location, output location, and location of the confounds file
    BOLD_LOCATION=$(sed -n "${i}p" $FILEDIRECTORY/niigzLocations.temp)
    # The following line has been changed to overwrite _preproc_bold, rather than obey standard BIDS rules.
    OUTPUT_LOCATION=$(sed -n "${i}p" $FILEDIRECTORY/niigzLocations.temp) # This is changed in this version.
    CONFOUNDS_FILE=$(sed -n "${i}p" $FILEDIRECTORY/confoundLocations.temp)
    JSON_LOCATION=$(sed -n "${i}p" $FILEDIRECTORY/jsonLocations.temp)

    # Specify the memory (set this the same as above).
    PYMEM=$(qstat -j $JOB_ID | grep -oE 'h_vmem=[^,]+' | cut -d'=' -f2)B

    # Read the JSON file as a string
    JSON_DATA=$(cat "$JSON_LOCATION")

    # Extract the RepetitionTime value using grep and sed
    TR=$(grep -o '"RepetitionTime": *[0-9.]*' <<< "$JSON_DATA" | sed 's/"RepetitionTime": *\([^ ]*\)/\1/')

    python $PYTHON_LOCATION "$BOLD_LOCATION" "$OUTPUT_LOCATION" "$CONFOUNDS_FILE" $TR "$PYMEM"

done

# Remove the leftover file name data.
#rm -r $FILEDIRECTORY
