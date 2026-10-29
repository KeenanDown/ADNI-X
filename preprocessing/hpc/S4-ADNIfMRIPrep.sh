#!/bin/bash
# Historical SGE example; local paths and notification address were removed.
# Read ../README.md and adapt paths before use; processing may overwrite files.
#$ -pe smp 1
#$ -l h_vmem=128G #You might need to change this a few times for redoing subjects that broke. Default = 8G
#$ -l h_rt=40:0:0
#$ -l centos
#$ -t 1-286

# Specify the project-run code
PROJCODE=20Nov23

# Round modifier
ROUNDMOD=-Round2

# Specify the subject list location.
SUBJECT_LIST=/path/to/project/Datasets/ADNI/$PROJCODE/SubjectsToProcess_Round2.txt          # Use this to specify a list of subjects.

# Specify the subject folder inside of ADNIBIDS.
SUBJECT=$(sed -n "${SGE_TASK_ID}p" $SUBJECT_LIST)

# Create a working directory for the subject.
workdir=/path/to/scratch/$USER/fMRIPrepWork/$SUBJECT/
mkdir -p $workdir

# Specify the fMRIPrep .sif location.
fMRIPrepLOCATION=/path/to/project/Software/fMRIPrep-24-1-1.sif

# Specify the input location.
INPUT_LOCATION=/path/to/project/Datasets/ADNI/$PROJCODE/ADNIBIDS-$PROJCODE

# Specify the output location.
OUTPUT_LOCATION=/path/to/project/Datasets/ADNI/20Nov23/S4Out-$PROJCODE

# Make sure this exists
mkdir -p $OUTPUT_LOCATION

# Assign the Freesurfer license file
FS_LICENSE=/path/to/project/Resources/license.txt

# Quickly grab the amount of memory available in MB.
	# Retrieve the memory limit in GB
	MEMORY_LIMIT_GB=$(qstat -j $JOB_ID | grep -oE 'h_vmem=[^,]+' | cut -d'=' -f2)

	# Extract the numeric part of the string
	GB_NUMBER="${MEMORY_LIMIT_GB%"G"}"

	# Perform the conversion and multiplication
	MB_AVAILABLE=$((GB_NUMBER * 1024))

# Do the processing for SUBJECT.
if ! [ -d "$OUTPUT_LOCATION/$SUBJECT" ]
then
	apptainer run --cleanenv $fMRIPrepLOCATION $INPUT_LOCATION $OUTPUT_LOCATION \
	-w $workdir --clean-workdir --n_cpus ${NSLOTS} --omp-nthreads 1 participant \
	--participant-label $SUBJECT --mem $MB_AVAILABLE --fs-license-file $FS_LICENSE \
	--output-spaces MNI152Lin:res-2 --skip_bids_validation
else
	echo Subject $SUBJECT appears to already have been processed.
fi

# Destroy the data afterwards, because fMRIPrep sucks at it.
rm -r $workdir
