"""
RegressPython.py: A Python script for preprocessing functional MRI (fMRI) data.

This script performs several key preprocessing steps on fMRI data, specifically targeting the
cleaning and processing of BOLD (Blood-Oxygen-Level-Dependent) signals. It utilizes libraries
like nilearn, nibabel, numpy, and pandas to handle neuroimaging data and perform statistical
operations.

Key Features:
1. Memory Configuration: Adjusts the Python environment to utilize specified memory resources.
2. Input and Output Handling: Accepts file locations for BOLD signals, output, and confounds as command-line arguments.
3. Confound Processing: Reads and preprocesses confound data from a specified file, handling missing values.
4. BOLD Signal Cleaning: Performs cleaning of BOLD signals using specified confounds, applying low pass and high pass filters.

The script expects command-line arguments in the following order:
1. Location of the BOLD signal file.
2. Desired output location for the processed file.
3. Location of the confound file.
4. Repetition time (TR) as a floating-point number.
5. Memory allocation setting for Python.

Dependencies:
- nilearn: For neuroimaging data manipulation and analysis.
- nibabel: For reading and writing neuroimaging data files.
- numpy: For numerical operations on arrays.
- pandas: For handling tabular data.

Note:
This script is part of a larger pipeline and is typically invoked via a shell script, passing the required arguments.
"""

# Perform imports.
import nilearn.image
import nibabel as nib
import numpy as np
import pandas as pd
import sys
import os

# Make sure that python can access all of the memory.
os.environ['PYTHONMEM'] = sys.argv[5]

# Get the input, output and confound location, passed as an argument in the .sh script.
bold_location = sys.argv[1]
output_location = sys.argv[2]
confound_file = sys.argv[3]
tr = float(sys.argv[4])

# Get the confound_file and put into a dataframe.
confound_df = pd.read_csv(confound_file, sep="\t")

# Remove nans in derivatives by copying the first non-nan.
confound_df = confound_df.fillna(confound_df.shift(-1))

# Specify the confounds to use.
confounds_to_use = ["c_comp_cor_00", "c_comp_cor_01", "c_comp_cor_02", "c_comp_cor_03", "c_comp_cor_04",
                    "w_comp_cor_00", "w_comp_cor_01", "w_comp_cor_02", "w_comp_cor_03", "w_comp_cor_04",
                    "trans_x", "trans_x_derivative1",
                    "trans_y", "trans_y_derivative1",
                    "trans_z", "trans_z_derivative1",
                    "rot_x", "rot_x_derivative1",
                    "rot_y", "rot_y_derivative1",
                    "rot_z", "rot_z_derivative1",
                    "non_steady_state_outlier00"
                   ]

# Get the bold signal as a nifti image in nibabel.
#bold_signal = nib.load(bold_location)

# Turn this into a numpy array, transpose it for cleaning.
#bold_array = bold_signal.get_fdata()
#transposed_bold_array = np.transpose(bold_array, (3,0,1,2))

# Calculate some numbers for reshaping.
#n_samples = transposed_bold_array.shape[0]
#n_voxels = transposed_bold_array.shape[1] * transposed_bold_array.shape[2]\
#	* transposed_bold_array.shape[3]

# Make the 2D version of the signal.
#bold_2d_transposed = transposed_bold_array.reshape(n_samples, n_voxels)

# Clean the 2D Bold signals.
cleaned_image = nilearn.image.clean_img(bold_location, confounds = confound_df.loc[:,confound_df.columns.isin(confounds_to_use)],\
	 low_pass = 0.05, high_pass = 0.0025, t_r = tr, standardize = False)

# Return the data to 4D.
# cleaned_bold_4d_transposed = cleaned_bold_2d_transposed.reshape(transposed_bold_array.shape)
#
# Move the axes around again.
# cleaned_bold_4d = np.transpose(cleaned_bold_4d_transposed, (1,2,3,0))
#
# Output the image.


#cleaned_image = nib.Nifti1Image(cleaned_bold_2d_transposed, bold_signal.affine, bold_signal.header) #4d was changed to 2d
nib.save(cleaned_image, output_location)
