from genericpath import exists
import nibabel as nib
import numpy as np
from regfusion import vol_to_fsaverage
import nibabel.processing
import os.path
from os import path

input_path = '/usr/local/freesurfer/7.2.0/subjects/fsaverage/mri/T1.mgz'
output_path = '/home/senthilp/caesar/wen/degree/T1_2mm.mgz'

if not path.exists(output_path):
    voxel_size = [2, 2, 2]
    input_img = nib.load(input_path)
    resampled_img = nibabel.processing.resample_to_output(input_img, voxel_size)
    nib.save(resampled_img, output_path)


# Shen registered atlas file name
atlas = '/home/senthilp/caesar/camcan/cc700/freesurfer_output/sub-CC110098/mri/shen_freesurfer.mgz'

# 268 array values
corr_values = np.random.uniform(0,1,268)

# Load shen atlas data
obj = nib.load(atlas)
data = nib.load(atlas).get_fdata()

# Uniques atlas ID's from 1 to 268
want = np.arange(1,269)

# Create new matrix
map_data = np.zeros((data.shape))

# Loop through the want ID's and copy corr values to the matrix
for i, id in enumerate(want):
    index = np.where(data == id)
    map_data[index] = corr_values[i]


# Create nifti brain image volume
output = 'mapped_array.nii.gz'
t1 = nib.load(output_path)
affine = t1.affine
hdr = t1.header
result_img = nib.Nifti1Image(map_data, affine, header=hdr)
result_img.to_filename(output)

# Map volume to surface
input_file = output
lh, rh = vol_to_fsaverage(input_file, 'mapped_vol/sub-CC110098')
