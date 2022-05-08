clear;
npymatlab_path = '/home/wen/Documents/MEG_ww/npy-matlab-master';
addpath(genpath(npymatlab_path))
BCT_path = '/home/wen/Documents/MEG_ww/2019_03_03_BCT';
addpath(genpath(BCT_path))

%%
leftVnp = '/home/senthilp/caesar/camcan/cc700/freesurfer_output/sub-CC110126/mne_files/sub-CC110126_True_30_12_vcLeft.npy';
rightVnp = '/home/senthilp/caesar/camcan/cc700/freesurfer_output/sub-CC110126/mne_files/sub-CC110126_True_30_12_vcRight.npy';
leftVmat = readNPY(leftVnp);
rightVmat = readNPY(rightVnp);
[r,p] = corr(leftVmat,rightVmat,'type','Spearman')
