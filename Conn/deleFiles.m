myDIr = '/home/senthilp/caesar/camcan/cc700/freesurfer_output';

for i = 1:height(MEGsum)
    subname = MEGsum.participant_id(i,:);
    subFolder = fullfile(myDIr,subname,'mri','shen_corr');
    if isfolder(subFolder)
        subfiles = dir(subFolder);

        del_idx = cellfun(@(x)contains(x,{'_betweeness_','_cluster_','_degree_','_shortPath_'}),{subfiles.name},'UniformOutput',true);
        del_File = cellfun(@(x,y)fullfile(x,y),{subfiles(del_idx).folder},{subfiles(del_idx).name},'UniformOutput',false);
        
        if ~isempty(del_File)
            delete(del_File{:})
            
        end
    end
end