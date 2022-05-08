clear;
rng('shuffle');

load subs_info.mat

MEGsum = tdfread(AllInfo.sumFile{1});
MEGsum = struct2table(MEGsum);
MEGsum.participant_id(:,1:4) = [];
if length(unique(cellstr(MEGsum.participant_id))) == height(MEGsum)
    fprintf('--------\nrest meg participants N:%d\n',height(MEGsum))
else
    error('repeated entries')
end

task_i = 4;%tot
curTask = AllInfo.task{task_i};
filename = AllInfo.sumFile{task_i};
tasksum = tdfread(filename,'\t');
tasksum = struct2table(tasksum);
tasksum = sortrows(tasksum,'Subject');

if length(unique(cellstr(tasksum.Subject))) == height(tasksum)
    fprintf('%s participants N:%d\n',curTask,height(tasksum))
else
    error('repeated entries')
end

MEGsum.TOTindx = nan(height(MEGsum),1);
for sub_i = 1:height(MEGsum)
    [lia,locb] = ismember(MEGsum.participant_id(sub_i,:),cellstr(tasksum.Subject));
    if lia
        MEGsum.TOTindx(sub_i) = tasksum.ToT_ratio(locb);
    end
end

MEGsum(isnan(MEGsum.TOTindx),:)  = [];
subN = height(MEGsum);
%% get the shared ppts in freesurfer_output
dataDir = '/home/senthilp/caesar/camcan/cc700/freesurfer_output';
fs_folder = dir(dataDir);
subFolder = contains({fs_folder.name},'sub-CC');
fs_folder = {fs_folder([fs_folder.isdir]&subFolder).name}';
fs_folder = cellfun(@(x) x(5:end),fs_folder,'UniformOutput',false);
[subsNOmegfs,ia]= setdiff(MEGsum.participant_id,fs_folder);
MEGsum(ia,:)=[];
fprintf('Setdiff removed %d ppts\n',length(subsNOmegfs));

%% remove
% CC410129 shen Right-Claustrum 139 out side of ppt's T1_brain
% CC420464 ['Right-Claustrum', 'Vitreous-Humor']
rmvSubs = {'CC410129','CC420464'};
rmvID = find(ismember(cellstr(MEGsum.participant_id),rmvSubs));
MEGsum(rmvID,:) = [];
fprintf('Registration removed %d ppts\n',length(rmvID));

%%

freqStr = {'Theta','Alpha','Beta'};%,'Gamma','HighGamma'};
FN = length(freqStr);
for sub_i = 1:subN
    subname = MEGsum.participant_id(sub_i,:);
    for fi = 1:FN
    corrFile = fullfile(dataDir,['sub-' subname],'mri','shen_corr',['shen_corr_',freqStr{fi},'.mat']);
    load(corrFile)

    MEGsum.corr{sub_i,fi} = squeeze(corMat(1,2:end,2:end));%first [0 0 0]
    end
end

%% cpm
MEGsum = MEGsum(1:31,:);

addpath(genpath('/home/wen/Documents/MEG_ww/CPM-master'))
all_behav = MEGsum.TOTindx;%Nx1
no_sub = length(all_behav);
y_predict = nan(no_sub,3,FN);

for fi = 1:FN
    all_mats = cat(3,MEGsum.corr{:,fi});%MxMxN

    no_node = size(all_mats,1);

    for type_i = 1:3 % network type: positive  negative , integrated
        [y_predict(:,type_i,fi), performance,cpmMask{type_i,fi},cpmFeature{type_i,fi}] = cpm_main(all_mats,all_behav,'pthresh',0.05,'kfolds',no_sub,'type',type_i);
        [cpmR(type_i,fi), cpmPval(type_i,fi)] = corr(y_predict(:,type_i,fi),all_behav,'type','Spearman');
    end
end
%%
networkNames = {'MFN', 'FPN', 'DMN', 'SC', 'MON','VisI', 'VisII', 'VA'};
%medial frontal network;frontoparietal network;default mode network ;subcortical and cerebellar regions ;motor network; visualI network ;visual II network;visualassociation network

shen_net_table = tdfread('shen_268_parcellation_networklabels.csv');
shen_net_table = array2table(shen_net_table.Node0x2CNetwork,'VariableNames',{'node','network'});
shen_net_table.label = networkNames([shen_net_table.network])';
%%
SetFigBasic
cmap = brewermap(max(MEGsum.age),'PuOr');
cmap = cmap([MEGsum.age],:);

cpmStr = {'Positive','Negative' , 'Integrated'};
figure;
for type_i = 1:3
    for fi = 1:FN
    subplot(3,FN,(type_i-1)*FN+fi);hold all;axis square

      scatter(y_predict(:,type_i,fi), all_behav, 20, cmap, 'filled')

lsline
    text(0.8,0.9,sprintf('subN = %d\nr = %.3f\np = %3.3f',no_sub,cpmR(type_i),cpmPval(type_i)),'sc')
    xlabel('CPM prediction');
    ylabel('Observed beha')
    title([cpmStr{type_i} '-' freqStr{fi}])
    end
end
%% trained network
figure('Position',[100 100 1500 400]);

for type_i = 1:3
    for fi = 1:FN
    tmp_mask = cpmMask{type_i,fi};
    tmp_mask = reshape(tmp_mask,no_sub,no_node,no_node);
    tmp_mask = squeeze(mean(tmp_mask));%no_node*no_node
    tmp_mask(abs(tmp_mask)~=1) = 0;% select edges survived in every validation
    cpm_mask{type_i,fi} = tmp_mask;

    tmpR = cpmFeature{type_i,fi};
    tmpR = reshape(tmpR,no_sub,no_node,no_node);
    tmpR = squeeze(mean(tmpR));%no_node*no_node

    if type_i ==3
        tmp_mask = abs(tmp_mask);
    end
    cpm_net{type_i,fi} = tmp_mask.*tmpR;

    subplot(3,FN,(type_i-1)*FN+fi);hold all;axis tight
    imagesc(cpm_net{type_i,fi});colorbar
    title([cpmStr{type_i} '-' freqStr{fi}])
    end
end
%% output  CPM network 268*268
load remain

idx209 = setdiff(1:268,remain);
for type_i = 1:3
    mat268 = zeros(268);

    mat268(idx209,idx209) = cpm_net{type_i};

    txtFile = ['TOT_',cpmStr{type_i}, '.txt'];
    dlmwrite(txtFile, mat268)
end

%%
