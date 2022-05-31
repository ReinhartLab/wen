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

for sub_i = 1:subN
    subname = MEGsum.participant_id(sub_i,:);
    corrFile = fullfile(dataDir,['sub-' subname],'mri','shen_corr','shen_corr.mat');
    load(corrFile)
corr = squeeze(corr);
    MEGsum.corr{sub_i} = corr(2:end,2:end);%first [0 0 0]
end

%% cpm
MEGsum = MEGsum(1:320,:);

all_mats = cat(3,MEGsum.corr{:});%MxMxN 
all_behav = MEGsum.TOTindx;%Nx1

no_sub = size(all_mats,3);
no_node = size(all_mats,1);

addpath(genpath('/home/wen/Documents/MEG_ww/CPM-master'))
clear corr
y_predict = nan(no_sub,3);
for type_i = 1:3 % network type: positive  negative , integrated
    [y_predict(:,type_i), performance,cpmMask{type_i},cpmFeature{type_i}] = cpm_main(all_mats,all_behav,'pthresh',0.05,'kfolds',no_sub,'type',type_i);
    [cpmR(type_i), cpmPval(type_i)] = corr(y_predict(:,type_i),all_behav,'type','Spearman');
end

save TOTcpm cpmPval cpmR cpmMask cpmFeature y_predict

%%
load TOTcpm
SetFigBasic
cpmStr = {'Positive','Negative' , 'Integrated'};

cmap = brewermap(max(MEGsum.age),'PuOr');
cmap = cmap([MEGsum.age],:);

figure('Position',[100 100 1500 400]);
for type_i = 1:3
    subplot(1,3,type_i);hold all;axis square

  scatter(y_predict(:,type_i), all_behav, 20, cmap, 'filled')

lsline
    text(0.8,0.9,sprintf('subN = %d\nr = %.3f\np = %3.3f',no_sub,cpmR(type_i),cpmPval(type_i)),'sc')
    xlabel('CPM prediction');
    ylabel('Observed beha')
    title(cpmStr{type_i})
end
%% trained network
figure('Position',[100 100 1500 400]);

for type_i = 1:3
    tmp_mask = cpmMask{type_i};
    tmp_mask = reshape(tmp_mask,no_sub,no_node,no_node);
    tmp_mask = squeeze(mean(tmp_mask));%no_node*no_node
    tmp_mask(abs(tmp_mask)~=1) = 0;% select edges survived in every validation
    cpm_mask{type_i} = tmp_mask;

    tmpR = cpmFeature{type_i};
    tmpR = reshape(tmpR,no_sub,no_node,no_node);
    tmpR = squeeze(mean(tmpR));%no_node*no_node

    if type_i ==3
        tmp_mask = abs(tmp_mask);
    end
    cpm_net{type_i} = tmp_mask.*tmpR;

    subplot(1,3,type_i);hold all;axis tight
    imagesc(cpm_net{type_i});colorbar
    title(cpmStr{type_i})
end

%% output  CPM network 268*268
load remain
ShenNodeN = 268;% Shen atlas has 268 parcellations, T1 and MEG registration removed 59ROIs
mat268 = zeros(3,ShenNodeN,ShenNodeN);

idx209 = setdiff(1:ShenNodeN,remain);
for type_i = 1:3
    
    mat268(type_i,idx209,idx209) = logical(cpm_mask{type_i});

    txtFile = ['TOT_',cpmStr{type_i}, '.txt'];%generate binary mask for visualization  https://bioimagesuiteweb.github.io/webapp/connviewer.html
    dlmwrite(txtFile, squeeze(mat268(type_i,:,:)))
end

%% output volume file
type_i = 3;
ShenAtlas = 'shen_2mm_268_parcellation.nii.gz';
V= niftiread();
edgeMat = squeeze(logical(mat268(type_i,:,:)));
edgeV = 0.5*(sum(edgeMat)+sum(edgeMat,2)')./sum(sum(edgeMat==1))*100;
for i = 1:ShenNodeN
    V(V==i) = edgeV(i);
end
niftiwrite(V,'Shen268_2mm_TOT_CPMedgeProportion') 
%%
volFile = fullfile(pwd,'TOT_268_edgeProportion.nii');

surfL = fullfile(pwd,'TOT_268_edgeProportionL');
surfR = fullfile(pwd,'TOT_268_edgeProportionR');
    sprintf('mri_vol2surf --src %s --srcreg register.dat --hemi lh --o %s',volFile,surfL)

    sprintf('mri_vol2surf --src %s --srcreg register.dat --hemi rh --o %s',volFile,surfR)
%%
networkNames = {'MFN', 'FPN', 'DMN', 'SubC', 'MON','VisI', 'VisII', 'VA'};
%medial frontal network;frontoparietal network;default mode network ;subcortical and cerebellar regions ;motor network; visualI network ;visual II network;visualassociation network

shen_net_table = tdfread('shen_268_parcellation_networklabels.csv');
shen_net_table = array2table(shen_net_table.Node0x2CNetwork,'VariableNames',{'node','network'});
shen_net_table.label = networkNames([shen_net_table.network])';

%%
%% 
type_i = 3;
edgeMat = squeeze(logical(mat268(type_i,:,:)));
netMat = zeros(8);
netN = zeros(8);
for r = 1:ShenNodeN
    nr = shen_net_table.network(r);
    for c = 1:ShenNodeN
        nc = shen_net_table.network(c);
        netN(nr,nc) = netN(nr,nc) +1;
        if edgeMat(r,c) == 1
            netMat(nr,nc) = netMat(nr,nc)+1;

        end

    end
end
%%
figure;
heatmap(networkNames,networkNames,netMat,'Colormap',copper)

%%
% thresh = 0.01;%threshold for feature selection
% cpm.thresh = thresh;
% 
% 
% 
% behav_pred_pos = zeros(no_sub,1);
% behav_pred_neg = zeros(no_sub,1);
% behav_pred = zeros(no_sub,1);
% 
% for leftout = 1:no_sub
%     fprintf('\n Leaving out subj # %4.0f',leftout);
%     
%     % leave out subject from matrices and behavior
%     
%     train_mats = all_mats;
%     train_mats(:,:,leftout) = [];
%     train_vcts = reshape(train_mats,[],size(train_mats,3));
%     
%     train_behav = all_behav;
%     train_behav(leftout) = [];
%     
%     % correlate all edges with behavior using robust regression
%     edge_no = size(train_vcts,1);
%     r_mat = zeros(1, edge_no);
%     p_mat = zeros(1, edge_no);
%     
%     for edge_i = 1: edge_no
%         [~, stats] = robustfit(train_vcts(edge_i,:)', train_behav);
%         cur_t = stats.t(2);
%         r_mat(edge_i) = sign(cur_t)*sqrt(cur_t^2/(no_sub-1-2+cur_t^2));
%         p_mat(edge_i) = 2*(1-tcdf(abs(cur_t), no_sub-1-2));  %two tailed
%     end
% 
%     %     % correlate all edges with behavior using partial correlation
% %     [r_mat, p_mat] = partialcorr(train_vcts', train_behav, age);
% %     
% %        
% %     % correlate all edges with behavior using rank correlation
% %     [r_mat, p_mat] = corr(train_vcts', train_behav, 'type', 'Spearman');
%     
%         
%     r_mat = reshape(r_mat,no_node,no_node);
%     p_mat = reshape(p_mat,no_node,no_node);
%     
%     % set threshold and define masks 
%     pos_mask = zeros(no_node, no_node);
%     neg_mask = zeros(no_node, no_node);
%     
%     
%     pos_edge = find( r_mat >0 & p_mat < thresh);
%     neg_edge = find( r_mat <0 & p_mat < thresh);
%     
%     pos_mask(pos_edge) = 1;
%     neg_mask(neg_edge) = 1;
% 
%     %     %-----------------sigmoidal weighting---------------------------%
% %     pos_edges = find(r_mat > 0 );
% %     neg_edges = find(r_mat < 0 );
% %     
% %     % covert p threshold to r threshold
% %     T = tinv(thresh/2, no_sub-1-2);
% %     R = sqrt(T^2/(no_sub-1-2+T^2));
% %     
% %     % create a weighted mask using sigmoidal function
% %     % weight = 0.5, when correlation = R/3;
% %     % weight = 0.88, when correlation = R;
% %     pos_mask(pos_edges) = sigmf( r_mat(pos_edges), [3/R, R/3]);
% %     neg_mask(neg_edges) = sigmf( r_mat(neg_edges), [-3/R, R/3]);
% %     %---------------sigmoidal weighting-----------------------------%
%     
%     % get sum of all edges in TRAIN subs (divide by 2 to control for the
%     % fact that matrices are symmetric)
%     
%     train_sumpos = zeros(no_sub-1,1);
%     train_sumneg = zeros(no_sub-1,1);
%     
%     for ss = 1:size(train_sumpos)
%         train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*pos_mask))/2;
%         train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*neg_mask))/2;
%     end
%     
%     % build model on TRAIN subs
%     % combining both postive and negative features
%     b = regress(train_behav, [train_sumpos, train_sumneg, ones(no_sub-1,1)]);
%     
%     % run model on TEST sub
%     
%     test_mat = all_mats(:,:,leftout);
%     test_sumpos = sum(sum(test_mat.*pos_mask))/2;
%     test_sumneg = sum(sum(test_mat.*neg_mask))/2;
%     
%     behav_pred(leftout) = b(1)*test_sumpos + b(2)*test_sumneg + b(3);
%     
%     
% end
% clear corr
% % compare predicted and observed scores
% [cpmR, cpmPval] = corr(behav_pred,all_behav)
% 
% % [R_pos, P_pos] = corr(behav_pred_pos,all_behav)
% % [R_neg, P_neg] = corr(behav_pred_neg,all_behav)
% figure; plot(behav_pred,all_behav,'r.'); lsline
% % figure(1); plot(behav_pred_pos,all_behav,'r.'); lsline
% % figure(2); plot(behav_pred_neg,all_behav,'b.'); lsline