clear;
rng('shuffle'); 
addpath(genpath('/home/wen/Documents/MATLAB/Toolbox/Cbrewer'))

%%

load subs_info.mat

MEGsum = tdfread(AllInfo.sumFile{1});
MEGsum = struct2table(MEGsum);
MEGsum.participant_id(:,1:4) = [];
if length(unique(cellstr(MEGsum.participant_id))) == height(MEGsum)
    fprintf('--------\nrest meg participants N:%d\n',height(MEGsum))
else
    error('repeated entries')
end

%% get the shared ppts in freesurfer_output
dataDir = '/home/senthilp/caesar/camcan/cc700/freesurfer_output';
fs_folder = dir(dataDir);
subFolder = contains({fs_folder.name},'sub-CC');
fs_folder = {fs_folder([fs_folder.isdir]&subFolder).name}';
fs_folder = cellfun(@(x) x(5:end),fs_folder,'UniformOutput',false);
[subsNOmegfs,uniqIA]= setdiff(MEGsum.participant_id,fs_folder);
MEGsum(uniqIA,:)=[];
fprintf('Setdiff T1 removed %d ppts\n',length(subsNOmegfs));

%%
subN = height(MEGsum);
% freqList = [4,8,12,20];
freqList = [4,4];
freqStr = num2cell(freqList);

freqStr = cellfun(@(x)[num2str(x),'Hz'],freqStr,'UniformOutput',false);
FN = length(freqStr);

for sub_i = 1:subN
    subname = MEGsum.participant_id(sub_i,:);

    for fi = 1:FN
        corrFile = fullfile(dataDir,['sub-' subname],'mri','shen_corr',['shen_corr_',freqStr{fi},'_lcmv.mat']);
        if isfile(corrFile)
            load(corrFile)
            MEGsum.corr{sub_i,fi} = atanh(corMat); % fisherZtransform
        end
    end
end

tmp = cellfun(@isempty,MEGsum.corr(:,1),'UniformOutput',true);
MEGsum(tmp,:) = [];

tmp = cellfun(@length,MEGsum.corr(:,1),'UniformOutput',true);
MEGsum(tmp<268,:) = [];
%%

networkNames = {'MFN','FPN','DMN','SC', 'MON','Vis1', 'Vis2', 'VA'};

%medial frontal network;frontoparietal network;default mode network ;
% subcortical and cerebellar regions ;motor network; visualI  ;visual II ;visual association

shen_net_table = tdfread('shen_268_parcellation_networklabels.csv');
shen_net_table = array2table(shen_net_table.Node0x2CNetwork0x2Chemisphere,'VariableNames',{'node','network','hemi'});
shen_net_table.label = networkNames([shen_net_table.network])';

nodeN = height(shen_net_table);

shen_net_table_hemi = shen_net_table;
tmp_idx = shen_net_table_hemi.hemi==2;
shen_net_table_hemi.label(tmp_idx) = cellfun(@(x)[x '_R'],shen_net_table_hemi.label(tmp_idx),'UniformOutput',0);
tmp_idx = shen_net_table_hemi.hemi==1;
shen_net_table_hemi.label(tmp_idx) = cellfun(@(x)[x '_L'],shen_net_table_hemi.label(tmp_idx),'UniformOutput',0);
[shen_net_table_hemi, sort_index] = sortrows(shen_net_table_hemi,{'hemi','label'});

[uniqName,uniqIA] = unique(shen_net_table_hemi.label,'stable');

tmp = diff([uniqIA;nodeN+1]);
a = [];
for i = 1:length(uniqName)
    a = [a 1:tmp(i)];
end
shen_net_table_hemi.withinID = a';

%% generate the network_general.csv
fileDir = '/home/senthilp/caesar/wen/CircosPlot/circos_data';
nodeLength = 100;
gntable = table;
gntable.network = unique(shen_net_table_hemi.label,'stable');
gntable.nodeN = diff([uniqIA;nodeN+1]);

gntable.start = ones(height(gntable),1);
gntable.end = nodeLength*gntable.nodeN;

gntable(1+height(gntable)/2:height(gntable),:) = flip(gntable(1+height(gntable)/2:height(gntable),:));% make it symmetrical in circle

gntable = flip(gntable);% R on the right side

csvname = fullfile(fileDir,'NodeNetwork_general.csv');
writetable(gntable,csvname,'WriteRowNames',true)
%% generate the dot_plot

dots_table = table;
for i = 1:height(gntable)
    tmp_Table = table;

    tmp_Table.network(1:gntable.nodeN(i)) = gntable.network(i);
    tmp_Table.center(1:gntable.nodeN(i)) = nodeLength*(1:gntable.nodeN(i))-nodeLength/2;

    dots_table = [dots_table;tmp_Table];
end
dots_table.value = ones(height(dots_table),1);

csvname = fullfile(fileDir,'NodeNetwork_dots.csv');
writetable(dots_table,csvname,'WriteRowNames',true)

%% threshold based on hipp 2012
ageGroup.Bin = [18 30 40 55 88];
% ageGroup.Bin = [18 35 50 65 88];

[ageGroup.N,ageGroup.edges] = histcounts(MEGsum.age,ageGroup.Bin);


for ageID = [1 length(ageGroup.N)]
    tmpID = discretize(MEGsum.age,ageGroup.Bin)==ageID;
    for fi = 1:FN
        tmp_Data = cat(3,MEGsum.corr{tmpID,fi});%268*268*subs
        tmp_Data = tmp_Data(sort_index,sort_index,:); % sort based one table

        p = nan(size(tmp_Data,1),size(tmp_Data,2),size(tmp_Data,3));
        parfor subj = 1:size(tmp_Data,3)
            p(:,:,subj) = hipp_threshold_ttest(squeeze(tmp_Data(:,:,subj)),0.001);
            fprintf('ageID =%d,fi = %d,subj =%d, done.\r',ageID,fi,subj)
        end
        hipp.h{ageID,fi} = squeeze(mean(p,3,"omitnan")>.8);% XX% ppts have this connection
        hipp.evR{ageID,fi} = hipp.h{ageID,fi}.*squeeze(mean(tmp_Data,3)); %envelop correlation value
        save(sprintf('hipp_tmp_A%d_F%d.mat',ageID,fi),'hipp','-v7.3')

        % generate line excel

        lineTable = shen_net_table_hemi(:,["label","withinID"]);

        lineTable.end = lineTable.withinID*nodeLength;
        lineTable.start = lineTable.end-nodeLength;
        lineTable.value1 = squeeze(mean(mean(tmp_Data,3,'omitnan'),'omitnan'))';
        lineTable.value2 = squeeze(mean(mean(tmp_Data,3,'omitnan'),2,'omitnan'));

        csvname = fullfile(fileDir,sprintf('NodeConn_line_%s_%dto%d_hipp.csv',freqStr{fi},ageGroup.Bin(ageID),ageGroup.Bin(ageID+1)));
        writetable(lineTable,csvname,'WriteRowNames',true)

        % generate circos link excel
        link_width= 40;

        dat = hipp.h{ageID,fi};

        a = tril(dat,-1);
        b = triu(dat,1);
        dat = ((a+b')./2)==1;

        x = repmat(shen_net_table_hemi,nodeN,1);
        x(:,{'network','hemi'}) = [];
        x.Properties.VariableNames{'label'} = 'network';
        y = sortrows(x,'network');
        y = renamevars(y,["node","network","withinID"], ["node2","network2","withinID2"]);
        linktable  = [x y] ;

        linktable.value = reshape(dat,size(dat,1)*size(dat,2),1);
        linktable(isnan(linktable.value) | linktable.value==0 | linktable.node == linktable.node2,:) = [];% remove insignificant connections

        linktable.start1 = linktable.withinID*nodeLength - nodeLength/2-link_width/2;
        linktable.end1 = linktable.start1+link_width;

        linktable.start2 = linktable.withinID2*nodeLength - nodeLength/2-link_width/2;
        linktable.end2 = linktable.start2+link_width;

        linktable = linktable(:,["network","withinID","start1","end1","network2","withinID2","start2","end2"]);

        csvname = fullfile(fileDir,sprintf('NodeConn_linkplot_%s_%dto%d_hipp.csv',freqStr{fi},ageGroup.Bin(ageID),ageGroup.Bin(ageID+1)));
        writetable(linktable,csvname,'WriteRowNames',true)

    end
end
save('hipp_thresh.mat','hipp','-v7.3')
