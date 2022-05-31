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
[nofs,uniqIA]= setdiff(MEGsum.participant_id,fs_folder);
MEGsum(uniqIA,:)=[];
fprintf('Setdiff T1 removed %d ppts\n',length(nofs));

%%
subN = height(MEGsum);
% freqList = [4,8,12,16,20];
freqList = [4,8];

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

tmp = cellfun(@length,MEGsum.corr(:,2),'UniformOutput',true);
MEGsum(tmp<268,:) = [];


%%
networkNames = {'MFN','FPN','DMN','SC', 'MON','Vis1', 'Vis2', 'VA'};
networkN = length(networkNames);
%medial frontal network;frontoparietal network;default mode network ;
% subcortical and cerebellar regions ;motor network; visualI  ;visual II ;visual association

shen_net_table = tdfread('shen_268_parcellation_networklabels.csv');
shen_net_table = array2table(shen_net_table.Node0x2CNetwork0x2Chemisphere,'VariableNames',{'node','network','hemi'});
shen_net_table.label = networkNames([shen_net_table.network])';

[shen_net_table, sort_index] = sortrows(shen_net_table,{'label'});
[uniqName,uniqIA] = unique(shen_net_table.label,'stable');
%% plot age correlation
r_thresh = 0.01;
for fi = 1:FN
    tmp_Data = cat(3,MEGsum.corr{:,fi});
    tmp = size(tmp_Data);
    tmp_array = reshape(tmp_Data,tmp(1)*tmp(2),tmp(3));

    p = nan(length(tmp_array),1);
    rho = p;
    parfor i = 1:length(tmp_array)
        [rho(i),p(i)] = corr(tmp_array(i,:)',MEGsum.age,'type','Spearman');
    end
    AgeCorr.p{fi} = reshape(p,tmp(1),tmp(2));
    AgeCorr.h{fi} = AgeCorr.p{fi}<r_thresh;

    AgeCorr.rho{fi} = reshape(rho,tmp(1),tmp(2));
end
%%
MyMap = brewermap(100, 'RdBu'); %with 100 being the number of steps in your colormap
myMap = flipud(MyMap); %Since we want blue for negative and red for positive values
SetFigBasic

figure('name','Corr age and edge strength')
for  fi = 1:FN
    subplot(2,ceil(FN/2),fi);
    axis square;hold all

    dat = AgeCorr.rho{fi}.*AgeCorr.h{fi};
    dat = dat(sort_index,sort_index);

    uniqIA2 = [uniqIA;height(shen_net_table)+1];
    for i = 1:networkN
        for j = 1:networkN
            netCorr(fi,i,j) = mean(mean(dat(uniqIA2(i):uniqIA2(i+1)-1,uniqIA2(j):uniqIA2(j+1)-1),'omitnan'),'omitnan');
        end
    end

    imagesc(squeeze(netCorr(fi,:,:)));
    clim([-0.1 0.1])

    colormap(myMap)
    cb = colorbar;
    cb.Title.String = 'Rho';
    cb.Ticks = -0.5:0.1:0.5;
    box on
    set(gca,'ydir','normal','xlim',[1 networkN],'ylim',[1 networkN],'xtick',1:networkN,'xticklabel',uniqName,'ytick',1:networkN,'yticklabel',uniqName)
    xlabel('Network');ylabel('Network')
    title(sprintf('(N=%d), %s',height(MEGsum),freqStr{fi}))

end
%%
% https://stackoverflow.com/questions/46395193/create-stacked-2d-matrix-along-z-axis-with-surf-in-matlab

figure('Position',[200 200 800 1200]);
axis square;hold all
offset = 0.5;
for  fi = 1:FN
    dat = AgeCorr.rho{fi}.*AgeCorr.h{fi};
    dat = dat(sort_index,sort_index);

    h = surf(repmat(offset*(fi-1), size(dat)), dat);
    set(h,'edgecolor','none')
end
campos([-100 -170 90])
clim([-0.3 0.3])
view(45,-45)
grid off
rotate3d on
xlabel('Node');ylabel('Node');zlabel('Frequency (Hz)')
set(gca,'ydir','reverse','xlim',[0 height(shen_net_table)+1],'ylim',[0 height(shen_net_table)+5],'zlim',[0 FN*offset], ...
    'xtick',uniqIA+5,'xticklabel',uniqName,'ytick',uniqIA+2,'yticklabel',uniqName,'ztick',0:1:FN,'ZTickLabel',freqList(1:2:end))
colormap(myMap)
cb = colorbar;
cb.Title.String = 'Rho';
cb.Ticks = -0.5:0.1:0.5;

%%

