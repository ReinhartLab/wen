clear;
addpath(genpath('/home/wen/Documents/MATLAB/Toolbox/npy-matlab-master'))

load subs_info.mat

MEGsum = tdfread(AllInfo.sumFile{1});
MEGsum = struct2table(MEGsum);
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
[subsNOmegfs,uniqIA]= setdiff(MEGsum.participant_id,fs_folder);
MEGsum(uniqIA,:)=[];
fprintf('Setdiff T1 removed %d ppts\n',length(subsNOmegfs));

%%
% MEGsum = MEGsum(1:32,:);    
subN = height(MEGsum);

SUBDIR = '/home/senthilp/caesar/camcan/cc700/freesurfer_output';
% freqList = [4,8,12,20];

G.degree = nan(subN,length(freqList),268);
G.between = G.degree;
G.path = G.between;
G.cluster = G.between;

for sub_i = 1:subN
    
    subname  = MEGsum.participant_id(sub_i,:);
    for fi = 1:length(freqList)
        freq = freqList(fi);

        npyName = fullfile(SUBDIR,subname,'mri','shen_corr',[subname '_degree_' num2str(freq) 'Hz.npy']);
      
        if isfile(npyName)
        G.degree(sub_i,fi,:) = readNPY(npyName);

        npyName = fullfile(SUBDIR,subname,'mri','shen_corr',[subname '_betweeness_' num2str(freq) 'Hz.npy']);
        G.between(sub_i,fi,:) = readNPY(npyName);
 %(betweenness centrality);the relative frequency of how often that node is part of shortest paths

        npyName = fullfile(SUBDIR,subname,'mri','shen_corr',[subname '_shortPath_' num2str(freq) 'Hz.npy']);
        G.path(sub_i,fi,:) = readNPY(npyName);

        npyName = fullfile(SUBDIR,subname,'mri','shen_corr',[subname '_cluster_' num2str(freq) 'Hz.npy']);
        G.cluster(sub_i,fi,:) = readNPY(npyName);
      end

    end

end
rmvID = isnan(G.degree(:,1,1));
MEGsum(rmvID,:) = [];
G.degree(rmvID,:,:) = [];
G.between(rmvID,:,:) = [];
G.path(rmvID,:,:) = [];
G.cluster(rmvID,:,:) = [];
%%
ageGroup.Bin = 17:2:29;
[ageGroup.N,ageGroup.edges] = histcounts(MEGsum.age,ageGroup.Bin);
ageGroup.mid = (ageGroup.edges(1:end-1)+ageGroup.edges(2:end))./2;
plotMat = [];

for ageID = 1:length(ageGroup.N)
    tmpID = discretize(MEGsum.age,ageGroup.Bin)==ageID;
    for fi = 1:size(G.degree,2)
        freq = freqList(fi);
        plotMat(1,ageID,fi,:) = mean(G.degree(tmpID,fi,:));
        plotMat(2,ageID,fi,:) = mean(G.between(tmpID,fi,:));
        plotMat(3,ageID,fi,:) = mean(G.path(tmpID,fi,:));
        plotMat(4,ageID,fi,:) = mean(G.cluster(tmpID,fi,:));

    end
end
%%
addpath(genpath('/home/wen/Documents/MATLAB/Toolbox/Cbrewer'))
MyMap = brewermap(100, 'RdBu'); %with 100 being the number of steps in your colormap
myMap = flipud(MyMap); %Since you want blue for negative and red for positive values

figure;

titleStr = {'Degree','Betweenness','Path','cluster'};
for i = 1:size(plotMat,1)
subplot(1,size(plotMat,1),i);
hold all;axis square
dat = squeeze(mean(plotMat(i,:,:,:),4))';
imagesc(dat);


%  clim([-0.1 0.4])

    colormap(myMap)
    cb = colorbar;
%     cb.Title.String = 'Rho';
    %     caxis([-0.5 0.5])

xlabel('Age');ylabel('Frequency(Hz)');
set(gca,'YDir','normal','ytick',1:size(plotMat,3),'YTickLabel',freqList,'xtick',1:size(plotMat,2),'XTickLabel',ageGroup.mid)
title(titleStr{i})
end
%% age correlation
threshR = .01;

dat = reshape(G.degree,size(G.degree,1),size(G.degree,2)*size(G.degree,3));
[r,p] = corr(MEGsum.age,dat);
G.r.degree = reshape(r,size(G.degree,2),size(G.degree,3));
G.p.degree = reshape(p,size(G.degree,2),size(G.degree,3));
degree = G.r.degree.*(G.p.degree<threshR);

dat = reshape(G.between,size(G.between,1),size(G.between,2)*size(G.between,3));
[r,p] = corr(MEGsum.age,dat);
G.r.between = reshape(r,size(G.between,2),size(G.between,3));
G.p.between = reshape(p,size(G.between,2),size(G.between,3));
between = G.r.between.*(G.p.between<threshR);

dat = reshape(G.path,size(G.path,1),size(G.path,2)*size(G.path,3));
[r,p] = corr(MEGsum.age,dat);
G.r.path = reshape(r,size(G.path,2),size(G.path,3));
G.p.path = reshape(p,size(G.path,2),size(G.path,3));
path = G.r.path.*(G.p.path<threshR);

dat = reshape(G.cluster,size(G.cluster,1),size(G.cluster,2)*size(G.cluster,3));
[r,p] = corr(MEGsum.age,dat);
G.r.cluster = reshape(r,size(G.cluster,2),size(G.cluster,3));
G.p.cluster = reshape(p,size(G.cluster,2),size(G.cluster,3));
cluster = G.r.cluster.*(G.p.cluster<threshR);

save Graph_cluster.mat cluster
save Graph_degree.mat degree
save Graph_between.mat between
save Graph_path.mat path