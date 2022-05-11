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

%% get the shared ppts in freesurfer_output
dataDir = '/home/senthilp/caesar/camcan/cc700/freesurfer_output';
fs_folder = dir(dataDir);
subFolder = contains({fs_folder.name},'sub-CC');
fs_folder = {fs_folder([fs_folder.isdir]&subFolder).name}';
fs_folder = cellfun(@(x) x(5:end),fs_folder,'UniformOutput',false);
[subsNOmegfs,uniqIA]= setdiff(MEGsum.participant_id,fs_folder);
MEGsum(uniqIA,:)=[];
fprintf('Setdiff removed %d ppts\n',length(subsNOmegfs));

%% remove
% CC410129 shen Right-Claustrum 139 out side of ppt's T1_brain
% CC420464 ['Right-Claustrum', 'Vitreous-Humor']
rmvSubs = {'CC410129','CC420464'};
rmvID = find(ismember(cellstr(MEGsum.participant_id),rmvSubs));
MEGsum(rmvID,:) = [];
fprintf('Registration removed %d ppts\n',length(rmvID));

%%
subN = height(MEGsum);

freqStr = {'5Hz','5Hz'};%,'Beta','Gamma','HighGamma'};
FN = length(freqStr);
for sub_i = 1:subN
    subname = MEGsum.participant_id(sub_i,:);

    for fi = 1:FN
        corrFile = fullfile(dataDir,['sub-' subname],'mri','shen_corr',['shen_corr_',freqStr{fi},'_lcmv.mat']);
        if isfile(corrFile)
            load(corrFile)
            MEGsum.corr{sub_i,fi} = squeeze(corMat(2:end,2:end));%first [0 0 0]
            MEGsum.corr{sub_i,fi} = atanh(MEGsum.corr{sub_i,fi}); % fisherZtransform
        end
    end
end

tmp = cellfun(@isempty,MEGsum.corr(:,1),'UniformOutput',false);
MEGsum(cell2mat(tmp),:) = [];

%%
% ageGroup.Bin = [18 35 50 65 90];
ageGroup.Bin = [18:3:30];
[ageGroup.N,ageGroup.edges] = histcounts(MEGsum.age,ageGroup.Bin)
for ageID = 1:length(ageGroup.N)
    tmpID = discretize(MEGsum.age,ageGroup.Bin)==ageID;
    for fi = 1:FN
        tmp_Data = cat(3,MEGsum.corr{tmpID,fi});
        tmp = size(tmp_Data);
        tmp_array = reshape(tmp_Data,tmp(1)*tmp(2),tmp(3));

        p = nan(length(tmp_array),1);
        t = p;
        parfor i = 1:length(tmp_array)
            [~,p(i),~,stat] = ttest(tmp_array(i,:),0,'Alpha',.01,'tail','both');
            t(i) = stat.tstat;
        end

        p = reshape(p,tmp(1),tmp(2));
        fdR.t{ageID,fi} = reshape(t,tmp(1),tmp(2));

        [fdR.h{ageID,fi}, fdR.crit_p{ageID,fi}, fdR.adj_ci_cvrg{ageID,fi}, fdR.adj_p{ageID,fi}]=fdr_bh(p,0.01,'pdep','yes');
    end
end

%%
load removedNodes.mat
networkNames = { 'DMN','FPN', 'MFN', 'MON','SC', 'VA','VisI', 'VisII'};
%medial frontal network;frontoparietal network;default mode network ;subcortical and cerebellar regions ;motor network; visualI network ;visual II network;visualassociation network

shen_net_table = tdfread('shen_268_parcellation_networklabels.csv');
shen_net_table = array2table(shen_net_table.Node0x2CNetwork,'VariableNames',{'node','network'});
shen_net_table.label = networkNames([shen_net_table.network])';

shen_net_table(rmvdNodes,:) = [];

[shen_net_table, sort_index] = sortrows(shen_net_table,"label");
nodeN = height(shen_net_table);

[uniqName,uniqIA]=unique(shen_net_table.label);

%% generate the network_general.csv
fileDir = '/home/senthilp/caesar/wen/Conn/CircosPlot/circos_data';
nodeLength = 100;
gntable = table;
gntable.network = uniqName;
gntable.start = ones(length(uniqName),1);
gntable.end = [diff([uniqIA;nodeN])]*nodeLength;

csvname = fullfile(fileDir,'NodeNetwork_general.csv');
writetable(gntable,csvname,'WriteRowNames',true)

%%  node link plot file

for ageID = 1:length(ageGroup.N)
    for  fi = 1:FN
        for dirFlag = 1:2
            dat = fdR.h{ageID,fi};
            dat = dat(sort_index,sort_index);

            a = tril(dat,-1);
            b = triu(dat,1);

            if dirFlag ==1         % unidirectional connectivity
                dat = logical(a+b');
                csvname = fullfile(fileDir,sprintf('NodeConnUni_linkplot_%s_%dto%d.csv',freqStr{fi},ageGroup.Bin(ageID),ageGroup.Bin(ageID+1)));
            else
                dat = (a+b')==2; %bidirectional connectivity
                csvname = fullfile(fileDir,sprintf('NodeConnBi_linkplot_%s_%dto%d.csv',freqStr{fi},ageGroup.Bin(ageID),ageGroup.Bin(ageID+1)));
            end

            dat = tril(dat);

            x = repmat(shen_net_table,nodeN,1);
            x(:,"network") = [];
            x.Properties.VariableNames{'label'} = 'network';
            y = sortrows(x,'node');
            y = renamevars(y,["node","network"], ["node2","network2"]);

            linktable  = [x y] ;
            linktable.value = reshape(dat,size(dat,1)*size(dat,2),1);
            linktable(isnan(linktable.value) | linktable.value==0,:) = [];% remove insignificant connections

            linktable.start1 = linktable.node*100;
            linktable.end1 = linktable.start1+nodeLength;

            linktable.start2 = linktable.node2*100;
            linktable.end2 = linktable.start2+nodeLength;

            linktable = linktable(:,[2 6 7 4 8 9 1 3 5]);

            writetable(linktable,csvname,'WriteRowNames',true)

        end
    end
end
% %% plotting network correlation
% 
% figure;
% for ageID = 1:length(ageGroup.N)
%     for  fi = 1:FN
%         subplot(4,2,(ageID-1)*2+fi);
%         axis square;hold all
% 
%         dat  = squeeze(netCon(ageID,fi,:,:));
%         imagesc(dat);
%         cb = colorbar;
%         cb.Title.String = 't-value';
%         caxis([8 16])
% 
%         set(gca,'ydir','normal','xlim',[1 length(uniqName)],'ylim',[1 length(uniqName)],'xticklabel',uniqName,'yticklabel',uniqName)
%         xlabel('Network');ylabel('Network')
%         xtickangle(45)
%         title(sprintf('Age %d~%d(N=%d),%s',ageGroup.edges(ageID),ageGroup.edges(ageID+1),ageGroup.N(ageID),freqStr{fi}))
%     end
% end
% 
% 
% %% plot age correlation
% 
% for fi = 1:FN
%     tmp_Data = cat(3,MEGsum.corr{:,fi});
%     tmp = size(tmp_Data);
%     tmp_array = reshape(tmp_Data,tmp(1)*tmp(2),tmp(3));
% 
%     p = nan(length(tmp_array),1);
%     rho = p;
%     parfor i = 1:length(tmp_array)
%         [rho(i),p(i)] = corr(tmp_array(i,:)',MEGsum.age,'type','Spearman');
%     end
%     AgeCorr.p{fi} = reshape(p,tmp(1),tmp(2));
%     AgeCorr.h{fi} = AgeCorr.p{fi}<.05;
% 
%     AgeCorr.rho{fi} = reshape(rho,tmp(1),tmp(2));
% %     AgeCorr.hPN{fi} = (AgeCorr.p{fi}<.05 & AgeCorr.rho{fi}>0) - (AgeCorr.p{fi}<.05 & AgeCorr.rho{fi}<0);
% 
% end
% 
% 


