clear;
load subs_info.mat
taskROI = AllInfo.task;
taskROI(contains(taskROI,'MEG_rest')) = [];

%% MEG_Rest

MEGsum = tdfread(AllInfo.sumFile{1});
MEGsum = struct2table(MEGsum);
MEGsum.participant_id(:,1:4) = [];
if length(unique(cellstr(MEGsum.participant_id))) == height(MEGsum)
    fprintf('--------\nrest meg participants N:%d\n',height(MEGsum))
else
    error('repeated entries')
end
%% fluid IQ

task_i = 2;
curTask = AllInfo.task{task_i};
filename = AllInfo.sumFile{task_i};
tasksum = tdfread(filename,'\t');
tasksum = struct2table(tasksum);
tasksum = sortrows(tasksum,'CCID');

if length(unique(cellstr(tasksum.CCID))) == height(tasksum)
    fprintf('%s participants N:%d\n',curTask,height(tasksum))
else
    error('repeated entries')
end
MEGsum.FluidIQ = nan(height(MEGsum),1);
for sub_i = 1:height(MEGsum)
    [lia,locb] = ismember(MEGsum.participant_id(sub_i,:),cellstr(tasksum.CCID));
    if lia && tasksum.Ntrials(locb)==tasksum.Nexpected(locb)
        MEGsum.FluidIQ(sub_i)= tasksum.TotalScore(locb);       
    end
end

% QA:a total score of <12 would be unusual and warrant a double check of the raw data.
fprintf('IQ<12 N = %d\n',sum(MEGsum.FluidIQ<12)) %CCID = 'CC712027', 87 year-old fIQ = 11

%% TOT

task_i = 4;
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
%% VSTM
				
% VARIABLE    DESCRIPTION                           Units/Type[Range]				
% ----------  ------------------------------------  -----------------				
% K      	Number of reportable items 		  Items/Double[0:set-size]	
% 	       (set-size * probability of 			
%              reporting an item)				
% Precision   Accuracy of reported hue              degrees^-1/Double[0:Inf]				
%              (reciprocal of standard deviation 				
%              of fitted von-mises)				
% MSE	    Mean squared error			  degrees^2/Double[0 32400]
% RT	    Median reaction time		  ms/Double[0 Inf]	

% "At loads >1, the model returns the additional parameter:"				
% 				
% VARIABLE    DESCRIPTION                           Units/Type[Range]				
% ----------  ------------------------------------  -----------------				
% NonTarg     The probability of reporting a 	  probability/double[0 1]			
%              memorised but uncued item					

if exist('VSTM.mat','file')
    load VSTM.mat
else
    genVSTMsum
end

task_i = 3;
curTask = AllInfo.task{task_i};
% filename = AllInfo.sumFile{task_i};
% tasksum = tdfread(filename,'\t');
% tasksum = struct2table(tasksum);
% tasksum = sortrows(tasksum,'CCID');

fprintf('%s participants N:%d\n',curTask,height(WMs))

MEGsum.wmk = nan(height(MEGsum),4);
MEGsum.wmp = nan(height(MEGsum),4);
MEGsum.wmd = nan(height(MEGsum),4);
MEGsum.wmNT = nan(height(MEGsum),3);

MEGsum.wmk0 = nan(height(MEGsum),1);
MEGsum.wmp0 = nan(height(MEGsum),1);

for sub_i = 1:height(MEGsum)    
    [lia,locb] = ismember(MEGsum.participant_id(sub_i,:),cellstr(WMs.CCID));
    if lia
        MEGsum.VSTM{sub_i}= table2struct(tasksum(locb,:));
        MEGsum.wmp(sub_i,1) = WMs.Prcsn_ss1(locb);
        MEGsum.wmp(sub_i,2) = WMs.Prcsn_ss2(locb);
        MEGsum.wmp(sub_i,3) = WMs.Prcsn_ss3(locb);
        MEGsum.wmp(sub_i,4) = WMs.Prcsn_ss4(locb);
        MEGsum.wmp0(sub_i) = WMs.Prcsn_PerceptionTest(locb);

        MEGsum.wmk(sub_i,1) = WMs.K_ss1(locb);
        MEGsum.wmk(sub_i,2) = WMs.K_ss2(locb);
        MEGsum.wmk(sub_i,3) = WMs.K_ss3(locb);
        MEGsum.wmk(sub_i,4) = WMs.K_ss4(locb);
        MEGsum.wmk0(sub_i) = WMs.K_PerceptionTest(locb);

        MEGsum.wmd(sub_i,1) = WMs.Doubt_ss1(locb);
        MEGsum.wmd(sub_i,2) = WMs.Doubt_ss2(locb);
        MEGsum.wmd(sub_i,3) = WMs.Doubt_ss3(locb);
        MEGsum.wmd(sub_i,4) = WMs.Doubt_ss4(locb);

        MEGsum.wmNT(sub_i,1) = WMs.NonTarg_ss2(locb);
        MEGsum.wmNT(sub_i,2) = WMs.NonTarg_ss3(locb);
        MEGsum.wmNT(sub_i,3) = WMs.NonTarg_ss4(locb);
    end
end
MEGsum.wmkavg = mean(MEGsum.wmk,2,'omitnan'); 
MEGsum.wmpavg = mean(MEGsum.wmp,2,'omitnan'); 
MEGsum.wmdavg = mean(MEGsum.wmd,2,'omitnan'); 
MEGsum.wmNTavg = mean(MEGsum.wmNT,2,'omitnan'); 

%%
SetFigBasic
figure(Position=[100 100 600 250])
subplot(1,2,1);hold all
axis square;
valid_idx = ~isnan(MEGsum.TOTindx);
x = MEGsum.age(valid_idx);
y = MEGsum.TOTindx(valid_idx);
scatter(x,y);
xlabel('Age')
ylabel('TOT\_ratio')
% [rho,pval] = corr(x,y,'type','Spearman');
[rho,pval] = corr(x,y,'type','Pearson');

title(sprintf('N = %d\nrho = %.3f, p = %.3f',sum(valid_idx),rho,pval))


% [h,p,k,c] = lillietest(x)

subplot(1,2,2);hold all
axis square;
valid_idx = ~isnan(MEGsum.FluidIQ);
x = MEGsum.age(valid_idx);
y = MEGsum.FluidIQ(valid_idx);
scatter(x,y);
xlabel('Age')
ylabel('FluidIQ')
[rho,pval] = corr(x,y,'type','Spearman');
[rho,pval] = corr(x,y,'type','Pearson');

title(sprintf('N = %d\nrho = %.3f, p = %.3f',sum(valid_idx),rho,pval))
%%

figure('Position',[100 100 1500 250])
subplot(1,4,1);hold all
axis square;
valid_idx = ~isnan(MEGsum.wmkavg);
x = MEGsum.age(valid_idx);
y = MEGsum.wmkavg(valid_idx);
scatter(x,y);
xlabel('Age')
ylabel('Capacity K')
[rho,pval] = corr(x,y,'type','Spearman');
title(sprintf('N = %d\nrho = %.3f, p = %.3f',sum(valid_idx),rho,pval))

subplot(1,4,2);hold all
axis square;
valid_idx = ~isnan(MEGsum.wmpavg);
x = MEGsum.age(valid_idx);
y = MEGsum.wmpavg(valid_idx);
scatter(x,y);
xlabel('Age')
ylabel('Precision')
[rho,pval] = corr(x,y,'type','Spearman');
title(sprintf('N = %d\nrho = %.3f, p = %.3f',sum(valid_idx),rho,pval))

subplot(1,4,3);hold all
axis square;
valid_idx = ~isnan(MEGsum.wmdavg);
x = MEGsum.age(valid_idx);
y = MEGsum.wmd(valid_idx,1);
scatter(x,y);
xlabel('Age')
ylabel('Doubt')
[rho,pval] = corr(x,y,'type','Spearman');
title(sprintf('N = %d\nrho = %.3f, p = %.3f',sum(valid_idx),rho,pval))

subplot(1,4,4);hold all
axis square;
valid_idx = ~isnan(MEGsum.wmNTavg);
x = MEGsum.age(valid_idx);
y = MEGsum.wmNT(valid_idx,1);
scatter(x,y);
xlabel('Age')
ylabel('Non-target')
[rho,pval] = corr(x,y,'type','Spearman');
title(sprintf('N = %d\nrho = %.3f, p = %.3f',sum(valid_idx),rho,pval))