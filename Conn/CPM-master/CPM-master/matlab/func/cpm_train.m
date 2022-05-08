function [r,p,pmask,mdl]=cpm_train(x,y,pthresh,type)
% Train a Connectome-based Predictive Model
% x            Predictor variable
% y            Outcome variable
% pthresh      p-value threshold for feature selection
% r            Correlations between all x and y
% p            p-value of correlations between x and y
% pmask        Mask for significant features
% mdl          Coefficient fits for linear model relating summary features to y
% type         1=postive, 2=negative,3=both,default is 3 (by ww 20220228)


% Select significant features
[r,p]=corr(x',y);
if type ==1
    pmask=(+(r>0));
elseif type ==2
    pmask=-(+(r<0));
elseif type ==3
    pmask=(+(r>0))-(+(r<0));
end

pmask=pmask.*(+(p<pthresh));
if sum(pmask)==0
    error('no enough significant edges')%(by ww 20220228)
end

% For each subject, summarize selected features
for i=1:size(x,2)
    if type ==1
        summary_feature(i)=nanmean(x(pmask>0,i));
    elseif type ==2
        summary_feature(i)=nanmean(x(pmask<0,i));
    elseif type ==3
        summary_feature(i)=nanmean(x(pmask>0,i))-nanmean(x(pmask<0,i));
    end
end

% Fit y to summary features
mdl=robustfit(summary_feature,y');

