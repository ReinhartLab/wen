function thresh_H = hipp_threshold_ttest(dat,alpha)
% dat is a n-by-n matrix
% alpha, 0.05/0alpha
tic
thresh_H = zeros(size(dat,1),size(dat,2));
for i = 1:size(dat,1)
    dat_row = dat(i,:);
    dat_row(isnan(dat_row)) = [];

    parfor j = 1:size(dat,2)
        value = dat(i,j);
        dat_col = dat(:,j);
        dat_col(isnan(dat_col)) = [];

        if ttest(dat_row,value,"Alpha",alpha,'Tail','left') && ttest(dat_col,value,"Alpha",alpha,'Tail','left')
            thresh_H(i,j) = 1;
        end
    end
end
toc