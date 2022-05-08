% /home/senthilp/caesar/Power-Envelope-Correlation/cognitive/spss_generate_ww.py


%%
clear;

freqs = [2 3 4 6 8 12 16 24 32 48 64 96 128];
freqN = length(freqs);

for fi = 1:freqN
    curF = freqs(fi);
    xlsFile = fullfile('/home/senthilp/caesar/Power-Envelope-Correlation/cognitive/',['output_' num2str(curF),'.xlsx']);
    corrSubs = readtable(xlsFile);
    allData{fi} = corrSubs(:,3:end);
    variableStr = corrSubs.Properties.VariableNames;
    variableStr = variableStr(4:end);

    for t = 1:length(variableStr)
        y = table2array(corrSubs(:,variableStr{t}));
        [rho(t,fi),pval(t,fi)] = corr(corrSubs.age,y,'type','Spearman');
    end
end

freqStr = num2cell(freqs);
freqStr = cellfun(@(x)[num2str(x) 'Hz'],freqStr,'UniformOutput',false);
rho = array2table(rho,'RowNames',variableStr,'VariableNames',freqStr);
pval = array2table(pval,'RowNames',variableStr,'VariableNames',freqStr);

allData = cellfun(@table2array,allData,'UniformOutput',false);
allData = cat(3,allData{:});%subN*var*freq
allData(:,1,:) = [];%remove age
%%
edges = [18 30 40 50 60 70 80 89];
[N,edges,bin] = histcounts(corrSubs.age,edges);

clear avgInfo
for ri = 1:length(variableStr)
    for a = 1:length(edges)-1
        ageIdx = bin==a;
        tmp_freqDat = squeeze(allData(ageIdx,ri,:));%sub*freq
        for fi = 1:freqN
            avgInfo{ri,a}(fi) = datastats(tmp_freqDat(:,fi));
        end
    end
end
%%
SetFigBasic;
ageGroup = [edges(1:end-1)' edges(2:end)'-1];
ageGroup = num2str(ageGroup);
ageGroup(:,3) = '~';
ageGroup(:,4) = '';

myColors = brewermap(max(bin),'PuOr');
figure('Position',[100 100 2900 1200]);
for ri = 1:length(variableStr)
    subplot(2,4,ri);hold all;axis square
    for a = 1:length(edges)-1
        subN = sum(bin==a);
        dat = [avgInfo{ri,a}.mean];
        se = [avgInfo{ri,a}.std]./sqrt(subN);
        shadedErrorBar(1:freqN,dat,se,{'color',myColors(a,:)});        
    end
    for fi = 1:freqN
        if table2array(pval(ri,fi))<.001
            text(fi,0.002*(mod(fi,2)+1),'***','FontSize',12,'Color','r','HorizontalAlignment','center','HandleVisibility','off')
        elseif table2array(pval(ri,fi))<.01
            text(fi,0.002*(mod(fi,2)+1),'**','FontSize',12,'color','g','HorizontalAlignment','center','HandleVisibility','off')
        elseif table2array(pval(ri,fi))<.05
            text(fi,0.002*(mod(fi,2)+1),'*','FontSize',12,'HorizontalAlignment','center','HandleVisibility','off')
        end
    end
    set(gca,'xtick',1:freqN,'XTickLabel',freqStr,'YLim',[0 inf])
    title(['Interhemispheric correlation of ' variableStr{ri}],'Interpreter','none')
    legend(ageGroup)
end