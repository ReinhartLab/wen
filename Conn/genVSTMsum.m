fileDir = '/home/senthilp/caesar/camcan/cognitive/cc700-scored/VSTMcolour/release001/data';
allTxt = dir(fileDir);
allTxt = allTxt(~[allTxt.isdir]);

%% QUALITY CONTROL

% Suffix of output file will be _ScoredPartial if the expected number of VSTM trials (224) were not found. (Data will still be analysed although fits may be less robust.)
% Error messages will be written to final column if:
%  - VSTM blocks were not complete
%  - percetual control test was not completed
%  - mean absolute error on perceptual control test >30 degrees (the volunteer may be colour-blind).

partialID = contains({allTxt.name},'Partial');
allTxt = allTxt(~partialID);

WMs = table;
sub_i = 0;
for i = 1:length(allTxt)
    txtFile = fullfile(fileDir,allTxt(i).name);
    WMinfo = tdfread(txtFile,'\t');
    WMinfo = struct2table(WMinfo);

    if isnan(WMinfo.ErrorMessages)
        sub_i = sub_i+1;
        tmp_id = strfind(allTxt(i).name,'CC');
        WMinfo.CCID = allTxt(i).name(tmp_id:tmp_id+7);
        WMs(sub_i,:) = WMinfo;
    end
end
WMs = movevars(WMs,'CCID','Before','Prcsn_ss1');
%%

if length(unique(cellstr(WMs.CCID))) == height(WMs)
    save VSTM.mat WMs
else
    error('repeated entries')
end
