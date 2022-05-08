clear;
AllInfo(1).task = 'MEG_rest';
AllInfo(1).sumFile = '/home/senthilp/caesar/camcan/cc700/meg/pipeline/release004/BIDS_20190411/meg_rest_raw/participants.tsv';

AllInfo(2).task = 'FluidIQ';
AllInfo(2).sumFile = '/home/senthilp/caesar/camcan/cognitive/cc700-scored/Cattell/release001/summary/Cattell_summary_ww.txt';
% ww removed summary information and generated a new file ended with '_ww'

AllInfo(3).task = 'VSTM';
AllInfo(3).sumFile = '/home/senthilp/caesar/camcan/cognitive/cc700-scored/VSTMcolour/release001/summary/VSTMcolour_summary_ww.txt';
% ww removed summary information and generated a new file ended with '_ww'
% However, this file wasn't used because it doesn't have updated measures 
% see /home/senthilp/caesar/camcan/cognitive/cc700-scored/VSTMcolour/release001/README.txt

AllInfo(4).task = 'TOT';
AllInfo(4).sumFile = '/home/senthilp/caesar/camcan/cognitive/cc700-scored/TOT/release001/summary/TOT_summary.txt';

AllInfo = struct2table(AllInfo);
save('subs_info.mat','AllInfo')

