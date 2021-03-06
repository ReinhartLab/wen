
set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 1.5, ...
'DefaultAxesXColor', 'k', ...
'DefaultAxesYColor', 'k', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 12, ...
'DefaultAxesFontName', 'Arial', ...
'DefaultaxesFontWeight', 'bold',...
'DefaultLineLineWidth', 1.5, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 12, ...
'DefaultTextFontName', 'Arial', ...
'DefaultAxesXGrid','off',...
'DefaultAxesYGrid','off',...
'DefaultAxesBox', 'off', ...
'DefaultAxesTickLength', [0.0 0.0]);
 
% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');
