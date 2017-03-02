function ToolbarExportVectorsBtn(hTools, hLines)
% Add a button for exporting vector to a toolbar
%
% Usage:
%   ToolbarExportVectorsBtn(T), where T is the handle for the toolbar where
%   the button should be added and H is a vector of line handles that will
%   be exported to the Matlab global workspace.
%
%

[mCData, mMap] = imread(fullfile(matlabroot,'/toolbox/matlab/icons/HDF_object01.gif'));
mMap(sum(mMap, 2) == 3, :) = NaN;
uipushtool(hTools, 'cdata', ind2rgb(mCData, mMap), ...
    'TooltipString', 'Export vectors', ...
    'userdata', hLines, ...
    'ClickedCallback', @ExportLines);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExportLines(hObj, varargin)
% Export traces to a variable in the global workspace

hLines = get(hObj, 'userdata');
mX = [];
mY = [];
for h = hLines'
    mX(:, end+1) = get(h, 'xdata')';
    mY(:, end+1) = get(h, 'ydata')';
end
assignin('base', 'mETA_x', mX);
assignin('base', 'mETA_y', mY);
disp('Lines exported to base workspace as mETA_x and mETA_y. Plot with ''plot(mETA_x, mETA_y)''')

return
