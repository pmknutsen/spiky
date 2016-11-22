function varargout = spiky(varargin)
% SPIKY Analysis software for physiological and behavioral
%
% Usage:
%   < spiky >
%   Runs the Spiky GUI
%
%   < spiky(FUN) >
%   Runs the internal function FUN in Spiky. Ex:m
%     spiky, e.g. spiky('LoadTrial(''A1807004.daq'')')
%
%   Alternative syntax for running subroutines:
%     global Spiky
%     Spiky.SUB(arg1, ..., argn)
%   where SUB is the subroutine to run.
%

% Notes:
%   When adding help entries to local functions, last comment line should
%   end with a single percentage and then newline.
%

% Spiky - Spike sorting GUI for Matlab
% Copyright (C) 2005-2016 Per M Knutsen <p.m.knutsen@medisin.uio.no>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 

% Abort if Spiky is already running
hPlot = findobj('Tag', 'Spiky');
if ~isempty(hPlot)
    disp('Spiky is already running')
    return
end

clc
disp('Spiky - Spike-sorting and analysis in Matlab')
disp('Copyright (C) 2005-2016 Per M Knutsen <p.m.knutsen@medisin.uio.no>')
disp('This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are')
fprintf('are welcome to redistribute it under certain conditions; see LICENSE for details.\n\n')

% Set path to spike-sorting software
disp('Initializing paths...')
sPath = which('spiky');
sPath = sPath(1:end-7);
addpath(genpath(sPath), '-end');

% Get list of all internal sub-routines
disp('Initializing function handles...')
csStr = mlintmex('-calls', which(mfilename));
[~,~,~,~,subs] = regexp(csStr, '[S]\d* \d+ \d+ (\w+)\n');
cSubs = [subs{:}]';

% Generate function handles for all sub-routines
global Spiky;
Spiky = struct([]);
Spiky(1).main = struct([]);
Spiky(1).import = struct([]);
Spiky(1).export = struct([]);
Spiky(1).help = eval('@GetHelp');

Spiky.analysis = struct('discrete', struct([]), 'continuous', struct([]));
for i = 1:length(cSubs)
    Spiky.main(1).(cSubs{i}) = eval(['@' cSubs{i}]);
end

% Create handles to all /analysis/discrete functions
sThisPath = [sPath filesep 'analysis' filesep 'discrete' filesep];
tFiles = dir(sThisPath);
for f = 1:length(tFiles)
    if ~isempty(strfind(tFiles(f).name, '.m')) && isempty(strfind(tFiles(f).name, '.m~'))
        sName = tFiles(f).name(1:end-2);
        Spiky.analysis.discrete(1).(sName) = eval(['@' sName]);
    end
end

% Create handles to all /analysis/discrete functions
sThisPath = [sPath filesep 'analysis' filesep 'continuous' filesep];
tFiles = dir(sThisPath);
for f = 1:length(tFiles)
    if ~isempty(strfind(tFiles(f).name, '.m')) && isempty(strfind(tFiles(f).name, '.m~'))
        sName = tFiles(f).name(1:end-2);
        Spiky.analysis.continuous(1).(sName) = eval(['@' sName]);
    end
end

% Create handles to all /import
sThisPath = [sPath filesep 'import' filesep];
tFiles = dir(sThisPath);
for f = 1:length(tFiles)
    if ~isempty(strfind(tFiles(f).name, '.m')) && isempty(strfind(tFiles(f).name, '.m~'))
        sName = tFiles(f).name(1:end-2);
        Spiky.import(1).(sName) = eval(['@' sName]);
    end
end

% Create handles to all /export
sThisPath = [sPath filesep 'export' filesep];
tFiles = dir(sThisPath);
for f = 1:length(tFiles)
    if ~isempty(strfind(tFiles(f).name, '.m')) && isempty(strfind(tFiles(f).name, '.m~'))
        sName = tFiles(f).name(1:end-2);
        Spiky.export(1).(sName) = eval(['@' sName]);
    end
end

% Make Spiky available immediately in base workspace
evalin('base', 'global Spiky')

fprintf('\nSpiky API:\nTo interface with Spiky from your own code, load the global variable Spiky.\n')
fprintf('Functions are called with the syntax Spiky.[function](). Example:\n\n')
fprintf('\tglobal Spiky;Spiky.help(''ZoomReset'')\n\tSpiky.main.ZoomReset().\n\n')
fprintf('To access the documentation of API function use Spiky.help([function]).\n\n')

global g_hSpike_Viewer
g_hSpike_Viewer = figure;
hGUI = g_hSpike_Viewer;

% Set FV default settings
FV = SetFVDefaults();
SetStruct(FV);
ThemeObject([]); % set default line colors

% Set GUI window properties
vPos = get(hGUI, 'position');
ThemeObject(hGUI, 'Name', 'Spiky', 'MenuBar', 'none', 'UserData', FV, ...
    'Tag', 'Spiky', 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', ...
    'PaperPosition', [.05 .05 .9 .9], 'InvertHardcopy', 'off', ...
    'position', [vPos(1)-200 vPos(2)-50 vPos(3)+350 vPos(4)+100], ...
    'closeRequestFcn', @ExitSpiky, 'visible', 'off', ...
    'WindowButtonMotionFcn', @GUIMouseMotion );

if isfield(get(hGUI),'SizeChangedFcn')
    set(hGUI, 'SizeChangedFcn', @GUIResize )
end

movegui(hGUI, 'center')
set(hGUI, 'visible', 'on')

% Create GUI toolbar and menu
GUIRefreshToolbar();
GUIRefreshMenu();

% Set Spiky to by default NOT be in MERGE MODE
global g_bMergeMode g_bBatchMode
g_bMergeMode = 0;
g_bBatchMode = false;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GUIResize(varargin)
% Execute code when GUI is resized
%

% Refresh GUI when it is resized
% Resizing the GUI in Linux has the effect of shuffling menu items around
%GUIRefresh()

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GUIMouseMotion(varargin)
% Execute code when mouse cursor moves over GUI
%

% Create tooltip if it does not exist
hFig = GetGUIHandle();
hTxt = findobj(hFig, 'Tag', 'AxesCursorLocationTip');
if isempty(hTxt)
    hTxt = uicontrol(hFig, 'Style', 'text', 'Position', [0 0 200 15], ...
        'HorizontalAlignment', 'left', 'Tag', 'AxesCursorLocationTip', ...
        'fontsize', 8, 'backgroundcolor', get(hFig, 'color'));
end
ThemeObject(hTxt)

% Get handle of axes cursor is above
% Note that overobj() is an undocumented function and may changed in the
% future
hAx = overobj('axes');
if isempty(hAx)
    set(hTxt(1), 'string', '')
else
    mPnt = get(hAx, 'CurrentPoint');
    vXY = mPnt(1, 1:2);
    % Update tooltip
    if abs(vXY(2)) <= 1
        sYFormat = '%.2f';
    else
        sYFormat = '%.1f';
    end
    sUnit = get(get(hAx, 'ylabel'), 'userdata');
    set(hTxt(1), 'string', sprintf([' x = %.1f  y = ' sYFormat ' ' sUnit], vXY(1), vXY(2)), ...
        'backgroundcolor', get(hFig, 'color'))
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GUIRefresh(varargin)
% Refresh Spiky GUI
%
global g_bBatchMode

hGUI = GetGUIHandle();
ViewTrialData();
GUIRefreshToolbar();
GUIRefreshMenu();
g_bBatchMode = 0;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GUIRefreshToolbar()
% Refresh toolbar in Spiky GUI
% 
% This is an alias function for sp_refreshtoolbar(). Usage is identical,
% ie simply: GUIRefreshToolbar()
%
sp_refreshtoolbar();
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GUIRefreshMenu()
% Refresh menu in Spiky GUI
% 
% This is an alias function for sp_refreshmenu(). Usage is identical,
% ie simply: GUIRefreshMenu()
%
sp_refreshmenu();
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GetScriptHelp(varargin)
% Get help entry of a script file
% 
% Usage:
%   GetScriptHelp() to select script manually
%   GetScriptHelp('script.m')
%
sScriptPath = which('spiky');
sScriptPath = [sScriptPath(1:end-7) 'scripts' filesep];

if nargin == 1
    sScript = varargin{1};
else
    % Select script file manually
    [sScript, sScriptPath] = uigetfile(sprintf('%s*.m', sScriptPath), 'Select script file');
end

% Strip .m from filename
[~, sScript, ~] = fileparts(sScript);

% Get full path to script file
sScriptPath = [sScriptPath sScript '.m'];

if exist(sScriptPath, 'file')
    disp(sprintf('Script:\t%s', sScript))
    disp(help(sScriptPath))
else
    sp_disp(sprintf('No help entry for script %s', sScript));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GetHelp(varargin)
% Get help entry of local function inside spiky.m
%
if nargin < 1
    sp_disp('Which help entry do you want?');
    return;
else
    sFunction = varargin{1};
    % Remove @ sign in case a function handle was passed, eg:
    % Spiky.help(Spiky.main.MergeChannels)
    if isa(sFunction, 'function_handle')
        sFunction = func2str(sFunction);
    end
end

% Read spiky.m
sPath = which('spiky');
sSpiky = fileread(sPath);

% Find the string 'function FUNCTIONNAME'
sExpr = ['[^\n]*function ' sFunction '([^\n]*'];
[~, nStart, nEnd] = regexp(sSpiky, sExpr, 'match');
if isempty(nStart)
    sp_disp(sprintf('No help entry for %s()', sFunction));
    return
end

disp(sprintf('Function call:\t%s', sSpiky(nStart+9:nEnd)))

% Find comments that follow
sExpr = '%[\r\n|\n|\r]';
[~, nStartC, nEndC] = regexp(sSpiky(nEnd+2:nEnd+5000), sExpr, 'match');
disp(sSpiky(nEnd:(nEnd+nStartC(1))))

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CreateAnalysisMenu(hParent, sType, varargin)
% Create and return a uimenu with callbacks to continuous analysis
% functions. Use this function to create uimenus for GUIs and figures.
% 
% Usage:
%   CreateAnalysisMenu(H, TYPE)
% 
% where H is the parent of the uimenu and TYPE is the type of analysis
% functions included in the menu. Valid values of TYPE are 'continuous' or
% 'discrete'.
%

% Load list of continuous analysis functions
sSpikyPath = which('spiky');
sContPath = [sSpikyPath(1:end-7) filesep 'analysis' filesep sType filesep ];
tFiles = dir(sContPath);

% Create menu with label
if nargin > 2
    sLabel = varargin{1};
else
    sLabel = ['&' regexprep(sType, '(\<[a-z])','${upper($1)}')];
end
hMenu = uimenu('Parent', hParent, 'Label', sLabel);

for f = 1:length(tFiles)
    if ~isempty(strfind(tFiles(f).name, '.m')) && isempty(strfind(tFiles(f).name, '.m~'));
        sName = strrep(tFiles(f).name(1:end-2), '_', ' ');
        vIndx = strfind(sName, ' ');
        sName([1 vIndx+1]) = upper(sName([1 vIndx+1]));        
        uimenu('Label', sName, 'Parent', hMenu, ...
            'Callback', {@RunAnalysis, tFiles(f).name(1:end-2), sType});
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CheckUpdate(varargin)
% Compare commit SHA hash on disk with latest version online
%
hWin = msgbox('Checking for a newer version of Spiky...', 'Spiky Update');
try
    sResp = urlread('https://api.github.com/repos/pmknutsen/spiky/commits');
catch
    close(hWin)
    warndlg(lasterr, 'Spiky')
    return
end
[~,~,~,~,cSHA] = regexp(sResp, '[{"sha":"(\w+)"');
sSHA = cell2mat(cSHA{1});
close(hWin)
if strcmp(GetGitHash(), sSHA(1:10))
    msgbox('You have the latest version of Spiky.', 'Spiky Update')
else
    msgbox(sprintf('A new version of Spiky is available at:\nhttps://github.com/pmknutsen/spiky'), 'Spiky Update')
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetAmplitudeUnit(varargin)
% Set the voltage amplitude unit on a per-channel basis.
% SetAmplitudeUnit(S)
%   sets a common unit S for all channels with an undefined unit
%   (i.e. channels where the unit is assumed to be volts)
% 
% SetAmplitudeUnit(S, C)
%   sets the voltage unit S for channel C
% 
% where     S can be:
%               'V'  (volts)
%               'mV' (millivolts)
%               'uV' (microvolts)
%

if isobject(varargin{1})
    sUnit = varargin{3};
else
    sUnit = varargin{1};
end

switch sUnit
    case 'V',  nFactor = 1;
    case 'mV', nFactor = 10^3;
    case 'uV', nFactor = 10^6;
    otherwise
        sp_disp(sprintf('%s is an invalid unit', sUnit))
        return
end
tAmplitudeUnit = struct('sUnit', sUnit, 'nFactor', nFactor);

[FV, ~] = GetStruct();
bReplot = 1;
if strcmp(FV.tAmplitudeUnit.sUnit, tAmplitudeUnit.sUnit)
    bReplot = 0;
end

FV.tAmplitudeUnit = tAmplitudeUnit;
SetStruct(FV);

if bReplot
    ViewTrialData();
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PanRight(varargin)
% Pan the displayed data range towards right
%
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
vChild = findobj(get(findobj('Tag', 'Spiky'), 'children'), 'Type', 'axes'); % axes handles
vX = get(vChild(end), 'xlim');
FV.vXlim = vX+(diff(vX)/4); % pan right by 25%
SetStruct(FV); ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PanLeft(varargin)
% Pan the displayed data range towards left
%
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
vChild = findobj(get(findobj('Tag', 'Spiky') ,'children'), 'Type', 'axes'); % axes handles
vX = get(vChild(end), 'xlim');
FV.vXlim = vX-(diff(vX)/4); % pan left by 25%
SetStruct(FV); ViewTrialData();
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomReset(varargin)
% Reset all zoom levels so the entire data range is shown
%
[FV, ~] = GetStruct();
if ~isfield(FV, 'tData') return, end
FV.vXlim = [];
cFieldnames = fieldnames(FV.tYlim);
for f = 1:length(cFieldnames)
    FV.tYlim.(cFieldnames{f}) = [];
end
SetStruct(FV);
ViewTrialData();
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomIn(varargin)
% Increase horizontal zoom by 25%
%
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
vChild = findobj(get(findobj('Tag', 'Spiky'), 'children'), 'Type', 'axes'); % axes handles
vX = get(vChild(end), 'xlim');
FV.vXlim = [mean(vX)-(diff(vX)/4) mean(vX)+(diff(vX)/4)];
SetStruct(FV); ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomOut(varargin)
% Decrease horizontal zoom by 25%
%
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
vChild = findobj(get(findobj('Tag', 'Spiky') ,'children'), 'Type', 'axes'); % axes handles
vX = get(vChild(end), 'xlim');
FV.vXlim = [mean(vX)-(diff(vX)*2) mean(vX)+(diff(vX)*2)];
SetStruct(FV); ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomRange(varargin)
% Select a horizontal zoom range manually
%
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
if FV.bPanOn, pan off, end
hold on;
hLin = plot(NaN,NaN);

% If panning is currently enabled we need to turn that off in order
% to detect mouse presses in window
set(findobj(get(GetGUIHandle(), 'children'), 'type', 'axes'), 'ButtonDownFcn', '')

% Point 1
waitforbuttonpress
global pnt1
pnt1 = get(gca,'CurrentPoint'); % button down detected

% Callback that draws the line between points
set(GetGUIHandle(), 'WindowButtonMotionFcn', @ZoomRangeUpdateLine)
try
    waitforbuttonpress;
catch
    return;
end
pnt2 = get(gca,'CurrentPoint'); % button down detected

set(GetGUIHandle(), 'windowButtonMotionFcn', ''); % reset motion callback

hLine = findobj('Tag', 'SpikyXZoomIndicator');
if isempty(hLine), return, end

FV.vXlim = sort(get(hLine(3), 'xdata'));
delete(hLine)

SetStruct(FV);
ViewTrialData()
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomRangeUpdateLine(varargin)
% Callback to update display of zoom indicator during horizontal zoom
%

global pnt1
pnt2 = get(gca, 'CurrentPoint'); % button down detected
hLine = findobj('Tag', 'SpikyXZoomIndicator');
hLine = flipud(hLine);
if isempty(hLine)
    for i = 1:6
        hLine(i) = line(0,0);
    end
    vYPos = get(gca, 'ylim');
    vXPos = get(gca, 'xlim');
    nYPos = mean(vYPos);
    nH = diff(vYPos) * .2;
    % 'Shadow' line
    set(hLine(1), 'Tag', 'SpikyXZoomIndicator', 'ydata', [nYPos nYPos])
    set(hLine(2), 'Tag', 'SpikyXZoomIndicator', 'ydata', [nYPos-nH nYPos+nH], 'xdata', [pnt1(1,1) pnt1(1,1)])
    set(hLine(3), 'Tag', 'SpikyXZoomIndicator', 'ydata', [nYPos-nH nYPos+nH], 'xdata', [pnt1(1,1) pnt1(1,1)])
    % 'Highlight' line
    set(hLine(4), 'Tag', 'SpikyXZoomIndicator', 'ydata', [nYPos nYPos])
    set(hLine(5), 'Tag', 'SpikyXZoomIndicator', 'ydata', [nYPos-nH nYPos+nH], 'xdata', [pnt1(1,1) pnt1(1,1)])
    set(hLine(6), 'Tag', 'SpikyXZoomIndicator', 'ydata', [nYPos-nH nYPos+nH], 'xdata', [pnt1(1,1) pnt1(1,1)])
    set(hLine([1:3]), 'linewidth', 3, 'color', 'k')
    set(hLine([4:6]), 'linewidth', 1, 'color', 'w')
end
set(hLine([1 4]), 'xdata', [pnt1(1,1) pnt2(1,1)])
set(hLine([3 6]), 'xdata', [pnt2(1,1) pnt2(1,1)])
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PCACleaning(varargin)
% 
%
if ~IsDataLoaded, return, end
[FV, hWin] = GetStruct();
global g_bBatchMode
% Build a matrix containing the continuous data vectors of selected
% electrode. Note that the selected electrode must contain at least 3
% channels (i.e. be a tri-trode, tetrode or more. Stereotrodes and
% individual electrodes cannot be cleaned!)
%
if ~g_bBatchMode
    if isfield(FV, 'tExperimentVariables')
        if isfield(FV, 'sVariable')
            nIndx = find(strcmp({FV.tExperimentVariables.sVariable}, 'bPCACleanedChannels'));
            if ~isempty(nIndx)
                if FV.tExperimentVariables(nIndx).sValue
                    switch questdlg('This file has already been PCA cleaned. If it is necessary to redo the cleaning it is recommended you first re-load the original files from disk.', ...
                            'Spiky', 'Continue', 'Abort', 'Abort')
                        case 'Abort', return
                    end
                end
            end
        end
    end
end

csChannels = FV.csDisplayChannels;
if length(csChannels) < 3
    waitfor(warndlg('A minimum of 3 channels must be selected for PCA cleaning.'))
    return
end

% Initialize matrix to be cleaned
mRaw = []; % create matrix that should be filtered
for i = 1:length(csChannels)
    vRaw = FV.tData.(csChannels{i})';
    % Check that vector is same length as the others
    if ~isempty(mRaw)
        if length(vRaw) ~= size(mRaw, 1)
            waitfor(warndlg(['All channels must have the same vector length. Try de-selecting channels that should not be included for PCA cleaning.']))
            return
        end
    end
    mRaw = [mRaw vRaw];
end

% Initialize global variables used in pca_cleaning()
global ptspercut;
ptspercut = FV.tData.([csChannels{1} '_KHz']) * 1000; % equal to sampling frequency

% Clean
vIndx = [strfind(FV.sLoadedTrial, '/') strfind(FV.sLoadedTrial, '\')];
if isempty(vIndx)
    sp_disp(sprintf('Start PCA cleaning of %s ...', FV.sLoadedTrial))
else
    sp_disp(sprintf('Start PCA cleaning of %s ...', FV.sLoadedTrial(vIndx(end)+1:end)))
end
mClean = pca_cleaning(mRaw, 0); % run PCA cleaning
sp_disp('Done PCA cleaning!')

% Plot cleaned on top of raw data
persistent p_bNeverPlotPCA
if isempty(p_bNeverPlotPCA), p_bNeverPlotPCA = 0; end
if ~p_bNeverPlotPCA && ~g_bBatchMode
    sAns = questdlg('Do you want to plot and compare PCA cleaned traces?', 'Spiky', 'No', 'Yes', 'Never plot', 'Yes');
else
    sAns = 'No';
end
if strcmpi(sAns, 'Never plot'), p_bNeverPlotPCA = 1; end

if strcmpi(sAns, 'Yes') && ~p_bNeverPlotPCA && ~g_bBatchMode
    hFig = figure;
    ThemeObject(hFig)
    set(hFig, 'name', 'Spiky PCA Cleaning')
    
    for c = 1:size(mRaw,2)
        hAx = axes('position', [.05 (1/size(mRaw,2))*(c-1) .95 1/size(mRaw,2)]);
        ThemeObject(hAx)
        set(hAx, 'xtick', [])
        hold
        plot(mRaw(:,c), 'g')
        plot(mClean(:,c), 'r')
        ylabel(csChannels{c})
    end
    linkaxes(get(gcf,'children'), 'xy')
    uiwait(hFig)
end

if g_bBatchMode
    sAns = 'Yes';
else
    sAns = questdlg('Keep cleaned data? Keeping will overwrite the data vectors in this .spb file. Original, raw data is not over-written.', 'Spiky', 'Yes', 'No', 'Cancel', 'Yes');
end
switch sAns
    case 'Yes'
        % Replace raw data with cleaned data in FV
        for i = 1:length(csChannels)
            FV.tData.(csChannels{i}) = mClean(:,i)';
        end

        % Log this change in FV.tExperimentVariables
        if isfield(FV, 'tExperimentVariables')
            if isfield(FV, 'sVariable')
                nIndx = find(strcmp({FV.tExperimentVariables.sVariable}, 'bPCACleanedChannels'));
                if isempty(nIndx)
                    FV.tExperimentVariables(end+1).sVariable = 'bPCACleanedChannels';
                    FV.tExperimentVariables(end).sValue = '1';
                    % Save names of channels used for PCA cleaning
                    FV.tExperimentVariables(end+1).sVariable = 'sPCACleanedChannels';
                    FV.tExperimentVariables(end).sValue = mat2str(cell2mat(csChannels'));
                else
                    FV.tExperimentVariables(nIndx).sValue = '1';
                    nIndx = strcmp({FV.tExperimentVariables.sVariable}, 'sPCACleanedChannels');
                    FV.tExperimentVariables(nIndx).sValue = mat2str(cell2mat(csChannels'));
                end
            end
        else
            % TODO
        end

        SetStruct(FV)
        ViewTrialData
        % Remember task for batching
        if ~g_bBatchMode BatchRedo([], 'PCACleaning'); end
    case 'No', return
    case 'Cancel', return
end
figure(hWin)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NewExperimentVariable(sVarName, sVarValue)
% Create a new experiment variable.
% 
% Values of existing variables are replaced.
% 
% Usage:    NewExperimentVariable(VARNAME, VARVALUE)
%
[FV, ~] = GetStruct();
if ~isfield(FV, 'tExperimentVariables')
    FV(1).tExperimentVariables = struct();
end

if ~isfield(FV.tExperimentVariables, 'sVariable')
    nIndx = 1;
else
    nIndx = find(strcmp({FV.tExperimentVariables.sVariable}, sVarName));
    if isempty(nIndx)
        nIndx = length(FV.tExperimentVariables) + 1;
    end
end

FV.tExperimentVariables(nIndx).sVariable = sVarName;
FV.tExperimentVariables(nIndx).sValue = sVarValue;

SetStruct(FV);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveResults(varargin)
% Save data to an SPB file in default path
% 
% Usage:    SaveResults()
%
global g_hSpike_Viewer
if ~IsDataLoaded, return, end
sPointer = get(g_hSpike_Viewer,'Pointer');
set(g_hSpike_Viewer,'Pointer','watch')
[FV, ~] = GetStruct();
sPath = [FV.sLoadedTrial(1:end-4) '.spb'];
save(sPath, 'FV', '-v7.3')

% Update the title of the Spiky window to reflect that settings were saved
set(GetGUIHandle(), 'Name', sprintf('%s - Spiky', GetCurrentFile()))
set(GetGUIHandle(),'Pointer', sPointer)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OpenSettings(varargin)
% OpenSettings loads data and settings from an SPB file
%
persistent p_bAlwaysUseDuplicate
if isempty(p_bAlwaysUseDuplicate)
    p_bAlwaysUseDuplicate = 0;
end
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
tData = FV.tData; % original raw data and paths
sDirectory = FV.sDirectory;
sLoadedTrial = FV.sLoadedTrial;
bSelect = 0;
if ~isempty(varargin)
    if ~isstr(varargin{1})
        bSelect = 1;
        bRefresh = 1;
    else
        sPath = varargin{1};
        bRefresh = 0;
    end
end
if bSelect
    sPath = [FV.sLoadedTrial(1:end-4) '.spb'];
    [sFile, sPath] = uigetfile(sPath, 'Select data file');
    sPath = [sPath sFile];
end

try
    load(sPath, 'FV', '-MAT')
catch
    [sLastMsg, ~] = lasterr();
    % Check for 'corrupt' string in error message
    if ~isempty(strfind(sLastMsg, 'corrupt'))
        warndlg('Failed reading file. The file may be corrupt. Deleting it will remove this error.', 'Spiky')
    end
end
    
% If data already exists in tData, possibilities are;
% 1) Fields are duplicates of those on disk (filtered, inverted, PCA cleaned etc)
% 2) Additional fields to those on disk have previously been imported

% TODO:
% If fields are real duplicates, ask what to do (keep/replace)
% If fields have been imported, keep

% FV.tData  ->  saved/processed data from .spb file
% tData     ->  raw data from .daq file only (all kept by default)
% tDataTemp ->  Temporary structure we will store data to keep in

% Keep all saved data by default
tDataTemp = FV.tData;

% Find duplicate or additional fields in DAQ file
%  - Additional (imported) fields are kept by default
%  - User is asked whether to keep duplicate fields or revert to original (raw)
if isfield(FV, 'tData') % is there saved data?
    if ~isempty(FV.tData)
        sAns = '';
        cFieldnamesDAQ = fieldnames(tData);
        cFieldnamesSPB = fieldnames(FV.tData);
        for f = 1:length(cFieldnamesDAQ)
            if any(strcmp(cFieldnamesDAQ{f}, cFieldnamesSPB)) % duplicate field (ask)
                if p_bAlwaysUseDuplicate, sAns = 'Always use duplicate';
                else
                    if isempty(sAns) % ask whether to keep or replace
                        % Ask only once
                        sAns = questdlg('A duplicate of the raw data exists in an existing Spiky generated file. Do you wish to use this data, or reload all raw data from the .daq file? Note that reloading raw data will erase imported fields.', ...
                            'Duplicate data warning', 'Use duplicate', 'Always use duplicate', 'Reload data', 'Use duplicate');
                    end
                end
                switch sAns
                    case 'Use duplicate'
                        % No need to do anything if saved fields are to be
                        % kept since these are included by default above
                    case 'Always use duplicate'
                        p_bAlwaysUseDuplicate = 1;
                    case 'Reload data' % retrieve and keep original field
                        %cFieldnamesDAQ{f}
                        tDataTemp.(cFieldnamesDAQ{f}) = tData.(cFieldnamesDAQ{f});
                end
            else % Keep by default new/imported fields that do not exist in .daq file
                if isfield(FV.tData, cFieldnamesDAQ{f})
                    tDataTemp.(cFieldnamesDAQ{f}) = FV.tData.(cFieldnamesDAQ{f});
                end
            end
        end
    end
end

FV.tData = tDataTemp;

FV.sDirectory = sDirectory;
FV.sLoadedTrial = sLoadedTrial;
if bRefresh, ViewTrialData(); end
SetStruct(FV)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AutoloadNewFiles(varargin)
% ** NOT IMPLEMENTED IN THIS VERSION **
%
[FV, ~] = GetStruct();
if ~isfield(FV, 'sDirectory')
    warndlg('You must first select the directory that should be monitored in File->Open Directory', 'Spiky')
else
    map_autoload([FV.sDirectory filesep])
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetFilterChannels(varargin)
% Interactively select which channels should be filtered
%
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
global g_hSpike_Viewer g_vSelCh
% create list of selectable channels
hFig = figure;
cFields = fieldnames(FV.tData);
csDisplayChannels = {};

csChannelDescriptions = {};
for i = 1:length(cFields)
    if ~isempty(strfind(cFields{i}, '_TimeBegin'))
        if isfield(FV.tData, cFields{i}(1:end-10))
            csDisplayChannels{end+1} = cFields{i}(1:end-10);
            nIndx = find(strcmpi({FV.tChannelDescriptions.sChannel}, csDisplayChannels{end}));
            if isempty(nIndx)
                csChannelDescriptions{end+1} = '';
            else
                csChannelDescriptions{end+1} = FV.tChannelDescriptions(nIndx).sDescription;
            end
        end
    end
end
nCL = length(csDisplayChannels)*25+5;

vHW = get(g_hSpike_Viewer, 'position');
set(hFig, 'position', [vHW(1:2) 200 nCL+25], 'menu', 'none', 'Name', 'Select channels', 'NumberTitle', 'off')

% Get currently selected filter channels
if isfield(FV, 'tFilteredChannels')
    csFiltCh = {FV.tFilteredChannels.sChannel};
else
    csFiltCh = {};
end

for i = 1:length(csDisplayChannels)
    sBoxStr = sprintf('%s (%s)', csChannelDescriptions{i}, csDisplayChannels{i});
    if any(strcmp(csDisplayChannels{i}, csFiltCh))
        uicontrol(hFig, 'Style', 'checkbox', 'Position', [10 nCL 150 20], 'String', sBoxStr, 'HorizontalAlignment', 'left', 'backgroundcolor', get(hFig, 'color'), 'value', 1);
    else
        uicontrol(hFig, 'Style', 'checkbox', 'Position', [10 nCL 150 20], 'String', sBoxStr, 'HorizontalAlignment', 'left', 'backgroundcolor', get(hFig, 'color'), 'value', 0);
    end
    nCL = nCL - 25;
end
uicontrol(hFig, 'Style', 'pushbutton', 'Position', [50 nCL 50 20], ...
    'Callback', 'global g_vSelCh; g_vSelCh=flipud(get(get(gcf,''children''),''value'')); g_vSelCh=[g_vSelCh{:}];close(gcf)', ...
    'String', 'OK' ); % OK button
uiwait % wait for user to close window or click OK button
g_vSelCh = logical(g_vSelCh(1:end-1));

% Create the tFilteredChannels structure (will replace old)
tFilteredChannels = struct('sChannel', [], 'vBandpass', [], 'bRectify', []);
csCh = csDisplayChannels(g_vSelCh);
for i = 1:length(csCh)
    tFilteredChannels(i).sChannel = csCh{i};
    tFilteredChannels(i).vBandpass = [NaN NaN]; % there is no default
    tFilteredChannels(i).bRectify = 0;
    if isfield(FV, 'tFilteredChannels')
        if isfield(FV.tFilteredChannels, 'sChannel')
            iCh = strcmp(csCh{i}, {FV.tFilteredChannels.sChannel});
            if any(iCh)
                if isfield(FV.tFilteredChannels(iCh), 'vBandpass')
                    if ~isempty(FV.tFilteredChannels(iCh).vBandpass)
                        tFilteredChannels(i).vBandpass = FV.tFilteredChannels(iCh).vBandpass;
                    end
                end
                if isfield(FV.tFilteredChannels(iCh), 'bRectify')
                    if ~isempty(FV.tFilteredChannels(iCh).bRectify)
                        tFilteredChannels(i).bRectify = FV.tFilteredChannels(iCh).bRectify;
                    end
                end
            end
        end
    end
end
FV.tFilteredChannels = tFilteredChannels;

% Force all filtered channels to be displayed
FV.csDisplayChannels = unique([FV.csDisplayChannels csCh]);
clear global g_vSelCh
SetStruct(FV)

% Display filters options dialog of any of the settings are NaNs
if any( isnan([FV.tFilteredChannels.vBandpass FV.tFilteredChannels.bRectify]) )
    FilterOptions(varargin)
else
    ViewTrialData
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FilterOptions(varargin)
% Set the low and high pass frequencies of filtered channels
% 
% FilterOptions() opens a table to enter parameters manually
%
[FV, ~] = GetStruct();
if ~isfield(FV, 'tFilteredChannels'); return; end
if isempty(FV.tFilteredChannels(1).sChannel)
    waitfor(warndlg('You first need to select the channels that should be filtered.', 'Spiky'))
    return
end

% Initialize figure
hFig = figure('closeRequestFcn', 'set(gcbf,''userdata'',1)', ...
    'Name', 'Spiky Filtered Channels', 'ToolBar', 'none', ...
    'menuBar','none', 'visible', 'off');
ThemeObject(hFig)
centerfig(hFig, GetGUIHandle());
set(hFig, 'visible', 'on')
cColumnNames = {'Channel', 'Highpass (Hz)' 'Lowpass (Hz)' 'Rectify (1/0)'};
cColumnFormat = {{'numeric' 'Adjustable'}, {'numeric' 'Adjustable'}, {'numeric' 'Adjustable'}, {'numeric' 'Adjustable'}};
cColumnEditable =  [false true true true];

% Create data cell for table
cData = {};
for i = 1:length(FV.tFilteredChannels)
    cData{i, 1} = GetChannelDescription(FV.tFilteredChannels(i).sChannel);
    if isempty(cData{i, 1})
        cData{i, 1} = FV.tFilteredChannels(i).sChannel;
    end
    cData{i, 2} = FV.tFilteredChannels(i).vBandpass(1);
    cData{i, 3} = FV.tFilteredChannels(i).vBandpass(2);
    cData{i, 4} = FV.tFilteredChannels(i).bRectify;
end

% Create table
hTable = uitable('Units', 'normalized','Position', [0 0 1 1], 'Data', cData, ...
    'ColumnName', cColumnNames, 'ColumnEditable', cColumnEditable, ...
    'ColumnWidth', {150, 100, 100 100});

% Wait for figure to be closed
waitfor(hFig, 'userdata')
cData = get(hTable, 'data');
delete(hFig)

% Assign new values to structure
for i = 1:length(FV.tFilteredChannels)
    FV.tFilteredChannels(i).vBandpass(1) = cData{i, 2};
    FV.tFilteredChannels(i).vBandpass(2) = cData{i, 3};
    FV.tFilteredChannels(i).bRectify = cData{i, 4};
end

SetStruct(FV)
ViewTrialData()
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetChannelCalculator(varargin)
% 
% 

if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
[sCh, ~] = SelectChannelNumber(FV.csChannels);

if ~isfield(FV, 'tChannelCalculator')
    FV.tChannelCalculator = struct([]);
    sEval = '';
else
    if isfield(FV.tChannelCalculator, sCh)
        sEval = FV.tChannelCalculator.(sCh);
    else
        sEval = '';
    end
end
cAnswer = inputdlg(['Enter a string that will be evaluated every time this channel is used. ' ...
                    'Note that not all operations currently support channel calculations. ' ...
                    'Substitute signal vector with ''a'' in your string:'], ...
                    'Channel Calculator', 3, {sEval});
if isempty(cAnswer), return; end % cancel button pressed
FV.tChannelCalculator(1).(sCh) = cAnswer{1};
SetStruct(FV); ViewTrialData();
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vCont = ChannelCalculator(vCont, sCh)
% Apply custom calculation to channel data.
% 
% Usage:
%   ChannelCalculator(V, C) where V is the vector data and C is the channel
%   name.
%
[FV, ~] = GetStruct();

if ~isfield(FV, 'tChannelCalculator'), return; end
if ~isfield(FV.tChannelCalculator, sCh), return; end

% Re-calculate
sEval = FV.tChannelCalculator.(sCh);
if isempty(sEval); return; end
a = vCont;

try
    if ~isempty(strfind(sEval, '='))
        sEval = [sEval ';']; % suppress output
        eval(sEval);
    else
        a = eval(sEval);
    end
catch
    sErr = sprintf('An error occurred when evaluating the expression ''%s'' for channel %s', sEval, sCh);
    waitfor(warndlg(sErr, 'Spiky'))
    return
end

if length(a) ~= length(vCont)
    sErr = sprintf('An error occurred when evaluating the expression ''%s'' for channel %s:\nSize of matrix after evaluation must be same as before.', sEval, sCh);
    waitfor(warndlg(sErr, 'Spiky'))
    return
else
    vCont = a;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetDirectory(varargin)
% Select a new path interactively as the Spiky working directory. All DAQ
% files contained in the selected directory are added to the Files list in
% the File menu.
% 
% Usage: SetDirectory()
%

FV = SetFVDefaults(); % initialize FV with default settings
global g_bBatchMode

% Path to session directory (without last slash at end)
if isfield(FV, 'sDirectory')
    sDirectory = uigetdir([FV.sDirectory filesep]);
else
    sDirectory = uigetdir;
end
if sDirectory == 0, return, end
cd(sDirectory)

% Get list of DAQ files
sFiles = dir(sprintf('%s%s*.daq', sDirectory, filesep));
if isempty(sFiles)
    warndlg('No DAQ files were found in the chosen directory.', 'No files found')
    return
end
hFileList = findobj(GetGUIHandle(), 'Tag', 'MenuFileList');

% Clear current file list (children of MenuFileList)
hFileListChildren = findobj(GetGUIHandle(), 'Parent', hFileList);
delete(hFileListChildren)
for f = 1:size(sFiles, 1)
    uimenu(GetGUIHandle(), 'Parent', hFileList, 'Label', sFiles(f).name, ...
        'UserData', [sDirectory filesep sFiles(f).name], ...
        'Callback', sprintf('Spiky.main.OpenFile([], %d);', f) );
end

% Update FV structure and refresh GUI
FV.sDirectory = sDirectory;
SetStruct(FV)
OpenFile([], 1);
g_bBatchMode = false;
ViewTrialData();

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetDirectoryTree(varargin)
% Load list of all DAQ files from all descending directories below a top directory
%

% Set the top directory of the tree
FV = SetFVDefaults();

% Path to session directory (without last slash at end)
if isfield(FV, 'sDirectory')
    sDirectory = uigetdir([FV.sDirectory filesep]);
else
    sDirectory = uigetdir;
end
if sDirectory == 0, return, end

% Iterate recursively over all sub-directories and extract all .DAQ files
[cPaths, cFiles] = GetFilePaths(sDirectory, '.daq');

% Display error if no files were found
if isempty(cFiles)
    warndlg('No DAQ files were found in the chosen directory.', 'Spiky')
    return
end

% Reset current file list
hFileList = findobj(GetGUIHandle(), 'Tag', 'MenuFileList');
hFileListChildren = findobj(GetGUIHandle(), 'Parent', hFileList);
delete(hFileListChildren)

% Update GUI
for f = 1:length(cFiles)
    uimenu(GetGUIHandle(), 'Parent', hFileList, 'Label', cFiles{f}, ...
        'UserData', [cPaths{f} filesep cFiles{f}], ...
        'Callback', sprintf('Spiky.main.OpenFile([], %d);', f) );
end

% Update FV structure and refresh GUI
FV.sDirectory = sDirectory;
SetStruct(FV)
OpenFile([],1);
ViewTrialData();

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function csExt = GetImportFilters(varargin)
% Get list of import filters
%
% Usage:
%   GetImportFilters(), returns list of filters
%   GetImportFilters('*.mat'), returns filter list with '*.mat' at top
%   GetImportFilters(DIR) returns filter list with an extension found in
%   the directory DIR placed at top
% 
% The second calling syntax can be used to select a default filter. The
% third syntax is useful when default should match a file type on disk.
%
if ~isempty(varargin)
    if exist(varargin{1}, 'dir')
        % Choose default file types from selected directory
        tDir = dir(varargin{1});
        cExt = cell(length(tDir), 1);
        for f = 1:length(tDir)
            [~, ~, cExt{f}] = fileparts(tDir(f).name);
            cExt{f} = ['*', cExt{f}];
        end
        cDef = unique(cExt);
    else
        sDef = varargin{1}; % default from input
    end
else
    sDef = []; % no default
end

sDir = which(mfilename);
tDir = dir([sDir(1:end-7) 'import']);
csExt = {};
for i = 3:length(tDir)
    if ~isempty(strfind(tDir(i).name, '~'))
        continue % ignore backup files
    end
    csExt{end+1,1} = ['*.' tDir(i).name(8:end-2)];
    sDescr = eval(['help(''' tDir(i).name ''')']);
    csExt{end,2} = [sDescr(2:end-1) ' (' csExt{end,1} ')'];
end

% Remove invalid defaults
if exist('cDef', 'var')
    iRem = [];
    for e = 1:length(cDef)
        if ~any(strcmp(csExt(:, 1), cDef{e}))
            iRem(end+1) = e;
        end
    end
    cDef(iRem) = [];
    if ~isempty(cDef)
        sDef = cDef{1};
    end
end

% Re-order if a default filter was specified
if ~isempty(sDef)
    iDef = strcmp(csExt(:, 1), sDef);
    if any(iDef)
        csExt = [csExt(iDef, :); csExt(~iDef, :)];
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function csExt = GetExportFilters()
% Get list of import filters.
%
sDir = which(mfilename);
tDir = dir([sDir(1:end-7) 'export']);
csExt = {};
for i = 3:length(tDir)
    csExt{end+1,1} = ['*.' tDir(i).name(8:end-2)];
    sDescr = eval(['help(''' tDir(i).name ''')']);
    csExt{end,2} = [sDescr(2:end-1) ' (' csExt{end,1} ')'];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cPaths, cFiles] = GetFilePaths(sBaseDir, sSuffix)
% Get paths of videos in a directory and its sub-directories
% Recursive search for .bin files from a selected directory
% inputs: sSuffix    Find files this file extension (e.g. '.avi')
%         sBaseDir   Base directory
%
cPaths = {};
cFiles = {};
tDirList = dir(sBaseDir);
for t1 = 3:length(tDirList)
    if tDirList(t1).isdir % depth = 1
        % recursively call this function
        sBaseDirRecurs = [sBaseDir filesep tDirList(t1).name];
        [cPaths2, cFiles2] = GetFilePaths(sBaseDirRecurs, sSuffix);
        cPaths = {cPaths{:}, cPaths2{:}};
        cFiles = {cFiles{:}, cFiles2{:}};
    else
        % Add file names to paths cell
        if strcmp(tDirList(t1).name(end-3:end), sSuffix)
            cPaths{end+1} = sBaseDir;
            cFiles{length(cPaths)} = tDirList(t1).name;
        end
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vTime = GetTime(nBeginTime, nEndTime, vCont, nFs)
% Get time vectors for continuous 1- or 2-D data
% 
% Usage:
%   GetTime(BEGIN, END, CONT, FS)
%
%   where BEGIN is the start time, END is the end time (can be empty), CONT
%   is the continuous signal (only its length is used) and FS is the
%   sampling rate.
%
if all(size(vCont) > 1) % 2D trace (e.g. spectrogram)
    nEndTime = nBeginTime + (size(vCont, 2))*(1/nFs);
    vTime = linspace(nBeginTime, nEndTime, size(vCont, 2));
else % 1D trace (voltage etc)
    % TODO: Why don't I always calculate vTime from nFs
    % Compute end time if it is not provided
    if isempty(nEndTime)
        %length(vCont) nFs
    end
    vTime = linspace(nBeginTime, nEndTime, length(vCont));
    % Check that time intervals == nFs; if not, then recalculate
    % vTime from nFs
    nFs_i = 1 / diff(vTime(1:2));
    if nFs_i ~= nFs
        vTime = nBeginTime:(1/nFs):(nBeginTime+(1/nFs)*length(vCont)+1);
        vTime = vTime(1:length(vCont));
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bResult = IsMergeMode(varargin)
% Determine merge mode. Return 1 if Spiky is currently running in merge
% mode. Otherwise, returns 0.
%

% Get global merge mode status (generally correct)
global g_bMergeMode

% Assume we're in merge mode if the loaded file is called MergeFile
[FV, ~] = GetStruct();
[~, sFile] = fileparts(FV.sLoadedTrial);
if strcmpi(sFile, 'MergeFile')
    g_bMergeMode = 1;
end

if isempty(g_bMergeMode)
    g_bMergeMode = 0;
end
bResult = g_bMergeMode;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AttachAxisCrossHairs(varargin)
% Attach crosshairs callback to axis.
%
if nargin < 1
    sp_disp('No axis handle provided.');
    return;
end
if ~strcmp(get(varargin{1}, 'type'), 'axes')
    sp_disp('Invalid handle. Must be axis.');
    return;
end
hAx = varargin{1};
set(get(hAx, 'parent'), 'WindowButtonMotionFcn', @UpdateAxisCrossHairs);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateAxisCrossHairs(varargin)
% Draw crosshairs on current axis.
%

% Get axis handle
% Note: overobj() is an undocumented function
hAx = overobj('axes');
if isempty(hAx), return; end

% Get current axis coordinates
mPnt = get(hAx, 'CurrentPoint');
vXY = mPnt(1, 1:2);

% Get axis limits
vXlim = get(hAx, 'xlim');
vYlim = get(hAx, 'ylim');

% Create or update crosshairs
sTag = 'SpikyCrosshairs';
hX = findobj(hAx, 'tag', sTag);
bAxisHeld = ishold(hAx);
if isempty(hX)
    if ~bAxisHeld, hold(hAx, 'on'); end
    hX = plot(hAx, vXlim, vXY([2 2]), vXY([1 1]), vYlim);
    hX(end+1) = plot(hAx, vXY(1), vXY(2), 'ro');
    if ~bAxisHeld, hold(hAx, 'off'); end
    ThemeObject(hX(1:2), 'linestyle', '-.', 'tag', sTag);
    set(hX(3), 'color', 'r', 'tag', sprintf('%sMarker', sTag));
else
    set(hX(1), 'xdata', vXlim, 'ydata', vXY([2 2]));
    set(hX(2), 'xdata', vXY([1 1]), 'ydata', vYlim);
    hMarker = findobj(hAx, 'tag', sprintf('%sMarker', sTag));
    set(hMarker, 'xdata', vXY(1), 'ydata', vXY(2));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ViewTrialData(varargin)
% Update main window in GUI (channels and events) with current data and
% options.
%

% Make Spiky window current without raising it to the top
hFig = findobj('Tag', 'Spiky');
set(0, 'currentfigure', hFig)
ThemeObject(hFig);

if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();

% Destroy waitbar unless we are in Merge or Batch mode
global g_bBatchMode
if ~IsMergeMode() && ~g_bBatchMode
    SpikyWaitbar(1,1);
end

% Destroy waitbar if ViewTrialData was evoked from the GUI
if ~isempty(varargin)
    if ishandle(varargin{1})
        SpikyWaitbar(1,1);
    end
end

% Delete all current axes objects
delete(findobj(hFig, 'type', 'axes'))

% Delete all current uicontextmenu objects
delete(findobj(hFig, 'type', 'uicontextmenu'))

% Reset click callback (unless we're in pan mode)
hPan = pan(hFig);
if strcmp(get(hPan, 'Enable'), 'off')
    set(hFig, 'WindowButtonDownFcn', '');
end
zoom(hFig, 'off')

% Select channels
if isempty(FV.csDisplayChannels) && ~IsMergeMode(), SelectChannels; end
[FV, hWin] = GetStruct();

% Check that all channel names exist in FV.csChannels
cFields = fieldnames(FV.tData);
for i = 1:length(cFields)
    if ~isempty(strfind(cFields{i}, '_TimeBegin'))
        FV.csChannels{end+1} = cFields{i}(1:end-10);
    end
end
FV.csChannels = unique(FV.csChannels);

% Update y-axis limits for all channels
for i = 1:length(FV.csChannels)
    if isfield(FV.tData, FV.csChannels{i})
        vMinMax = [min(FV.tData.(FV.csChannels{i})) max(FV.tData.(FV.csChannels{i}))];
        [vMinMax, FV] = AdjustChannelGain(FV, vMinMax, FV.csChannels{i});
        if ~isfield(FV.tYlim, FV.csChannels{i})
            FV.tYlim(1).(FV.csChannels{i}) = vMinMax;
        end
    else
        if ~isfield(FV.tYlim, FV.csChannels{i})
            FV.tYlim(1).(FV.csChannels{i}) = []; % assume its digital
        end
    end
end

% Check that all DisplayChannels exist (and remove those that dont)
vRemIndx = [];
for i = 1:length(FV.csDisplayChannels)
    if ~isfield(FV.tData, FV.csDisplayChannels{i})
        vRemIndx = [vRemIndx i];
    end
end
FV.csDisplayChannels(vRemIndx) = [];
SetStruct(FV, 'nosaveflag')

hSubplots = [];
bShowDigitalEvents = strcmp(get(findobj(hWin, 'Label', 'Show &Events'),'checked'), 'on');
if isempty(FV.csDigitalChannels), bShowDigitalEvents = 0; end

if bShowDigitalEvents
    nSubEventHeight = 0.15; % event axis is fixed height
else nSubEventHeight = 0; end

% Iterate over channels alphabetically
for i = 1:length(FV.csDisplayChannels)
    sCh = FV.csDisplayChannels{end - i + 1};
    vCont = ChannelCalculator(FV.tData.(sCh), sCh); % continuous trace (V)
    if length(vCont) <= 1 && ~FV.bPlotRasters
        uiwait(warndlg(sprintf('Cannot display channel %s as it is not a vector.', sCh)))
        continue
    end
    
    % Get channel description
    sYLabel = GetChannelDescription(sCh);
    if isempty(sYLabel), sYLabel = sCh; end

    if ( ~FV.bPlotRasters || ~isfield(FV.tSpikes, sCh) ) || ~IsMergeMode()
        [vCont, ~] = AdjustChannelGain(FV, vCont, sCh);
        nContAllMinMax = [min(vCont(:)) max(vCont(:))];
        [FV, ~] = GetStruct(); % reload FV since channel gains may have changed in previous line
        
        nFs = FV.tData.([sCh '_KHz']) * 1000; % sampling frequency (Hz)
        nBeginTime = FV.tData.([sCh '_TimeBegin']); % start of sampling (sec)
        sEndField = [sCh '_TimeEnd'];
        if isfield(FV.tData, sEndField)
            nEndTime = FV.tData.(sEndField); % start of sampling (sec)
        else
            nEndTime = [];
        end
        vTime = GetTime(nBeginTime, nEndTime, vCont, nFs);

        % Limit vector to FV.vXlim
        if ~isempty(FV.vXlim)
            vXlim(1) = max([1 round((FV.vXlim(1)-nBeginTime)/diff(vTime(1:2)))]);
            vXlim(2) = min([length(vTime) round((FV.vXlim(2)-nBeginTime)/diff(vTime(1:2)))]);
            if ~(vXlim(2) <= vXlim(1))
                vTime = vTime(vXlim(1):vXlim(2));
                if all(size(vCont) > 1) % 2D
                    vCont = vCont(:, vXlim(1):vXlim(2));
                else % 1D
                    vCont = vCont(vXlim(1):vXlim(2));
                end
            end
        end
    end
    
    % Custom filter and decimate continuous channels
    if ~FV.bPlotRasters || ~isfield(FV.tSpikes, sCh)
        % Custom filter
        if isfield(FV, 'tFilteredChannels')
            iParams = strcmp(sCh, {FV.tFilteredChannels.sChannel});
            if any(iParams)
                nHiPass = FV.tFilteredChannels(iParams).vBandpass(1);
                nLoPass = FV.tFilteredChannels(iParams).vBandpass(2);
                bRectify = FV.tFilteredChannels(iParams).bRectify;                
                if nFs < 2000 || nLoPass < 2000
                    [vCont, vTime, ~] = FilterChannel(vCont, vTime, nFs, nLoPass, nHiPass, bRectify, 'none');
                else
                    [vCont, vTime, ~] = FilterChannel(vCont, vTime, nFs, nLoPass, nHiPass, bRectify, 'decimate');
                end
            end
            FV.tYlim.(sCh) = [min(vCont) max(vCont)];
        end
        % Decimate
        nLen = 20000;
        vShowIndx = unique(round(linspace(1,length(vTime),nLen)));
        
        if any(size(vCont) == 1)
            vTimeAllTrace = vTime;
            vTime = vTime(vShowIndx);
            vContAllTrace = vCont;
            vCont = vCont(vShowIndx);
        end
    end

    % Continuous trace
    nSubs = length(FV.csDisplayChannels);

    % Create axis
    nYBase = .06; % lower boundary of bottom axis
    nXBase = .07; % left-most boundary axes
    nYSep = 0.008; % separation between events and top axis
    nSubHeight = (((1-nYBase-nYSep)-(nSubs*.02+nSubEventHeight)) / nSubs);
    nSubY = nYBase + (.02*(i-1)) + (nSubHeight*(i-1));
    hSubplots(end+1) = axes('position', [nXBase nSubY .91 (nSubHeight)+.02], ...
        'tag', sCh, ...
        'nextplot', 'ReplaceChildren');
    
    % Plot spike rasters
    bShowSpikes = 0;
    if FV.bPlotRasters && isfield(FV.tSpikes, sCh)
        if ~isempty(FV.tSpikes.(sCh).waveforms)
            bShowSpikes = 1;
        end
    end
    
    if bShowSpikes
        tSpikes = FV.tSpikes.(sCh);
        nFs = tSpikes.Fs(1);
        nRow = 0;
        csUnits = {};

        if ~isfield(tSpikes, 'hierarchy'), vUnits = 0; % units not sorted
        else vUnits = unique(tSpikes.hierarchy.assigns); end % units are sorted

        for nU = 1:length(vUnits) % iterate over units
            if ~isfield(tSpikes, 'hierarchy')
                vIndx = 1:length(tSpikes.spiketimes);
            else
                vIndx = find(tSpikes.hierarchy.assigns == vUnits(nU));
            end
            vSpiketimes = tSpikes.spiketimes(vIndx); % samples
            vSpiketimes = vSpiketimes ./ nFs; % sec

            % Drop spikes that violate the deadtime parameter
            vSpiketimes = DropDeadtimeSpikes(vSpiketimes);
            
            % Drop spikes outside of visible range (+/- 50%)
            if ~isempty(FV.vXlim)
                vIndx = vSpiketimes < [min(FV.vXlim)] | vSpiketimes > [max(FV.vXlim)];
                vSpiketimes(vIndx) = [];
            end            
            if vUnits(nU) == 0, mCol = [.4 .4 .4];
            else mCol = FV.mColors(nU,:); end

            nRow = nRow + 1;

            if ~isfield(FV, 'bPlotInstSpikeRate'), FV.bPlotInstSpikeRate = 0; end
            if FV.bPlotInstSpikeRate
                % Round spiketimes to nearest ms
                if length(vSpiketimes) > 1
                    nStartSec = vSpiketimes(1);
                    vSpiketimes = vSpiketimes - nStartSec;
                    vSpiketimes = round(vSpiketimes .* 1000); % ms
                    
                    % Create impulse vector from spike times
                    [vN, vX] = hist(vSpiketimes, 1:vSpiketimes(end));
                    vImpVec = zeros(vSpiketimes(end), 1);
                    vImpVec(vX) = vN;
                    
                    nLen = 10; % base width 5 ms
                    vWin = [linspace(1,0,nLen)];
                    vWin = vWin/sum(vWin); % normalize to area 1
                    vC = conv(vImpVec, vWin);
                    vC = vC(1:end-nLen);
                    vT = (1:vSpiketimes(end)) ./ 1000 + nStartSec;
                    vT = vT(1:end-1);
                else vC = []; vT = []; end
                
                hLin = plot(hSubplots(end), vT, vC, 'color', mCol);
                mUserData = [vT' vC];
            else
                % Plot spike rasters
                
                % Add 'thickness' to spikes (1/1000 of the window width)
                if isempty(vSpiketimes) continue; end
                if isempty(FV.vXlim) nThickTime = diff(vSpiketimes([1 end])) / 2500;
                else nThickTime = diff(FV.vXlim) / 2500; end
                vX1 = vSpiketimes' - nThickTime;
                vX2 = vSpiketimes' + nThickTime;
                if isempty(vX1), continue, end
                
                % Set spike width
                vIndices = [vX1' vX2' vX2' [vX1(2:1:end) NaN]'];
                vIndices = reshape(vIndices', numel(vIndices), 1);
                vXlim = FV.vXlim;
                if isempty(FV.vXlim)
                    vXlim = [min(vSpiketimes) max(vSpiketimes)];
                end
                nSpikeWidth = GetEventWidth(vXlim);
                vIndices(2:2:end) = vIndices(1:2:end) + nSpikeWidth;
                
                % Rearrange indices
                vYs = repmat([nRow nRow NaN NaN], 1, length(vIndices)/4);
                
                % Get height of event markers
                set(hSubplots(end), 'units', 'pixels')
                vPos = get(hSubplots(end), 'position');
                nSpikeHeight = vPos(4) / length(vUnits) / 2;
                set(hSubplots(end), 'units', 'normalized')
                hLin = plot(vIndices, vYs, 'linewidth', nSpikeHeight, 'color', mCol);
                mUserData = [vIndices vYs'];
            end
            
            hold on
            hMenu = uicontextmenu;
            set(hMenu, 'userdata', mUserData, 'Tag', [sCh '-' num2str(vUnits(nU))]);
            set(hMenu, 'Tag', [sCh '-' num2str(vUnits(nU))]);
            uimenu(hMenu, 'Label', '&Copy to Figure', 'Callback', @CopyChannelToFig, 'Tag', sCh);
            uimenu(hMenu, 'Label', '&Hide', 'Callback', @SelectChannels, 'Tag', sCh);
            uimenu(hMenu, 'Label', '&Delete Section', 'Callback', @DeleteSection, 'Tag', sCh, 'Separator', 'on');
            uimenu(hMenu, 'Label', '&Shift Time', 'Callback', @ShiftBeginTime, 'Tag', sCh);
            uimenu(hMenu, 'Label', '&Replace Values', 'Callback', @ReplaceValues, 'Tag', sCh);
            uimenu(hMenu, 'Label', 'Show &Markers', 'Callback', 'set( get(gcbo, ''userdata''), ''marker'', ''o'')', 'userdata', hLin, 'Separator', 'on');
            uimenu(hMenu, 'Label', 'Meas&ure', 'Callback', @MeasureLine, 'Tag', sCh);
            uimenu(hMenu, 'Label', 'Digiti&ze Channel (Manual)', 'Callback', @DigitizeChannel, 'Tag', sCh);
            uimenu(hMenu, 'Label', '&Properties', 'Separator', 'on', 'Callback', @ChannelProperties, 'Tag', sCh);
            uimenu(hMenu, 'Label', '&Set Electrode Position', 'Callback', @SetElectrodePosition, 'Tag', sCh);
            uimenu(hMenu, 'Label', '&Set Receptive Field', 'Callback', @SetReceptiveField, 'Tag', sCh);
            set(hLin, 'uicontextmenu', hMenu)
            
            csUnits{end+1} = num2str(vUnits(nU));
        end
    else
        if length(find(size(vCont)>1)) == 1
            % Plot a continuous trace
            if ~exist('nContTrace', 'var')
                nContTrace = 1;
            end
            nContTrace = nContTrace + 1;
            hLin = plot(hSubplots(end), vTime, vCont, 'color', FV.mColors(nContTrace,:));
            vY = [(floor(min(vCont) / 10) * 10 - 5) (ceil(max(vCont) / 10) * 10 + 5)];
        else
            % Plot a 3D matrix (e.g. spectrogram)
            vY = [];
            if length(vY) ~= size(vCont, 1)
                vY = 1:size(vCont, 1);
            end
            if isfield(FV.tData, [sCh '_Scale'])
                vY = FV.tData.([sCh '_Scale']);
            end
            
            % Insert the min/max values of the entire matrix into displayed
            % segment so that colormap is scaled identically regardless of
            % zoom level.
            vCont(1,1) = nContAllMinMax(1);
            vCont(end,end) = nContAllMinMax(2);
            hLin = imagesc(vTime, vY, real(vCont), 'parent', hSubplots(end));
            set(hSubplots(end), 'ydir', 'normal')
        end
        % attach context menu to line to enable additional options
        hMenu = uicontextmenu;
        if exist('vTimeAllTrace')
            set(hMenu, 'userdata', [vTimeAllTrace(:) vContAllTrace(:)], 'Tag', sCh);
        end
        set(hMenu, 'Tag', sCh);
        CreateAnalysisMenu(hMenu, 'continuous', 'Analysis')
        uimenu(hMenu, 'Label', '&Copy to Figure', 'Callback', @CopyChannelToFig, 'Tag', sCh, 'Separator', 'on');
        uimenu(hMenu, 'Label', '&Hide', 'Callback', @SelectChannels, 'Tag', sCh);
        uimenu(hMenu, 'Label', 'Move Up', 'Callback', {@ChangeDisplayOrder, sCh, 1});
        uimenu(hMenu, 'Label', 'Move Down', 'Callback', {@ChangeDisplayOrder, sCh, -1});
        uimenu(hMenu, 'Label', '&Delete Section', 'Callback', @DeleteSection, 'Tag', sCh, 'Separator', 'on');
        uimenu(hMenu, 'Label', '&Shift Time', 'Callback', @ShiftBeginTime, 'Tag', sCh);
        uimenu(hMenu, 'Label', '&Replace Values', 'Callback', @ReplaceValues, 'Tag', sCh);
        uimenu(hMenu, 'Label', 'Show &Markers', 'Callback', 'set( get(gcbo, ''userdata''), ''marker'', ''o'')', 'userdata', hLin, 'Separator', 'on');
        uimenu(hMenu, 'Label', 'Meas&ure', 'Callback', @MeasureLine, 'Tag', sCh);
        uimenu(hMenu, 'Label', 'Digiti&ze Channel', 'Callback', @DigitizeChannelAuto, 'Tag', sCh, 'Separator', 'on');
        uimenu(hMenu, 'Label', 'Digitize Manually', 'Callback', @DigitizeChannel, 'Tag', sCh);
        uimenu(hMenu, 'Label', 'Digitize Crossing', 'Callback', @DigitizeChannelCrossing, 'Tag', sCh);
        uimenu(hMenu, 'Label', '&Properties', 'Separator', 'on', 'Callback', @ChannelProperties, 'Tag', sCh);
        uimenu(hMenu, 'Label', 'Delete Channel', 'Callback', {@DeleteChannel, sCh}, 'Tag', sCh, 'Separator', 'on');
        set(hLin, 'uicontextmenu', hMenu)
    end
    
    % Set axis labels
    if i == 1, xlabel(hSubplots(end), 'Time (sec)'); end
    if isfield(FV.tData, [sCh '_Unit'])
        sUnit = FV.tData.([sCh '_Unit']);
    else
        sUnit = FV.tAmplitudeUnit.sUnit;
    end
    
    if FV.bPlotRasters && isfield(FV.tSpikes, sCh)
        hLabel = ylabel(hSubplots(end), hSubplots(end), sYLabel);
    else
        if isempty(sUnit)
            hLabel = ylabel(sprintf('%s', strrep(sYLabel, '_', ' ')));
        else
            hLabel = ylabel(sprintf('%s', strrep(sYLabel, '_', ' ')));
        end
        set(hLabel, 'userdata', sUnit, ...
            'units', 'normalized', ...
            'position', [-.03 0.5 0]);
    end
    set(hLabel, 'Interpreter', 'tex')
    ThemeObject(hSubplots(end))
    
    if IsMergeMode()
        if isempty(FV.vXlim), axis(hSubplots(end), 'tight')
        else set(hSubplots(end), 'xlim', FV.vXlim); end
    else
        nX_s = vTime(1); nX_e = vTime(end);
        set(hSubplots(end), 'xlim', [nX_s nX_e]);
    end
    
    % Plot spike threshold line(s)
    hCntxtMenu  = uicontextmenu;
    if isfield(FV.tSpikeThresholdsPos, sCh) && ~FV.bPlotRasters % positive
        nThresh = FV.tSpikeThresholdsPos.(sCh);
        hold on;
        hThresh = plot([vTime(1) vTime(end)], [nThresh nThresh], ':'); hold off
        set(hThresh, 'UIContextMenu', hCntxtMenu, 'Tag', [sCh '_pos']);
        ThemeObject(hThresh)
        if isfield(FV.tYlim, sCh)
            if isempty(FV.tYlim.(sCh))
                FV.tYlim.(sCh) = [-nThresh nThresh].*1.1;
            else
                FV.tYlim.(sCh)(2) = max([FV.tYlim.(sCh)(1) nThresh*1.1]);
            end
        end
    end
    if isfield(FV.tSpikeThresholdsNeg, sCh) && ~FV.bPlotRasters % negative
        nThresh = FV.tSpikeThresholdsNeg.(sCh);
        hold on;
        hThresh = plot([vTime(1) vTime(end)], [nThresh nThresh], ':'); hold off
        set(hThresh, 'UIContextMenu', hCntxtMenu, 'Tag', [sCh '_neg']);
        ThemeObject(hThresh)
        if isfield(FV.tYlim, sCh)
            if isempty(FV.tYlim.(sCh))
                FV.tYlim.(sCh) = [nThresh -nThresh].*1.1;
            else
                FV.tYlim.(sCh)(1) = min([FV.tYlim.(sCh)(1) nThresh*1.1]);
            end
        end
    end

    % Plot threshold line for event detection (if set for channel)
    if isfield(FV, 'tEventThresholds')
        nInd = find(strcmp(sCh, {FV.tEventThresholds.sChannel}));
        if ~isempty(nInd)
            nThresh = FV.tEventThresholds(nInd).nThreshold;
            hold on;
            hThresh = plot([vTime(1) vTime(end)], [nThresh nThresh], ':'); hold off
            ThemeObject(hThresh)
            set(hThresh, 'color', FV.mColors(1,:))
        end
    end
    
    uimenu(hCntxtMenu, 'Label', 'Delete', 'Callback', @RemoveSpikeThreshold)
    uimenu(hCntxtMenu, 'Label', 'Copy', 'Callback', @CopyPasteThreshold)
    uimenu(hCntxtMenu, 'Label', 'Paste', 'Callback', @CopyPasteThreshold)
    uimenu(hCntxtMenu, 'Label', 'Set Threshold', 'Callback', @SetThresholdManually, 'separator', 'on')
    
    % Plot spikes on top of continuous trace (only if rasters are not displayed)
    if isfield(FV.tSpikes, sCh) && ~FV.bPlotRasters
        vSpiketimes = FV.tSpikes.(sCh).spiketimes; % nFs resolution
        if ~isempty(vSpiketimes)
            vSpiketimes = vSpiketimes ./ nFs;
            vSpiketimes(vSpiketimes < vTime(1) | vSpiketimes > vTime(end)) = [];
            vVolts = repmat(nThresh, length(vSpiketimes), 1);
            hold on; hDot = plot(hSubplots(end), vSpiketimes', vVolts', '.');
            ThemeObject(hDot);
            set(hDot, 'ButtonDownFcn', @IdentifySpike)
        end
    end
    
    % Set Y axis limits
    if ~isfield(FV, 'bPlotInstSpikeRate'), FV.bPlotInstSpikeRate = 0; end
    if FV.bPlotInstSpikeRate && isfield(FV.tSpikes, sCh)
        % do nothing for now
    elseif FV.bPlotRasters && isfield(FV.tSpikes, sCh)
        % Indicate unit number and quality in y-labels
        if isfield(FV.tSpikes.(sCh), 'quality')
            tQuality = FV.tSpikes.(sCh).quality;
            for nu = 1:length(csUnits)
                nIndx = [tQuality.unit] == str2num(csUnits{nu});
                if isempty(tQuality(nIndx))
                    csUnits{nu} = [csUnits{nu} ' (Q?)'];
                else
                    csUnits{nu} = [csUnits{nu} ' (Q' num2str(tQuality(nIndx).score) ')'];
                end
            end
        end
        set(hSubplots(end), 'ylim', [.5 nRow+.5], 'ytick', 1:1:nRow, 'yticklabel', csUnits);
    elseif isfield(FV.tYlim, sCh)
        if size(vCont, 1) == 1
            if ~isempty(FV.tYlim.(sCh)) && ~any(isnan(FV.tYlim.(sCh)))
                nYd = (diff(FV.tYlim.(sCh)) * .1) .* [-1 1];
                vYLim = FV.tYlim.(sCh) + nYd;
                if vYLim(1) == vYLim(2)
                    nD = vYLim(1) * .05 - (10^-6);
                    set(hSubplots(end), 'ylim',  vYLim - [-nD nD]); % use ylim +/- 5%
                else
                    set(hSubplots(end), 'ylim',  vYLim); % use saved ylim * 5%
                end
                % Reduce number of y-ticks
                vYTicks = get(hSubplots(end), 'ytick');
                set(hSubplots(end), 'ytick', vYTicks(1:2:end));
            end
        else
            % 2D plot
            if exist('vY', 'var')
                set(hSubplots(end), 'ylim', [vY(1) vY(end)])
            end
        end
    end
end

% Plot digital events
if bShowDigitalEvents
    % Create subplot
    hSubplots(end+1) = axes('position', [.07 .985-nSubEventHeight .91 nSubEventHeight]);
    hold on
    
    % Get number of events
    nEvents = length(FV.csDigitalChannels); % number of events
    if any(strcmp('DAQ_Start', FV.csDigitalChannels))
        nEvents = nEvents - 1;
    end
    if any(strcmp('DAQ_Stop', FV.csDigitalChannels))
        nEvents = nEvents - 1;
    end
    if any(strcmp('DAQ_Trigger', FV.csDigitalChannels))
        nEvents = nEvents - 1;
    end
    
    % Get height of event markers
    set(hSubplots(end), 'units', 'pixels')
    vPos = get(hSubplots(end), 'position');
    nEventHeight = vPos(4) / nEvents / 2;
    set(hSubplots(end), 'units', 'normalized')
    
    % Iterate and plot digital events
    nY = 0;
    for nF = 1:length(FV.csDigitalChannels) % iterate over fieldnames
        % Ignore DAQ_Start, DAQ_Stop and DAQ_Trigger (default fields generated by DA Toolbox)
        if any(strcmp(FV.csDigitalChannels(nF), {'DAQ_Start', 'DAQ_Stop', 'DAQ_Trigger'}))
            continue
        end
        
        [vUpTimes, vDownTimes] = GetEventPairs(FV.csDigitalChannels{nF});

        % Do not display channels without events
        if isempty(vDownTimes) || isempty(vUpTimes)
            iIgnoreEvents(nF) = 1;
            continue
        end
        iIgnoreEvents(nF) = 0;
        nY = nY + 1;

        % Make sure events are arranged as a row index
        vUpTimes = vUpTimes(:)';
        vDownTimes = vDownTimes(:)';
        
        % NEW:
        % In previous versions, individual events were plotted with individual lines, like this:
        %  plot([vUpTimes;vDownTimes], repmat([nF;nF], 1, length(vUpTimes)), 'linewidth', 5, 'color', [1 1 0])
        % In the new code, 'gaps' will be produced with NaNs, drastically speeding up display updates
        % since only a single line now needs to be plotted for each event
        %
        % X values:
        % [1  2  2  3    3  4  4  5]
        % [U1 D1 D1 U2   U2 D2 D2 U3]
        %
        % Y values:
        % [1  1  0  0    1  1  0  0] where 0 is NaN
        %
        vIndices = [vUpTimes' vDownTimes' vDownTimes' [vUpTimes(2:1:end) NaN]'];
        vIndices = reshape(vIndices', numel(vIndices), 1);
        % Rearrange indices
        vYs = repmat([nY nY NaN NaN], 1, length(vIndices)/4);
        
        %if length(vIndices) > 100000000
        %    vPlotIndx = 1:round(length(vIndices)/1000):length(vIndices);
        %else
            vPlotIndx = 1:length(vIndices);
        %end
        
        % Set event width
        %vXlim = FV.vXlim;
        %if isempty(FV.vXlim)
        %    vXlim = [min([vUpTimes vDownTimes]) max([vUpTimes vDownTimes])];
        %end
        %nEventWidth = GetEventWidth(vXlim)/3;
        %iEnd = 2:2:length(vIndices);
        %iRep = [vIndices(iEnd) - vIndices(1:2:end)] < nEventWidth;
        %vIndices(iRep) = vIndices(iEnd(iRep)-1) + nEventWidth;
        
        plot(hSubplots(end), vIndices(vPlotIndx), vYs(vPlotIndx), 'linewidth', nEventHeight, 'color', FV.mColors(1,:))
    end

    % Plot file start (green) and end (red) markers
    if isfield(FV.tData, 'FileStart') % only in Merge Mode
        for i = 1:length(FV.tData.FileStart)
            nX = FV.tData.FileStart(i).Timestamp;
            if ~isempty(FV.vXlim)
                nRAdj = abs(diff(FV.vXlim)*.5); % adjust extended range by 50% for panning purposes
                if nX < (min(FV.vXlim)-nRAdj) || nX > (max(FV.vXlim)+nRAdj), continue, end
            end
            vY = [0 nF+.5];
            hLin = plot([nX nX],  vY, 'color', 'g');
            hFMMenu = uicontextmenu;
            uimenu(hFMMenu, 'Label', 'Open File', 'Callback', @FileMarkerOpenFile, 'UserData', FV.tData.FileStart(i).File);
            set(hLin, 'DisplayName', FV.tData.FileStart(i).File, 'ButtonDownFcn', @FileMarkerShowPath, ...
                'uicontextmenu', hFMMenu )
        end
    end

    if isfield(FV.tData, 'FileEnd') % only in Merge Mode
        for i = 1:length(FV.tData.FileEnd)
            nX = FV.tData.FileEnd(i).Timestamp;
            if ~isempty(FV.vXlim)
                nRAdj = abs(diff(FV.vXlim)*.5); % adjust extended range by 50% for panning purposes
                if nX < (min(FV.vXlim)-nRAdj) || nX > (max(FV.vXlim)+nRAdj), continue, end
            end
            vY = [0 nF+.5];
            hLin = plot([nX nX],  vY, 'color', 'r');
            hFMMenu = uicontextmenu;
            uimenu(hFMMenu, 'Label', 'Open File', 'Callback', @FileMarkerOpenFile, 'UserData', FV.tData.FileEnd(i).File);
            set(hLin, 'DisplayName', FV.tData.FileEnd(i).File, 'ButtonDownFcn', @FileMarkerShowPath, ...
                'uicontextmenu', hFMMenu )
        end
    end

    % Set axis properties
    ThemeObject(hSubplots(end));
    set(hSubplots(end), 'Tag', 'Events')
    
    if IsMergeMode()
        if isempty(FV.vXlim), axis tight
        else set(hSubplots(end), 'xlim', FV.vXlim); end
    else
        if isempty(FV.vXlim)
            if exist('vTime', 'var')
                nX_s = vTime(1);
                nX_e = vTime(end);
                set(hSubplots(end), 'xlim', [nX_s nX_e]);
            end
        else set(hSubplots(end), 'xlim', FV.vXlim); end
    end

    % Get y-tick label strings
    cYTickLabels = cell(1, length(FV.csDigitalChannels));
    for c = 1:length(FV.csDigitalChannels)
        % Ignore DAQ_Start, DAQ_Stop and DAQ_Trigger (default fields generated by DA Toolbox)
        if any(strcmp(FV.csDigitalChannels(c), {'DAQ_Start', 'DAQ_Stop', 'DAQ_Trigger'}))
            cYTickLabels{c} = '';
            continue
        end
        sDescr = GetChannelDescription(FV.csDigitalChannels{c});
        if isempty(sDescr)
            sDescr = FV.csDigitalChannels{c};
        end
        cYTickLabels{c} = sDescr;
    end

    if isempty(cYTickLabels)
        set(hSubplots(end), 'ylim', [0 1], 'ytick', [], 'yticklabel', [])
        hTxt = text(mean(get(hSubplots(1), 'xlim')), .5, 'No digital events found', 'horizontalalignment', 'center');
        ThemeObject(hTxt);
    else
        if all(ishandle(hSubplots))
            % Remove underscore in Y labels
            cYTickLabels = strrep(cYTickLabels, '_', '');
            set(hSubplots(end), 'ylim', [0.65 max([1, sum(~iIgnoreEvents)])+.35], ...
                'ytick', 1:sum(~iIgnoreEvents), 'yticklabel', cYTickLabels(~iIgnoreEvents))
        end
    end
    box on
    
    % Context menu in Events axis
    hMenu = uicontextmenu;
    uimenu(hMenu, 'Label', 'Event Statistics', 'Callback', @ShowEventStatistics);
    uimenu(hMenu, 'Label', 'Undigitize Channel', 'Callback', @UndigitizeChannel);
    uimenu(hMenu, 'Label', 'Edit Event Durations (B)', 'Callback', @ModifyEventDuration);
    uimenu(hMenu, 'Label', 'Edit Event Intervals (B)', 'Callback', @ModifyEventInterval);
    uimenu(hMenu, 'Label', 'Create Virtual Event', 'Callback', @CreateVirtualEvent);
    
    set(hSubplots(end), 'uicontextmenu', hMenu)
end

% Show x-axis grid
if ~isempty(hSubplots)
    if strcmpi(get(findobj(hFig, 'Tag', 'Spiky_Menu_ShowGrid'), 'Checked'), 'on')
        set(hSubplots(ishandle(hSubplots)), 'xgrid', 'on');
    end
end

% Get x-limits of all subplots and adjust so all are same
cXlim = get(hSubplots(ishandle(hSubplots)), 'xlim');

if size(cXlim,1) > 1
    mXlim = reshape([cXlim{:}],2,length([cXlim{:}])/2)';
else
    mXlim = cXlim;
end

% Remove [0 1] x-limits (likely an empty Events axes)
nRmIndx = find(ismember(mXlim, [0 1], 'rows'));
mXlim(nRmIndx, :) = [];

if ~isempty(mXlim)
    nXlim_MinMax = [min(mXlim(:,1)) max(mXlim(:,2))];
    set(hSubplots(ishandle(hSubplots)), 'xlim', nXlim_MinMax);
end

% Display normal time (i.e. zero is start of trial)
if strcmp(get(findobj(hFig, 'Tag', 'ShowNormalTime'), 'checked'), 'on') && ~isempty(hSubplots)
    for hSub = hSubplots(:)'
        if ~ishandle(hSub), continue, end
        vXlim = get(hSub, 'xlim');
        vXtick = get(hSub, 'XTick');
        if ~exist('nBeginTime', 'var'), nBeginTime = 0; end
        vXtickNormal = vXtick - nBeginTime;
        if diff(vXlim) > 10 % Round to sec if range > 10 sec
            nNormFact = 1;
        elseif diff(vXlim) > 1  % Round to 0.1 sec if range > 1 sec
            nNormFact = 10;
        elseif diff(vXlim) > .1  % Round to 0.01 sec if range > 0.1 sec
            nNormFact = 100;
        else  % Round to 0.001 sec if range < 0.1 sec
            nNormFact = 1000;
        end
        vXtickRound = round(vXtickNormal.*nNormFact)./nNormFact;
        set(hSubplots, 'XTicklabel', vXtickRound)
    end
end

% Remove xtick marks on all but bottom subplot
hSubplots(~ishandle(hSubplots)) = [];
set(hSubplots(2:end), 'XTicklabel', [])

% Re-arrange order of figure children such that top subplot is lowest
hHandles = get(hFig, 'children');
iAxes = [];
for h = 1:length(hHandles)
    if strcmp(get(hHandles(h), 'type'), 'axes')
        iAxes(end+1) = h;
    end
end

mPos = cell2mat(get(hHandles(iAxes), 'position')); % axes positions
[~, iOrder] = sort(mPos(:, 2)); % order by Y position
hHandles(iAxes) = hHandles(iAxes(iOrder));
set(hFig, 'children', hHandles)

SetStruct(FV, 'nosaveflag')
PanMode()
set(GetGUIHandle(), 'WindowButtonMotionFcn', @GUIMouseMotion);
return


function nBeginTime = GetBeginTime()
% Set the begin time of the experiment (absolute time, sec from midnight)
% 
% Usage:
%   GetBeginTime()
%
[FV, ~] = GetStruct();
csFieldnames = fieldnames(FV.tData);
cMatch = regexpi(csFieldnames, '^*_TimeBegin$');
vBeginTime = [];
for i = 1:length(cMatch)
    if ~isempty(cMatch{i})
        vBeginTime(end+1) = FV.tData.(csFieldnames{i});
    end
end
nBeginTime = mode(vBeginTime);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetTimeMarker(nT)
% Display a vertical time marker in all axis
%
% This function can be useful for external scripts that need to display, in
% the GUI, the time at which a signal is analysed.
%
% Usage:
%   SetTimeMarker(T) where T is the absolute time in seconds.
%

% Get axes handles
hGUI = GetGUIHandle();
hAx = findobj(hGUI, 'type', 'axes');

% If there are fewer markers than axes, then remove all
hMarkers = findobj(hAx, 'tag', 'SpikyTimeMarker');
if length(hMarkers) < length(hAx)
    delete(hMarkers)
end

% Create markers if they don't exist
if isempty(hMarkers)
    % Iterate over axes
    for ax = 1:length(hAx)
        hold(hAx(ax), 'on')
        hLin = plot(hAx(ax), [nT nT], get(hAx(ax), 'ylim'));
        set(hLin, 'YDataSource', 'get(get(gcbo,''parent''), ''ylim'')', ...
            'color', 'r', 'linewidth', 2, 'linestyle', '-', ...
            'tag', 'SpikyTimeMarker')
    end
else
    % Otherwise update the markers
    set(hMarkers, 'xdata', [nT nT])
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vUpTimes, vDownTimes] = GetEventPairs(csCh, varargin)
% Get event pairs (up and down times) of a digital channel
% 
% Usage:
%   GetEventPairs(S) where S is a string and the channel name
%   GetEventPairs(S, 'range/all') retrieves events in the currently
%   displayed range ('range' default) or all events ('all')
%
[FV, ~] = GetStruct();

vUpTimes = FV.tData.([csCh '_Up']); % UP times, sec
vUpTimes = sort(vUpTimes); % Sort

vDownTimes = FV.tData.([csCh '_Down']); % DOWN times, sec
vDownTimes = sort(vDownTimes); % Sort
if isempty(vDownTimes), vDownTimes = FV.tData.([csCh '_TimeEnd']); end

% Remove DOWN times that occur before UP events
if ~isempty(vUpTimes)
    vDownTimes(vDownTimes <= vUpTimes(1)) = [];
end

% If up and down times differ in length, crop the longest
if length(vUpTimes) > length(vDownTimes)
    vUpTimes = vUpTimes(1:length(vDownTimes));
elseif length(vUpTimes) < length(vDownTimes)
    vDownTimes = vDownTimes(1:length(vUpTimes));
end

% Remove events outside of axis limits
bGetEventsInRange = 1;
if nargin > 1
    if strcmp(varargin{1}, 'all')
        bGetEventsInRange = 0;
    end
end
if bGetEventsInRange && (~isempty(FV.vXlim) && ~isempty(vUpTimes) && ~isempty(vDownTimes))
    vIndx = find((vUpTimes < [min(FV.vXlim)] | vUpTimes > [max(FV.vXlim)]) ...
        & (vDownTimes < [min(FV.vXlim)] | vDownTimes > [max(FV.vXlim)]));
    vUpTimes(vIndx) = [];
    vDownTimes(vIndx) = [];
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get display width of events and spikes (based on time displayed interval)
function nW = GetEventWidth(vXlim)
if diff(vXlim) < 1
    nW = 0.001;
elseif diff(vXlim) < 5
    nW = 0.005;
elseif diff(vXlim) < 10
    nW = 0.01;
elseif diff(vXlim) < 30
    nW = 0.02;
elseif diff(vXlim) < 60
    nW = 0.04;
elseif diff(vXlim) < 150
    nW = 0.06;
elseif diff(vXlim) < 250
    nW = 0.12;
elseif diff(vXlim) < 500
    nW = 0.5;
else
    nW = 1;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change the display order of axes in main window
function ChangeDisplayOrder(hObject, eventdata, sCh, nMove)
[FV, ~] = GetStruct();
iCurr = find(strcmp(FV.csDisplayChannels, sCh));
iMove = iCurr + nMove;
iMove = max([iMove 1]);
iMove = min([iMove length(FV.csDisplayChannels)]);
csD = FV.csDisplayChannels;
csMove = csD(iCurr);
csD(iCurr) = [];
csD = [csD(1:(iMove-1)) csMove csD((iMove):end)];
FV.csDisplayChannels = csD;
SetStruct(FV, 'nosaveflag')
ViewTrialData();
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetTheme(hObject, varargin)
% Select theme to use
[FV, ~] = GetStruct();
switch get(hObject, 'selected')
    case 'on'
        set(hObject, 'Checked', 'on')
    case 'off'
        set(hObject, 'Checked', 'off')
end
set(setdiff(findobj('tag', 'Spiky_Menu_ShowTheme'), hObject), 'Checked', 'off')
FV.CurrentTheme = get(hObject, 'Label');
SetStruct(FV)
ViewTrialData();
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ThemeObject(hObj, varargin)
% Update properties of object based on currently selected theme.
% 
% Usage:
%   ThemeObject(H)
%   ThemeObject(H, 'PropertyName', 'PropertyValue')
% 
%   where H is a graphics object handle
% 
%
[FV, ~] = GetStruct();

% Load properties only when theme changes
persistent p
if ~isfield(FV, 'CurrentTheme')
    FV.CurrentTheme = 'spiky';
end
if isempty(p) || ~strcmp(FV.CurrentTheme, p.theme)
    % Load selected theme properties
    sDir = which(mfilename);
    sDir = [sDir(1:end-7) 'themes' filesep];
    tDir = dir(sDir);
    tDir = tDir(3:end);
    if ~isfield(FV, 'CurrentTheme'); FV.CurrentTheme = 'spiky'; end
    nIndx = strcmp([FV.CurrentTheme '.m'], {tDir.name});
    run([sDir tDir(nIndx).name])
    p.theme = FV.CurrentTheme;
    % Update line colors
    if isfield(p, 'axes')
        FV.mColors = p.colors;
    end
    SetStruct(FV)
end

% Set object properties
for nh = 1:length(hObj)
    if ~isprop(hObj(nh), 'type') && ~isprop(hObj(nh), 'type')
        return;
    end
    sType = get(hObj(nh), 'type');
    if isfield(p, sType) && ishandle(hObj(nh))
        set(hObj(nh), p.(sType))
    end
    % If Type is 'axes', then set backdrop
    if strcmp(sType, 'axes') && 0
        if isfield(p, 'backdrop')
            hAxes = handle(hObj(nh));
            try
                hAxes.Backdrop.Face.reset
                hAxes.Backdrop.Face.ColorData = p.backdrop.face.ColorData;
                hAxes.Backdrop.Face.ColorBinding = 'interpolated';
            catch
            end
        end
    end
end

% Set additional custom properties passed in varargin
if nargin >= 3
    for pp = 1:2:length(varargin)
        for nh = 1:length(hObj)
            set(hObj(nh), varargin{pp}, varargin{pp+1})
        end
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CopyChannelToFig(varargin)
% Copy right-clicked channel in main window to a separate, default Figure.
[FV, hWin] = GetStruct();
sCh = get(gcbo, 'Tag');
if ~isfield(FV.tData, sCh); return; end
vCont = ChannelCalculator(FV.tData.(sCh), sCh);
[vCont, FV] = AdjustChannelGain(FV, vCont, sCh);
nBegin = FV.tData.([sCh '_TimeBegin']); % sampling start, sec
nEnd = FV.tData.([sCh '_TimeEnd']); % sampling end, sec
nFs = FV.tData.([sCh '_KHz']) * 1000; % sampling frequency Hz
vTime = (nBegin+1/nFs):(1/nFs):(nBegin+length(vCont)/nFs); % absolute time, sec
if strcmp(get(findobj(hWin, 'Tag', 'ShowNormalTime'), 'checked'), 'on')
    vTime = vTime - nBegin;
end
hFig = figure;
ThemeObject(hFig)
hAx = axes();
if all(size(vCont) > 1) % 2D
    imagesc(vCont)
    if isfield(FV.tData, ([sCh '_Unit']))
        vY = FV.tData.([sCh '_Scale']);
    else
        vY = 1:size(vCont, 2);
    end
    hLin = imagesc(vTime, vY, real(vCont));
    set(gca, 'ydir', 'normal')
    xlabel('Time (s)')
    ylabel(FV.tData.([sCh '_Unit']))
else
    hLin = plot(vTime, vCont, 'k', 'linewidth', 1);
    xlabel('Time (s)')
    ylabel(FV.tAmplitudeUnit.sUnit)
end
ThemeObject(hLin);
ThemeObject(hAx)
hTit = title(['Channel ' sCh], 'interpreter', 'none');
ThemeObject(hTit)
axis tight
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FileMarkerShowPath(varargin)
hLeg = legend(gcbo);
hTxt = findobj(hLeg, 'type', 'text');
ThemeObject(hTxt)
set(hTxt,'interpreter','none')
legend boxoff;
cTxt = get(hLeg,'string');
if ~isempty(cTxt); sp_disp(cTxt{1}); end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FileMarkerOpenFile(varargin)
sPath = get(gcbo, 'userdata');
if exist(sPath, 'file'); OpenFile([], sPath); end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vSpiketimes = DropDeadtimeSpikes(vSpiketimes)
[FV, ~] = GetStruct();
% round off spiketimes to deadtime resolution and keep first of all unique values
vSpiketimesRound = round(vSpiketimes ./ (FV.nDeadTimeMs/1000)) .* (FV.nDeadTimeMs/1000); % sec
[~, vRoundIndx] = unique(vSpiketimesRound, 'first');
vSpiketimes = vSpiketimes(vRoundIndx);
vISI = diff(vSpiketimes); % sec
vDropIndx = find(vISI <= (FV.nDeadTimeMs / 1000)); % ms
vSpiketimes(vDropIndx+1) = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PanWindow(varargin)
[FV, hWin] = GetStruct();
global g_hSpike_Viewer % handle to Spiky window
if ~isfield(FV, 'tData')
    set(gcbo, 'state', 'off')
    return
end
% update menu item
hMenuItem = findobj(g_hSpike_Viewer, 'Label', 'Pan'); % handle of 'Plot Rasters'
switch get(hMenuItem, 'Checked')
    case 'on' % turn off pan mode
        set(hMenuItem, 'Checked', 'off');
        FV.bPanOn = 0;
    case 'off' % turn on pan mode
        set(hMenuItem, 'Checked', 'on');
        FV.bPanOn = 1;
end
SetStruct(FV)
PanMode()
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PanMode(varargin)
[FV, hWin] = GetStruct();
hPan = pan(hWin);
if FV.bPanOn
    set(hPan, 'enable', 'on', 'Motion', 'horizontal', 'ActionPostCallback', @UpdatePannedWindow)
else
    try
        set(hPan, 'enable', 'off')
    catch
        [~, sMsgID] = lastwarn;
        if ~strcmp(sMsgID, 'MATLAB:modes:mode:InvalidPropertySet')
            set(hWin, 'WindowButtonDownFcn', '')
        end
    end
end
linkaxes(findobj(get(GetGUIHandle(), 'children'), 'type', 'axes'), 'x')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdatePannedWindow(varargin)
% Wait until x-axis limits have changed
vOrigLims = get(gca, 'xlim');
% Get new x limits
[FV, ~] = GetStruct();
vCurrLim = get(gca, 'xlim');
FV.vXlim = vCurrLim; % pan left by 25%
SetStruct(FV); ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vCont, vTime, nNewFs] = FilterChannel(vCont, vTime, nFs, nLoPass, nHiPass, bRectify, sOption)
% Bandpass vector.
%
% Usage:
%   FilterChannel(C, T, Fs, LO, HI, R, OPT)
%       where   C   signal to filter
%               T   time vector (only required with 'decimate' option)
%               Fs  sampling rate (Hz)
%               LO  low-pass frequency (Hz)
%               HI  high-pass frequency (Hz)
%               R   1 = rectify, 0 = no rectify
%               OPT 'decimate' or 'none'
%
% See:
%   GetFilteredChannel()
%

% Check for sane inputs
nLoPass = min([nLoPass nFs]);
if isnan(nLoPass) || isnan(nHiPass) || (nHiPass >= nLoPass)
    nNewFs = nFs;
    return
end

% Interpolate NaN indices
vNaNIndx = isnan(vCont);
if ~all(vNaNIndx)
    vCont(vNaNIndx) = interp1(find(~vNaNIndx), vCont(~vNaNIndx), find(vNaNIndx), 'linear', 'extrap');
end

% High pass filter
if nHiPass > 0
    [b,a] = butter(3, nHiPass/(double(nFs)/2) , 'high'); % high-pass
    vCont  = single(filtfilt(b, a, double(vCont)));
end

% Low-pass filter
nLPPass = nLoPass/(double(nFs)/2);
if nLPPass < 1 && nLPPass > 0
    [b,a] = butter(3, nLoPass/(double(nFs)/2) , 'low'); % low-pass
    if length(vCont) > (3 * length(a))
        vCont  = single(filtfilt(b, a, double(vCont)));
    end
end

% Decimate signal (optional)
switch lower(sOption)
    case 'decimate'
        % decimate (i.e. resample signal)
        nNewFs = nLoPass * 2.5;
        nStep = ceil(nFs / nNewFs);
        nNewFs = nFs/nStep; % new Fs must be an integer
        if ~isempty(vTime)
            vTimeO = vTime;
            vTime = vTime(1:nStep:end);
            vNaNIndx = unique(round(interp1(vTime, 1:length(vTime), vTimeO(vNaNIndx))));
            vNaNIndx(find(vNaNIndx) > length(vTime)) = [];
        end
        vCont = vCont(1:nStep:end);
    otherwise
        nNewFs = nFs;
end

% Rectify
if bRectify, vCont = abs(vCont); end

% Return NaN's where they were removed above
if any(vNaNIndx)
    vCont(vNaNIndx(~isnan(vNaNIndx))) = NaN;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function csSelected = SelectChannels(varargin)
% Select the channels that should be displayed in the GUI
% 
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();

% If a channel was specific in varargin{1}, then toggle this in csDisplayChannels,
% i.e. remove if its there or insert if its not already displayed
if ~isempty(varargin)
    if ishandle(varargin{1})
        sCh = get(varargin{1}, 'Tag');
        % Check that tag of handle indeed references a known channel
        if ~(isempty(sCh) || ~any(strcmp(fieldnames(FV.tData), sCh)))
            nIndx = find(strcmp(FV.csDisplayChannels, sCh));
            if isempty(nIndx) % is not already displayed
                FV.csDisplayChannels{end+1} = sCh; % show
            else % is already displayed
                FV.csDisplayChannels(nIndx) = []; % remove
            end
            SetStruct(FV)
            ViewTrialData();
            return
        end
    end
end

global g_hSpike_Viewer
% create list of selectable channels
cFields = fieldnames(FV.tData);
csDisplayChannels = {};
csChannelDescriptions = {};
for i = 1:length(cFields)
    if ~isempty(strfind(cFields{i}, '_TimeBegin'))
        if isfield(FV.tData, cFields{i}(1:end-10))
            csDisplayChannels{end+1} = cFields{i}(1:end-10);
            if isfield(FV, 'tChannelDescriptions')
                if ~isempty(FV.tChannelDescriptions)
                    nIndx = find(strcmpi({FV.tChannelDescriptions.sChannel}, csDisplayChannels{end}));
                else nIndx = []; end
            else nIndx = []; end
            if isempty(nIndx)
                csChannelDescriptions{end+1} = '';
            else
                csChannelDescriptions{end+1} = FV.tChannelDescriptions(nIndx).sDescription;
            end
        end
    end
end

% Sort alphabetically
[csDisplayChannels, iOrder] = sort(csDisplayChannels);
csChannelDescriptions = csChannelDescriptions(iOrder);

global g_vSelCh;
if length(csDisplayChannels) == 1
    g_vSelCh = 1;
else
    if isempty(csDisplayChannels), return, end
    hFig = figure;
    nFigWidth = 300;
    nCL = length(csDisplayChannels)*25+35;
    vHW = get(g_hSpike_Viewer, 'position');
    set(hFig, 'position', [vHW(1)+40 (vHW(2)+vHW(4)-nCL-50) nFigWidth nCL+25], 'menu', 'none', 'Name', 'Select channels', 'NumberTitle', 'off')
    for i = 1:length(csDisplayChannels)
        sBoxStr = sprintf('%s (%s)', csChannelDescriptions{i}, csDisplayChannels{i});
        if any(strcmp(FV.csDisplayChannels, csDisplayChannels{i}))
            uicontrol(hFig, 'Style', 'checkbox', 'Position', [10 nCL nFigWidth-20 20], 'String', sBoxStr, 'HorizontalAlignment', 'left', 'backgroundcolor', get(hFig, 'color'), 'value', 1);
        else
            uicontrol(hFig, 'Style', 'checkbox', 'Position', [10 nCL nFigWidth-20 20], 'String', sBoxStr, 'HorizontalAlignment', 'left', 'backgroundcolor', get(hFig, 'color'), 'value', 0);
        end
        nCL = nCL - 25;
    end
    % 'Select all' checkbox
    uicontrol(hFig, 'Style', 'checkbox', 'Position', [10 nCL 170 20], 'String', 'Select all', ...
        'HorizontalAlignment', 'left', 'backgroundcolor', get(hFig, 'color'), 'value', 0, 'callback', ...
        'set(get(gcf,''children''),''value'', 1)');
    nCL = nCL - 25;
    uicontrol(hFig, 'Style', 'pushbutton', 'Position', [50 nCL 100 20], ...
        'Callback', 'global g_vSelCh; g_vSelCh=flipud(get(get(gcf,''children''),''value'')); g_vSelCh=[g_vSelCh{:}];close(gcf)', ...
        'String', '      OK      ' ); % OK button
    uiwait % wait for user to close window or click OK button
    g_vSelCh = logical(g_vSelCh(1:end-2));
end
% Assign channels to FV or just return them
bReturn = 0;
csSelected = csDisplayChannels(g_vSelCh);
csSelectedDescr = csChannelDescriptions(g_vSelCh);

if ~isempty(varargin)
    if strcmp(varargin{1}, 'return') || isempty(g_vSelCh)
        bReturn = 1;
    end
end
% Return selected channels back to calling function
if bReturn
    return
else % Assign channels to FV and update display
    
    % Sort channels in alphabetical order based on their Description
    % strings
    iEmpty = strcmp(csSelectedDescr, '');
    csSelectedDescr(iEmpty) = csSelected(iEmpty);
    [~, iSorted] = sort(csSelectedDescr);
    
    FV.csDisplayChannels = csSelected(iSorted); % selected channels
    clear global g_vSelCh
    SetStruct(FV)
    ViewTrialData
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExportData(varargin)
% Export data to a different file format
% 
% Usage:
%   ExportData()  invokes a dialog for selecting output file and format
%   ExportData(F) where F is the filename and filetype in a supported format
%                 (see /export directory). Example:
%                   ExportData('/mydir/myfile.mat')
%
if ~IsDataLoaded(), return; end

[FV, ~] = GetStruct();

sCurrFile = FV.sLoadedTrial;
sCurrFile(strfind(sCurrFile, '.'):end) = [];

% Validate input path
if ~isempty(varargin)
    if ~isstr(varargin{1})
        varargin = [];
    end
end

if isempty(varargin)
    [sFile, sPath] = uiputfile(GetExportFilters(), 'Select file to write', [FV.sDirectory filesep sCurrFile]);
    if sFile == 0; return; end
else
    [sPath, sFile] = fileparts(varargin{1});
end

% Run export filter
eval(sprintf('export_%s([sPath sFile], FV);', lower(sFile(end-2:end))))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bStop = StopAcceleratorRepeat
% StopAcceleratorRepeat() prevents the same function being repeated more
% than once in rapid succession when an accelerator key combination is
% pressed. In what appears to be a bug, Matlab occasionally registers a
% continuously held down accelerator combo as multiple key-presses, rather
% than a single press. When multiple accelerator key callbacks are
% executed, add this line to the beginning of your function:
%
%   if StopAcceleratorRepeat(), return; end
%
bStop = 0;
persistent tRepeats
nThresh = 0.75; % threshold, seconds

% Get the calling function
csDBStack = dbstack();
sCallingFunction = csDBStack(2).name;

% Confirm that spiky is the first callback function. Then, it is at least possibly
% we got here via an accelerator key. There seems to be no way to know for sure.
if length(csDBStack) < 3, return; end
if ~strcmpi(csDBStack(3).name, 'spiky'), return; end

if isfield(tRepeats, sCallingFunction)
    vLast = tRepeats(1).(sCallingFunction);
    if etime(clock, vLast) <= nThresh, bStop = 1; end
end
tRepeats(1).(sCallingFunction) = clock;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OpenFile(varargin)
% OpenFile loads both trial data and the associated settings/results file. This
% function is called from all menu items in the GUI and is a higher-level function
% than LoadTrial. Generally, OpenFile is the only function that calls LoadTrial.
% 
% OpenFile() is called from the GUI (File->Open or the Open button in the
% toolbar)
% 
% TODO Ignore first input when empty
% 
% Usage:
%   OpenFile([], -1)    Load previous file in current list
%   OpenFile([], 0)     Load next file in current list
%   OpenFile([], N)     Load N'th file in current list
%   OpenFile([], S)     Load the file with name S
%   OpenFile()          Select a file interactively
%
if StopAcceleratorRepeat(), return; end
[FV, hWin] = GetStruct();

global g_hSpike_Viewer g_bBatchMode;
persistent p_bNeverApply;

% Remove inputs if function was called from GUI
if ishandle(varargin{1})
    varargin = {[] []};
end

% Decide which file to load next
if isempty(varargin{2})
    % Interactively choose a file to open
    if isfield(FV, 'sDirectory')

        % TODO Use extension of currently loaded file as default
        if isfield(FV, 'sLoadedTrial')
            [~, ~, sDefExt] = fileparts(FV.sLoadedTrial);
            csExt = GetImportFilters(['*' sDefExt]);
        else
            csExt = GetImportFilters(FV.sDirectory);
        end
        [sFile, sPath] = uigetfile(csExt, 'Select data file', [FV.sDirectory filesep]);
    else
        [sFile, sPath] = uigetfile(GetImportFilters, 'Select data file');
        set(g_hSpike_Viewer, 'UserData', FV);
    end
    g_bBatchMode = false;
    if sFile == 0, return, end
    FV.sDirectory = sPath; % update current path
    
    % Load merge file
    if ~isempty(strfind(sFile, 'MergeFile')) && ~isempty(strfind(sFile, '.spb'))
        LoadMergeFile(sPath, sFile)
        return
    end
    
    % Get the index of this .mat file in its directory
    SetStruct(FV, 'nosaveflag'); % update current path
    sNextTrial = sFile;

    % Clear current file list
    hFileList = findobj(g_hSpike_Viewer, 'Tag', 'MenuFileList');
    delete(findobj(g_hSpike_Viewer, 'Parent', hFileList));
else
    % If a numeric input was supplied, do the following:
    %  -1  Load previous file in filelist
    %   0  Load next file in filelist
    %   n  Load the n'th file in the filelist
    % If varargin{2} is a string, then load that file
    if isstr(varargin{2})
        if exist(varargin{2}, 'file')
            sNextTrial = varargin{2};
        else
            sp_disp(['Failed to open ' varargin{2}])
            return;
        end
    else
        hFileList = findobj(g_hSpike_Viewer, 'Tag', 'MenuFileList');
        csPaths = flipud(get(get(hFileList, 'Children'), 'UserData'));
        switch varargin{2}
            case -1
                if isempty(csPaths)
                    warndlg('No file list has been loaded.', 'Spiky');
                    return
                end
                % Get index of current file
                nNewIndx = find(strcmp(csPaths, FV.sLoadedTrial)) - 1;
                if nNewIndx < 1
                    warndlg('You are already at the first file in the current file list.', 'Spiky')
                    return
                end
                sNextTrial = csPaths{nNewIndx};
            case 0
                if isempty(csPaths)
                    warndlg('No file list has been loaded.', 'Spiky');
                    return
                end
                nNewIndx = find(strcmp(csPaths, FV.sLoadedTrial)) + 1;
                if nNewIndx > length(csPaths)
                    warndlg('You are already at the last file in the current file list.', 'Spiky')
                    return
                end
                sNextTrial = csPaths{nNewIndx};
            otherwise
                if length(csPaths) < varargin{2}
                    warndlg('You have reached the last file in the current file list.', 'Spiky')
                    return
                end
                if iscell(csPaths)
                    sNextTrial = csPaths{varargin{2}};
                else
                    sNextTrial = csPaths(varargin{2},:);
                end
        end
    end
end

% Ask to save changes to current open file
if ~isempty(strfind(get(g_hSpike_Viewer, 'Name'), '*')) && ~isempty(FV.sLoadedTrial)
    switch questdlg('Save changes to current file?', 'Spiky', 'Yes', 'No', 'Cancel', 'Yes')
        case 'Yes', SaveResults
        case 'Cancel', return
    end
end

SetStruct(FV, 'nosaveflag'); % update current path and colormap

% Exit merge mode if applicable
if IsMergeMode(), StopMergeMode(); end

% Check for defined multitrodes to keep
csMultiTrodes = {};
csCh = fieldnames(FV.tSpikes);
for nCh = 1:length(csCh)
    if ~isempty(strfind(csCh{nCh}, '__')) % multitrode found
        csMultiTrodes{end+1} = csCh{nCh};
    end
end

% Keep copy of current settings
FV_old = FV; % keep copy of current settings
if isfield(FV_old, 'tData'), FV_old = rmfield(FV_old, 'tData'); end
if isfield(FV_old, 'tSpikes'), FV_old = rmfield(FV_old, 'tSpikes'); end
if isfield(FV_old, 'sLoadedTrial'), FV_old = rmfield(FV_old, 'sLoadedTrial'); end

% Ask to re-use current settings
if isempty(FV.sLoadedTrial)
    sAns = 'No';
else
    if isempty(p_bNeverApply), p_bNeverApply = 0; end
    if ~p_bNeverApply && ~g_bBatchMode
        sAns = questdlg('Should current settings be applied to the new file?', 'Spiky', 'Yes', 'No', 'Never apply', 'No');
    else
        sAns = 'No';
    end
    if strcmpi(sAns, 'Never apply'), p_bNeverApply = 1; sAns = 'No'; end
end

switch sAns
    case 'Yes'
        LoadTrial(sNextTrial); % load new trial data
        [FV, hWin] = GetStruct(); % get changed structure (new data/settings)
        % replace default settings with settings of previous file
        csFieldnames = fieldnames(FV_old);
        for fn = 1:length(csFieldnames)
            FV.(csFieldnames{fn}) = FV_old.(csFieldnames{fn});
        end
        FV.sLoadedTrial = sNextTrial;
        SetStruct(FV, 'nosaveflag');
        % run spike detection
        DetectSpikes();
        % concate spikes on multi-trodes
        for nCh = 1:length(csMultiTrodes)
            AddTetrode([], csMultiTrodes{nCh})
        end
    case 'No'
        % Reset to defaults
        FV = SetFVDefaults();
        if exist('sPath', 'var'), FV.sDirectory = sPath;
        else FV.sDirectory = pwd; end
        FV.sLoadedTrial = sNextTrial;
        SetStruct(FV, 'nosaveflag');
        % Load trial data. Note that the LoadTrial function will replace
        % the default settings with saved settings/results if they exist
        LoadTrial(sNextTrial);
    otherwise
        return;
end

if length(FV.sLoadedTrial ) > 70
    sLoadedTrial = ['...' FV.sLoadedTrial(end-69:end)];
else
    sLoadedTrial = FV.sLoadedTrial;
end
set(hWin, 'Name', sprintf('Spiky - %s', sLoadedTrial)) % update window name

ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImportFile(varargin)
% ImportFile retains the current data and imports data that can be
% appended. If no data is already loaded, ImportFile() runs OpenFile()
%
% ImportFile() is a batch compatible function
%
[FV, ~] = GetStruct();
global g_hSpike_Viewer g_bBatchMode;
if isempty(FV.sLoadedTrial), OpenTrial(varargin); return; end
persistent p_nDefFilter

% Decide which file to load next
if isempty(varargin{2})
    csFilters = GetImportFilters(p_nDefFilter);
    if isfield(FV, 'sDirectory')
        [sFile, sPath, nFiltIndx] = uigetfile(csFilters, 'Select data file', [FV.sDirectory filesep]);
    else
        [sFile, sPath, nFiltIndx] = uigetfile(csFilters, 'Select data file');
        set(g_hSpike_Viewer, 'UserData', FV);
    end
    if nFiltIndx == 0, return; end
    p_nDefFilter = csFilters{nFiltIndx};

    g_bBatchMode = false;
    if sFile == 0, return, end
    FV.sDirectory = sPath; % update current path
    
    SetStruct(FV, 'nosaveflag'); % update current path
    sNextTrial = sFile; % this will be the file to load
else
    % If varargin{2} is a string, then load that file
    if isstr(varargin{2})
        if exist(varargin{2}, 'file')
            sNextTrial = varargin{2};
        else
            sp_disp(['Failed to open ' varargin{2}])
            return;
        end
    else
        sp_disp(['Failed to open ' varargin{2}])
    end
end

% Exit merge mode if applicable
if IsMergeMode(), StopMergeMode(); end

% Load and store new data
LoadTrial(sNextTrial);
[FV, ~] = GetStruct();
SetStruct(FV, 'nosaveflag');
ViewTrialData();
if ~g_bBatchMode, BatchRedo([], 'ImportFile([],[])'); end % set up batch job
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV, hWin] = GetStruct(varargin)
% Fetch the FV structure from the Spiky GUI, or a sub-field of FV and the
% handle of the Spiky GUI window.
% 
% Usage:
%   [FV, hWin] = GetStruct()
%   [tData, hWin] = GetStruct('tData'), returns FV.tData
%
hWin = GetGUIHandle();
if ishandle(hWin)
    FV = get(hWin, 'UserData');
else
    FV = SetFVDefaults();
end

% Return only the requested field in FV (if any)
if nargin == 1
    if isfield(FV, varargin{1})
        FV = FV.(varargin{1});
    else
        error(['The requested field ' varargin{1} ' is not a field in FV.']);
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hWin = GetGUIHandle(varargin)
% Get the handle of the Spiky GUI window
% 
% Usage:
%   H = GetGUIHandle()
%
hWin = findobj('Tag', 'Spiky');
if length(hWin) > 1
    hWin = hWin(1);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetStruct(FV, varargin)
% varargin{1}    'nosaveflag'   Dont reset save flag (* in figure title)
global g_hSpike_Viewer

% set current working directory to that in opened file
if isfield(FV, 'sDirectory') && exist(FV.sDirectory, 'dir')
    cd(FV.sDirectory)
else return; end

if isempty(g_hSpike_Viewer), g_hSpike_Viewer = findobj('Tag', 'Spiky'); end
set(g_hSpike_Viewer, 'UserData', FV);
bSetSaveFlag = 1;
if ~isempty(varargin)
    if strcmpi(varargin, 'nosaveflag'), bSetSaveFlag = 0; end
end

iFlag = strfind(get(g_hSpike_Viewer, 'name'), '*');

if bSetSaveFlag || ~isempty(iFlag)
    set(g_hSpike_Viewer, 'Name', sprintf('%s - Spiky*', GetCurrentFile()))
else
    set(g_hSpike_Viewer, 'Name', sprintf('%s - Spiky', GetCurrentFile()))
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sLoadedFile = GetCurrentFile()
% Get name of the currently loaded file.
FV = GetStruct();
if isfield(FV, 'sLoadedTrial')
    iSep = strfind(FV.sLoadedTrial, filesep);
    if isempty(iSep)
        sLoadedFile = FV.sLoadedTrial;
    else
        sLoadedFile = FV.sLoadedTrial(iSep(end)+1:end);
    end
else sLoadedFile = ''; end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetGain(varargin)
% Manually set gain of channels in dialog window.
[FV, hWin] = GetStruct();
if ~IsDataLoaded, return, end

if isempty(FV.csChannels)
    waitfor(warndlg('There are no channels available. Cannot adjust gain.'))
    return
end

% Filter out duplicate channels
[cB, vI, vJ] = unique(FV.csChannels, 'first');
FV.csChannels = FV.csChannels(sort(vI)); % retain original order

% Set gain for each channel independently
global g_hSpike_Viewer;
if ~isfield(FV, 'tGain'), FV.tGain = struct([]); end
hFig = figure();
nCL = length(FV.csChannels)*25+15;
vHW = get(g_hSpike_Viewer, 'position');
set(hFig, 'position', [vHW(1:2) 260 nCL+25], 'menu', 'none', 'Name', 'Set Gains', 'NumberTitle', 'off', 'visible', 'off')
if exist('centerfig') centerfig(hFig, findobj('Tag', 'Spiky')); end
set(hFig, 'visible', 'on')

for nCh = 1:length(FV.csChannels)
    % channel name
    sName = FV.csChannels{nCh};
    % append channel description to its name
    if isfield(FV, 'tChannelDescriptions')
        if isfield(FV.tChannelDescriptions, 'sDescription')
            nIndx = find(strcmp({FV.tChannelDescriptions.sChannel}, sName));
            if ~isempty(nIndx)
                if ~isempty({FV.tChannelDescriptions(nIndx).sDescription})
                    sDescr = FV.tChannelDescriptions(nIndx).sDescription;
                    if ~isempty(sDescr)
                        sName = [sName ' (' sDescr ')'];
                    end
                end
            end
        end
    else
        FV.tChannelDescriptions = struct([]);
    end

    uicontrol(hFig, 'Style', 'text', 'Position', [10 nCL 180 20], ...
        'String', sName, 'HorizontalAlignment', 'left', ...
        'backgroundcolor', get(hFig, 'color'));
    if isfield(FV.tGain, FV.csChannels{nCh})
        nGainThis = FV.tGain.(FV.csChannels{nCh});
    else nGainThis = 1; end
    % input field
    uicontrol(hFig, 'Style', 'edit', 'Position', [180 nCL 70 20], 'String', num2str(nGainThis), 'HorizontalAlignment', 'center', 'backgroundcolor', [.8 .8 .8]);
    nCL = nCL - 25;
end
% OK / Cancel buttons
uicontrol(hFig, 'Style', 'pushbutton', 'Position', [120 10 55 20], 'Callback', ...
    'global g_cGains;csStrings = get(get(gcf,''children''), ''string''); g_cGains = csStrings(3:2:end-1); close(gcf)', 'String', 'OK' );
uicontrol(hFig, 'Style', 'pushbutton', 'Position', [190 10 55 20], 'Callback', ...
    'global g_mGains;g_cGains = {};close(gcf);', 'String', 'Cancel' );
uiwait % wait for user to close window or click OK/Cancel buttons
global g_cGains
if isempty(g_cGains), return, end
g_cGains = flipud(g_cGains)';
% assign gains
for nCh = 1:length(FV.csChannels)
    FV.tGain(1).(FV.csChannels{nCh}) = str2num(g_cGains{nCh});
end
clear g_cGains
SetStruct(FV)

% Plot (default)
if ~isempty(varargin)
    if strcmp(varargin{1}, 'noplot'), return
    else ViewTrialData
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetMaxSpikeJitter(varargin)
% Set the maximal duration (ms) that two spikes on different channels
% should be considered to be the same same
[FV, hWin] = GetStruct();
FV.nMaxCrossChannelSpikeJitter = GetInputNum('Set maximal spiketime jitter across electrodes (ms)', 'Spiky', FV.nMaxCrossChannelSpikeJitter);
SetStruct(FV)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EditPreTriggerTime(varargin)
% Change the duration of the signal to keep before threshold-crossings triggers
[FV, hWin] = GetStruct();
FV.nSpikeTrigPreMs = GetInputNum('Set pre-trigger time (ms)', 'Spiky', FV.nSpikeTrigPreMs);
SetStruct(FV)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetDeadtime(varargin)
% Set the length of time following threshold-crossings events that new
% spikes cannot be detected.
[FV, hWin] = GetStruct();
FV.nDeadTimeMs = GetInputNum(sprintf('Set post-spike deadtime (ms):\n(note that this parameter only affects visually displayed spikes)'), 'Spiky', FV.nDeadTimeMs);
SetStruct(FV)
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EditPostTriggerTime(varargin)
% Change the duration of the signal to keep after threshold-crossings triggers
[FV, hWin] = GetStruct();
FV.nSpikeTrigPostMs = GetInputNum('Set post-trigger time (ms)', 'Spiky', FV.nSpikeTrigPostMs);
SetStruct(FV)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nNum = GetInputNum(sMsg, sTitle, nDef)
% Get a number interactively from the user
% sMsg      Message to display
% sTitle    Window title
% nDef      Default value
sDef = {sprintf('%.3f', nDef)};
sNum = inputdlg(sMsg, sTitle, 1, sDef);
if isempty(sNum)
    nNum = nDef;
else
    sNum = cell2mat(sNum);
    nNum = str2num(sNum);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bResult = LoadTrial(snTrial)
% LoadTrial opens the selected file and imports data through the corresponding import
% filter in /import for the data format used. Additionally, previous settings and
% modifications is imported from the associated .spb file if it resides on disk.
%
% Syntax:   LoadTrial(snTrial) where snTrial is the string name of the file to load
%
[FV, ~] = GetStruct();
bResult = 0;

sFile = snTrial;
if ~exist(sFile, 'file')
    uiwait(warndlg('Data file does not exist'))
    bLoaded = 0;
    return
end

% Run import filter
iDot = strfind(sFile, '.');
eval(sprintf('FV = import_%s(sFile, FV);', lower(sFile(iDot(end)+1:end))))

% Check for errors that may have been reported during import
if isfield(FV, 'sImportError')
    sp_disp(FV.sImportError)
    uiwait(warndlg(FV.sImportError, 'Spiky::LoadTrial', 'modal'));
    return
end

% Update filename unless we are importing
sDBStack = dbstack;
bImportMode = any(strcmpi({sDBStack.name}, 'ImportFile'));
if ~bImportMode
    FV.sLoadedTrial = sFile; % uncommented 01/27/14 to fix bug #36
end

tData = FV.tData; % make backup (gets reinserted below)

% Names of all channels
cFields = fieldnames(FV.tData);
FV.csChannels = {};
for i = 1:length(cFields)
    if ~isempty(strfind(cFields{i}, '_TimeBegin'))
        FV.csChannels{end+1} = cFields{i}(1:end-10);
        if isfield(FV.tData, cFields{i}(1:end-10))
            FV.tYlim(1).(FV.csChannels{end}) = [min(FV.tData.(cFields{i}(1:end-10))) max(FV.tData.(cFields{i}(1:end-10)))];
        else FV.tYlim(1).(FV.csChannels{end}) = []; end % assume its digital
    end
end
FV.vXlim = [];
SetStruct(FV) % save settings obtained thus far

% Load associated .spb file (unless we are importing)
%sPath = [FV.sLoadedTrial(1:end-4) '.spb'];
if ~bImportMode
    sPath = [sFile(1:end-4) '.spb'];
    if exist(sPath, 'file')
        OpenSettings([sFile(1:end-4) '.spb'])
        [FV, hWin] = GetStruct();
    end
    
    % Put back tData as this was removed with the OpenSettings() call above
    if isempty(FV.tData)
        FV.tData = tData;
    end
end

% Some channels may have been digitized in another datafile and the setting
% propagated to this file via Distribute Settings. Therefore, check if any
% channels should be digitized now automatically
if ~isfield(FV, 'csDigitalChannels'), FV.csDigitalChannels = {}; end
csRemDigChannels = {};
for c = 1:length(FV.csDigitalChannels)
    sCh = FV.csDigitalChannels{c};
    if ~isfield(FV.tData, [sCh '_Up'])
        % If UP times dont exist, the try to digitize channel using stored
        % threshold value. If this fails, ignore this digital channel
        try
            nIndx = find(strcmpi({FV.tEventThresholds.sChannel}, sCh));
            nThresh = FV.tEventThresholds(nIndx).nThreshold; % threshold value
            
            vContData = FV.tData.(sCh);
            nFs = FV.tData.([sCh '_KHz']) * 1000; % sampling frequency (Hz)
            nTimeBegin = FV.tData.([sCh '_TimeBegin']); % start of sampling (sec)
            nTimeEnd = FV.tData.([sCh '_TimeEnd']); % start of sampling (sec)
            
            [vUpTimes vDownTimes] = DigitizeChannel(vContData, nThresh, nFs, nTimeBegin, nTimeEnd);
            FV.tData.([sCh '_Up']) = vUpTimes;      % UP sample times
            FV.tData.([sCh '_Down']) = vDownTimes;  % DOWN sample times
        catch csRemDigChannels{end+1} = sCh; end
    end
end
% Remove digital channels that could not be digitized above
if ~isempty(csRemDigChannels)
    for c = 1:length(csRemDigChannels)
        nIndx = find(strcmpi(FV.csDigitalChannels, csRemDigChannels{c}));
        FV.csDigitalChannels(nIndx) = [];
    end
end
SetStruct(FV) % update FV

bResult = 1;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AboutSpiky(varargin)
%
cAboutText = { 'Maintainer:', ...
    ' Per M Knutsen <p.m.knutsen@medisin.uio.no>', ...n
    ' Institute of Molecular Biology', ...
    ' Department of Physiology', ...
    ' University of Oslo, Norway', ...
    ' ', ...
    'Contributors:', ...
    ' Per M Knutsen, Weizmann Institute / UCSD / UiO', ...
    ' Idan Steinberg, Weizmann Institute', ...
    ' Amir Monovitch, Weizmann Institute', ...
    ' ', ...
    'Acknowledgements:', ...
    'Spiky includes modified code from the open-source Chronux project (http://www.chronux.org) for all of its spike-sorting routines. Chronux is released under the GNU General Public License.', ...
    ' ', ...
    'Musial PG, Baker SN, Gerstein GL, King EA, Keating JG (2002) Signal-to-noise ratio improvement in multiple electrode recording. J Neurosci Methods 115:29-43', ...
    ' ', ...
    'Fee MS, Mitra PP, Kleinfeld D (1996) Automatic sorting of multiple unit neuronal signals in the presence of anisotropic and non-Gaussian variability. J Neurosci Methods 69:175-188', ...
    };
msgbox(cAboutText, 'About Spiky')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bLoaded = IsDataLoaded(varargin)
% Check if data has been loaded into Spiky.
% 
% Usage:
%   IsDataLoaded()
%   bLoaded = IsDataLoaded()
%
[FV, hWin] = GetStruct();
if ~isfield(FV, 'tData') % check data has been loaded
    uiwait(warndlg('No data has been loaded'))
    bLoaded = 0;
    return
else bLoaded = 1; end % data has been loaded
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetSpikeThreshold(varargin)
%
[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end
zoom off
try
    [~, nY] = ginput(1);
catch
    return
end
zoom xon
sTag = get(gca, 'Tag');
vContData = FV.tData.(sTag);
[vContData, FV] = AdjustChannelGain(FV, vContData, sTag);
vContData = GetFilteredChannel(sTag, vContData);

nMed = median(vContData);
if nY == 0 || nY == nMed
    % Threshold line can't be zero!
    warndlg('Threshold cannot be zero.')
    return;
elseif nY < nMed
    FV.tSpikeThresholdsNeg(1).(sTag) = nY; % mV
else
    FV.tSpikeThresholdsPos(1).(sTag) = nY; % mV
end
SetStruct(FV)
DetectSpikes();
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vContOut, varargout] = GetFilteredChannel(sCh, vCont, varargin)
% Checks for filtering parameters and filters a channel if necessary
% 
% Usage:
%   vCont = GetFilteredChannel(sCh, vCont)
%   [vCont, vTime] = GetFilteredChannel(sCh, vCont, 'decimate')
%   [vCont, vTime, nNewFs] = GetFilteredChannel(sCh, vCont, 'decimate')
%    
% This function is a higher-level call to FilterChannel() and can
% optionally decimate. See usage.
% 
% See:
%   FilterChannel()
%

% Initialize outputs in case we abort early
vContOut = vCont;
if nargout > 1
    varargout(1:(nargout-1)) = {[]};
end

[FV, ~] = GetStruct();

if ~isfield(FV, 'tFilteredChannels'); return; end
iCh = strcmp(sCh, {FV.tFilteredChannels.sChannel});
if ~any(iCh); return; end

nHiPass = FV.tFilteredChannels(iCh).vBandpass(1);
nLoPass = FV.tFilteredChannels(iCh).vBandpass(2);
bRectify = FV.tFilteredChannels(iCh).bRectify;
nFs = FV.tData.([sCh '_KHz']) * 1000;

% Decimate
bDecimate = 0;
if ~isempty(varargin)
    if strcmpi(varargin{1}, 'decimate')
        bDecimate = 1;
    end
end

if bDecimate
    nBegin = FV.tData.([sCh '_TimeBegin']); % sampling start, sec
    nEnd = FV.tData.([sCh '_TimeEnd']); % sampling end, sec
    nFs = FV.tData.([sCh '_KHz']) * 1000; % sampling frequency Hz
    vTime = (nBegin+1/nFs):(1/nFs):(nBegin+length(vCont)/nFs); % absolute time, sec
    [vContOut, vTimeOut, nNewFs] = FilterChannel(vCont, vTime, nFs, nLoPass, nHiPass, bRectify, 'decimate');
    if nargout > 1, varargout{1} = vTimeOut; end
    if nargout > 2, varargout{2} = nNewFs; end
else
    [vContOut, ~, ~] = FilterChannel(vCont, [], nFs, nLoPass, nHiPass, bRectify, 'nodecimate');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RemoveSpikeThreshold(sCh, sDir)
%
[FV, ~] = GetStruct();
sTag = get(gco, 'Tag');
sCh = sTag(1:end-4);

% Ask whether to discard currently detected spikes
if isfield(FV.tSpikes, sCh)
    if strcmp(questdlg('Removing this threshold line will discard all detected spikes on this channel. Continue?', 'Spiky', 'Yes', 'No', 'No'), 'No')
        return;
    else FV.tSpikes = rmfield(FV.tSpikes, sCh); end
end

switch lower(sTag(end-2:end))
    case 'neg'
        FV.tSpikeThresholdsNeg = rmfield(FV.tSpikeThresholdsNeg, sCh);
    case 'pos'
        FV.tSpikeThresholdsPos = rmfield(FV.tSpikeThresholdsPos, sCh);
    otherwise
        return;
end
SetStruct(FV)
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetThresholdManually(sCh, sDir)
%
[FV, hWin] = GetStruct();
sTag = get(gco, 'Tag');
sCh = sTag(1:end-4);

% Ask for new threshold value
switch lower(sTag(end-2:end))
    case 'neg'
        cAns = inputdlg('New threshold value:', 'Spiky', 1, {num2str(FV.tSpikeThresholdsNeg.(sCh))});
        if isempty(cAns), return, end
        FV.tSpikeThresholdsNeg.(sCh) = str2num(cAns{1});
    case 'pos'
        cAns = inputdlg('New threshold value:', 'Spiky', 1, {num2str(FV.tSpikeThresholdsPos.(sCh))});
        if isempty(cAns), return, end
        FV.tSpikeThresholdsPos.(sCh) = str2num(cAns{1});
    otherwise
        return;
end
SetStruct(FV)
DetectSpikes();
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CopyPasteThreshold(sCh, sDir)
%
[FV, hWin] = GetStruct();
sTag = get(gco, 'Tag');
sCh = sTag(1:end-4);
if ~exist('clipboard') return; end
switch get(gcbo, 'label')
    case 'Copy'
        % Copy threshold into clipboard
        switch lower(sTag(end-2:end))
            case 'neg'
                nThresh = FV.tSpikeThresholdsNeg.(sCh);
            case 'pos'
                nThresh = FV.tSpikeThresholdsPos.(sCh);
        end
        clipboard('copy', nThresh);
    case 'Paste'
        % Paste threshold from clipboard
        nThresh = str2num(clipboard('paste'));
        if length(nThresh) == 1
            switch lower(sTag(end-2:end))
                case 'neg'
                    FV.tSpikeThresholdsNeg.(sCh) = nThresh;
                case 'pos'
                    FV.tSpikeThresholdsPos.(sCh) = nThresh;
            end
            % Update main window
            SetStruct(FV)
            DetectSpikes();
            ViewTrialData
        end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AutoDetectThresholds(varargin)
% Automatically select 4 x STD as negative spike threshold
%
if ~IsDataLoaded, return, end
[FV,~] = GetStruct();

% Select threshold to find thresholds on
csSelected = FV.csDisplayChannels;
hWait = waitbar(0, 'Automatically determining spike thresholds...');

for nCh = 1:length(csSelected)
    hWait = waitbar(nCh/length(csSelected), hWait);
    sCh = csSelected{nCh};
    vCont = ChannelCalculator(FV.tData.(sCh), sCh);
    [vCont, FV] = AdjustChannelGain(FV, vCont, sCh);

    % Filter channel
    vCont = GetFilteredChannel(sCh, vCont);
    
    % Negative threshold (values < 1th percentile)
    FV.tSpikeThresholdsNeg(1).(sCh) = nanmean(vCont) - (4*nanstd(vCont));
    FV.tSpikeThresholdsPos(1).(sCh) = nanmean(vCont) + (4*nanstd(vCont));
end
close(hWait)
SetStruct(FV)
ViewTrialData()
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RemoveThresholds(varargin)
% Remove existing thresholds on all channels
%
if ~IsDataLoaded, return, end
[FV,~] = GetStruct();

if isfield(FV, 'tSpikeThresholdsPos')
    FV.tSpikeThresholdsPos = struct([]);
end
if isfield(FV, 'tSpikeThresholdsNeg')
    FV.tSpikeThresholdsNeg = struct([]);
end

SetStruct(FV)
ViewTrialData()
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DetectSpikes(varargin)
% Find spikes in continuous data vector
% Function supports both negative and positive thresholds. Note, however,
% that positive thresholds are only used as an inclusion criteria;
% detection is based ONLY on the negative threshold. If a positive
% threshold is also set, then the negative threshold crossing waveform need
% to cross the positive threshold as well.
%
if ~IsDataLoaded, return, end

[FV, ~] = GetStruct();

% Iterate across all channels with threshold lines(s)
csChNeg = fieldnames(FV.tSpikeThresholdsNeg);
csChPos = fieldnames(FV.tSpikeThresholdsPos);
csChUnique = unique([csChNeg;csChPos]);
if length(csChUnique) > 1
    if ~IsMergeMode()
        bResult = SpikyWaitbar(0, length(csChUnique));
    end
end
for nCh = 1:length(csChUnique)
    sCh = csChUnique{nCh};
    if ~isfield(FV.tData, sCh) return; end
    
    % Negative threshold line
    if isfield(FV.tSpikeThresholdsNeg, sCh), vThresh(1) = FV.tSpikeThresholdsNeg.(sCh);
    else vThresh(1) = NaN; end
    
    % Positive threshold line
    if isfield(FV.tSpikeThresholdsPos, sCh), vThresh(2) = FV.tSpikeThresholdsPos.(sCh);
    else vThresh(2) = NaN; end

    % Filter continuous signal
    vCont = ChannelCalculator(FV.tData.(sCh), sCh);
    [vCont, FV] = AdjustChannelGain(FV, vCont, sCh);
    nFs = FV.tData.([sCh '_KHz']) * 1000; % channel sampling rate
    
    % Filter channel
    vCont = GetFilteredChannel(sCh, vCont);
    
    % Run spike detection
    nBeginTime = FV.tData.([sCh '_TimeBegin']) * 1000; % begin time
    tSpikes = GetWaveforms(vCont, vThresh, nFs, FV.nSpikeTrigPreMs, FV.nSpikeTrigPostMs, 0.01, nBeginTime);
    FV.tSpikes(1).(sCh) = tSpikes;
    FV.tSpikes.(sCh).dejittered = 0;

    if length(csChUnique) > 1
        if ~IsMergeMode()
            bResult = SpikyWaitbar(nCh, length(csChUnique));
            if ~bResult
                bResult = SpikyWaitbar(length(csChUnique), length(csChUnique));
                return;
            end
        end
    end
    
end
if length(csChUnique) > 1
    if ~IsMergeMode()
        bResult = SpikyWaitbar(length(csChUnique), length(csChUnique));
    end
end
SetStruct(FV)
ViewTrialData();
global g_bBatchMode
if ~g_bBatchMode, BatchRedo([], 'DetectSpikes'); end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tSpikes = GetWaveforms(vTrace, vThresh, nFS, nPreTrigMs, nPostTrigMs, nDeadTimeMs, nBeginTime)
% Extract waveforms (spikes) from a continuous signal.
% 
% Inputs:
%   vTrace        Continuous trace used for spike detection
%   vThresh       Threshold [NEG POS]
%   nFs           Sampling rate
%   nPreTrigMs    Time to keep before spike onset (ms)
%   nPostTrigMs   Time to keep after spike onset (ms)
%   nDeadTimeMs   Time after spike when no other spike can be detected (ms)
%   hWait         Absolute time of first sample (ms)
%

nPreTrig = round(nPreTrigMs * (nFS/1000)); % pre-trig duration (datapoints; rounded)
nPostTrig = round(nPostTrigMs * (nFS/1000)); % post-trig duration (datapoints; rounded)
nDeadtime = round(nDeadTimeMs * (nFS/1000)); % deadtime after spike (datapoints; rounded)

% times when signal was below threshold
vSpiketimes_Neg = find( diff(vTrace < vThresh(1)) == 1)'; % negative crossings
vSpiketimes_Pos = find( diff(vTrace > vThresh(2)) == 1)'; % positive crossings
vSpiketimes_Neg( [vSpiketimes_Neg-1] < 1) = [];
vSpiketimes_Pos( [vSpiketimes_Pos-1] < 1) = [];

% include negative crossing spikes that were lower after than before crossing thresholds
vTraceX = vTrace - vThresh(1);
vSpiketimes_Neg = vSpiketimes_Neg(vTraceX(vSpiketimes_Neg-1) > 0);

% include positive crossing spikes that were higher after than before crossing thresholds
vTraceX = vTrace - vThresh(2);
vSpiketimes_Pos = vSpiketimes_Pos(vTraceX(vSpiketimes_Pos+1) > 0);

% if positive threshold is used:
% exclude spikes that do not also cross the positive thresholds
% note that only the negative threshold is used for alignment
% the positive threshold is only an exclusive condition
vDropSpikes = [];
if ~isnan(vThresh(2)) && ~isempty(vSpiketimes_Pos)
    % condition: every negative threshold crossing need to be paired with a
    % positive threshold crossing within the window nPreTrig:nPostTrig
    mIntervals = [vSpiketimes_Neg - nPreTrig vSpiketimes_Neg + nPostTrig];
    for sn = 1:length(vSpiketimes_Neg)
        if ~any( vSpiketimes_Pos > mIntervals(sn,1) & vSpiketimes_Pos < mIntervals(sn,2) )
            vDropSpikes(end+1) = sn;
        end
    end
end

%vSpiketimes = sort([vSpiketimes_Neg; vSpiketimes_Pos]);
vSpiketimes_Neg(vDropSpikes) = [];
vSpiketimes = sort([vSpiketimes_Neg]);

% drop spiketimes close to edges (close is defined as 1.5x either the pre or
% post-threshold crossing duration)
vIndxDrop = vSpiketimes <= (nPreTrig*1.5) | vSpiketimes >= (length(vTrace)-1.5*nPostTrig);
vSpiketimes(vIndxDrop) = [];

vDeadTimeIndx = find( diff(vSpiketimes) <= nDeadtime) + 1; % crossings within deadtime
vSpiketimes(vDeadTimeIndx) = [];

% Drop spikes close to start/end of segment
vSpiketimes((vSpiketimes-nPreTrig) < 1 | (vSpiketimes+nPostTrig) > length(vTrace)) = [];

% Construct matric with all waveform indices
vRange = -nPreTrig:nPostTrig;
mIndx = [];
for i = vRange, mIndx = [mIndx vSpiketimes+i]; end
mWaveforms = vTrace(mIndx); % waveforms
tSpikes.waveforms = double(mWaveforms); % spikes must of double() type (required by sorting functions)
% Keep spiketimes after adjusting by nBeginTime
nBeginTime = nBeginTime * (nFS/1000); % samples
tSpikes.spiketimes = vSpiketimes + nBeginTime;
tSpikes.Fs = nFS;
tSpikes.threshV = [vThresh(1) inf]; % thresholds [hi lo]
tSpikes.threshT = nPreTrig + 2;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DejitterSpikes(varargin)
% Dejitter all waveforms. If spikes are unsorted, attempt to dejitter in
% bulk. TODO [[If spikes are already sorted, dejitter spikes instead on a
% unit-by-unit basis]].
%
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
csFieldnames = fieldnames(FV.tSpikes);
if length(csFieldnames) > 1
    hWait = waitbar(0, 'Dejittering spikes ...');
end
for i = 1:length(csFieldnames)
    if length(csFieldnames) > 1, waitbar(i/length(csFieldnames), hWait); end

    % Dont try to dejitter stereo/tetrodes
    if ~isempty(strfind(csFieldnames{i},'__')), continue, end
    if isempty(FV.tSpikes.(csFieldnames{i}).waveforms), continue, end
    tSpikesCh = FV.tSpikes.(csFieldnames{i});

    % Dejitter spikes unit-by-unit if already sorted
    if isfield(tSpikesCh, 'hierarchy')
        vAssigns = tSpikesCh.hierarchy.assigns;
    else
        vAssigns = zeros(size(tSpikesCh.spiketimes));
    end
    
    for a = unique(vAssigns)'
        iThis = vAssigns == a;
        tSpikes = struct([]);
        % Get waveforms of current unit
        tSpikes(1).waveforms = tSpikesCh.waveforms(iThis, :);
        tSpikes.spiketimes = tSpikesCh.spiketimes(iThis);
        tSpikes.Fs = tSpikesCh.Fs;
        tSpikes.threshV = tSpikesCh.threshV;
        tSpikes.threshT = tSpikesCh.threshT;
        
        % Skip if there are less than 10 spikes
        if any(size(tSpikes.waveforms) < 10), continue; end
        
        % in case several thresholds are given (eg. for merged files) check
        % they are all the same; if not, dejittering cannot be made
        if length(unique(tSpikes.threshV(:,1))) > 1
            sStr = ['Channel ' csFieldnames{i} ' contains more than one threshold value. Therefore, the channel cannot be dejittered.'];
            uiwait(warndlg(sStr))
            continue
        end
        tSpikes.threshV = tSpikesCh.threshV(1,:);
        tSpikes.threshT = tSpikesCh.threshT(1);
        
        % Dejitter
        try
            tSpikes = ss_dejitter(tSpikes, round(tSpikes.Fs(1)/10000)*2);
        catch
            tSpikes.waveforms = [];
        end

        if isempty(tSpikes.waveforms)
             % Dejitter failed
            sStr = sprintf('Failed to dejitter channel %s (unit %d). The error returned was:\n\n%s', ...
                csFieldnames{i}, a, lasterr);
            hFig = warndlg(sStr);
            uiwait(hFig)
        else
            % Dejitter succeeded
            if size(tSpikes.waveforms, 1) == size(tSpikesCh.waveforms, 1)
                tSpikesCh.waveforms = tSpikes.waveforms;
            else
                % Pad waveforms with NaNs
                nPad = size(tSpikesCh.waveforms, 2) - size(tSpikes.waveforms, 2);
                mPad = nan(size(tSpikes.waveforms, 1), nPad);
                tSpikesCh.waveforms(iThis, :) = [tSpikes.waveforms mPad];
            end
            sStr = ['Dejittered unit ' num2str(a) ' channel ' csFieldnames{i}];
        end
        sp_disp(sStr)
    end
    % Remove columns in tSpikesCh.waveforms that contain only NaNs
    iCols = all(isnan(tSpikesCh.waveforms), 1);
    tSpikesCh.waveforms(:, iCols) = [];

    % Replace remaining NaNs with zeros
    tSpikesCh.waveforms(isnan(tSpikesCh.waveforms)) = 0;

    % Return dejittered waveforms
    FV.tSpikes.(csFieldnames{i}) = tSpikesCh;
    FV.tSpikes.(csFieldnames{i}).dejittered = 1;
end
if length(csFieldnames) > 1, close(hWait); end
SetStruct(FV)
global g_bBatchMode
if ~g_bBatchMode BatchRedo([], 'DejitterSpikes'); end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RemoveOutlierSpikes(varargin)
% Remove waveforms that are considered outliers using K-means based outlier
% detection algorithm
%

if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
csFieldnames = fieldnames(FV.tSpikes);
if length(csFieldnames) > 1,
    hWait = waitbar(0, 'Removing outliers from all channels...');
    set(hWait, 'position', get(hWait, 'position')+[0 150 0 0])
end
for nCh = 1:length(csFieldnames)
    if length(csFieldnames) > 1, waitbar(nCh/length(csFieldnames), hWait), end
    tSpikes = FV.tSpikes.(csFieldnames{nCh});
    % Skip if there are less than 10 spikes
    if any(size(tSpikes.waveforms) < 10), continue; end
    tSpikes = ss_outliers(tSpikes, (1-1/(size(tSpikes.waveforms,2)*100))); % remove outliers
    FV.tSpikes.(csFieldnames{nCh}) = tSpikes;
end
if length(csFieldnames) > 1, close(hWait), end
SetStruct(FV)
global g_bBatchMode
if ~g_bBatchMode
    BatchRedo([], 'RemoveOutlierSpikes');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowAggregationTree(varargin)
% Show overclustered assignments and results after joining clusters
% TODO: Needs to be modified if users has joined/split clusters manually
%

if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
% Select channel from list of sorted units
csFieldnames = fieldnames(FV.tSpikes)';
vDel = [];
for fn = 1:length(csFieldnames)
    if ~isfield(FV.tSpikes.(csFieldnames{fn}), 'hierarchy'), vDel(end+1) = fn; end
end
csFieldnames(vDel) = [];
if ~CheckIfSorted, return, end
[sCh, bResult] = SelectChannelNumber(csFieldnames);
if ~bResult, return, end

hFig = figure;
ThemeObject(hFig)
set(hFig, 'name', 'Spiky Aggregration Tree', 'menubar', 'none')
aggtree(FV.tSpikes.(sCh));

hHeader = header(['Channel ' sCh ' - Spike Cluster Aggregation Tree'], 12);
ThemeObject(hHeader)
set(hHeader, 'interpreter', 'none')

SetStruct(FV)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bSorted = CheckIfSorted(varargin)
% Check whether spikes have been sorted on at least one channel/electrode
%
[FV, ~] = GetStruct();
% Select channel from list of sorted units
csFieldnames = fieldnames(FV.tSpikes)';
vDel = [];
for fn = 1:length(csFieldnames)
    if ~isfield(FV.tSpikes.(csFieldnames{fn}), 'hierarchy'), vDel(end+1) = fn; end
end
csFieldnames(vDel) = [];
if isempty(csFieldnames)
    warndlg('You must first run spike detection and sorting.')
    bSorted = 0;
else bSorted = 1; end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowSpikeClusters(varargin)
%
[FV, ~] = GetStruct();
% Select channel from list of sorted units
csFieldnames = fieldnames(FV.tSpikes)';
vDel = [];
for fn = 1:length(csFieldnames)
    if ~isfield(FV.tSpikes.(csFieldnames{fn}), 'hierarchy'), vDel(end+1) = fn; end
end
csFieldnames(vDel) = [];
if isempty(varargin{2})
    if ~CheckIfSorted, return, end
    [sCh, bResult] = SelectChannelNumber(csFieldnames);
    if ~bResult, return, end
else sCh = varargin{2}; end

% Check if this is a multitrode
if ~isempty(strfind(sCh,'__')) % is
    vIndx = strfind(sCh,'__');
    nFs = FV.tData.([sCh(1:vIndx(1)-1) '_KHz']) * 1000; % sampling frequency of 1st ch (Hz)
else
    nFs = FV.tData.([sCh '_KHz']) * 1000; % sampling frequency (Hz)
end
tSpikes = FV.tSpikes.(sCh);
vUnits = unique(tSpikes.hierarchy.assigns);
nUnits = length(vUnits);
mWaveforms = tSpikes.waveforms;
vSpiketimes = tSpikes.spiketimes; % samples

p = 10; % units per plot
nCols = min([p nUnits]);
hMenu = [];
hMainWin = figure;
ThemeObject(hMainWin)
set(hMainWin, 'name', 'Spiky Spike Clusters')

% Channel selection drop-down menu
csDescriptions = cell(size(csFieldnames));
for c = 1:length(csFieldnames)
    csDescriptions{c} = GetChannelDescription(csFieldnames{c});
end
uicontrol(hMainWin, 'Style', 'popupmenu', 'units', 'normalized', ...
    'Position', [0 .95 .2 .05], 'String', ...
    csDescriptions, 'userdata', csFieldnames, 'callback', ...
    'cUD=get(gcbo,''userdata'');nV=get(gcbo,''value'');close(gcf);Spiky.main.ShowSpikeClusters([],cUD{nV})')
nGroup = 0;

% Scale window width to number of columns
vPos = get(hMainWin, 'position');
if nUnits >= 10, set(hMainWin, 'position', vPos .* [1 1 2 1]), end

% Calculate number of spikes to show per unit (total number should not exceed 5,000)
nMaxSpikesTotal = 5000;
nDisplaySpikesPercent = min([1 nMaxSpikesTotal / length(tSpikes.hierarchy.assigns)]); % percent spikes to show from each unit

for u = 1:nUnits % one row per unit
    if ~ishandle(hMainWin) return; end
    
    % Create new uipanel for next group of units
    if p == 10
        p = 1;
        hCurrPanel = uipanel('Position', [0 0 1 1], 'parent', hMainWin);
        ThemeObject(hCurrPanel)

        % Panel selection button
        nGroup = nGroup + 1;
        sCmd = 'h=findobj(gcf,''type'',''uipanel'');set(h,''visible'',''off'');set(get(gcbo,''userdata''),''visible'',''on'')';
        uicontrol(hMainWin, 'Style', 'pushbutton', 'units', 'normalized', ...
            'Position', [.2+((nGroup-1)*.125) .95 .125 .05], ...
            'String', sprintf('Group %d', nGroup), 'callback', sCmd, ...
            'userdata', hCurrPanel )
    else p = p + 1; end
    
    % Create context menu that will attach to each subplot for current unit
    hMenu(end+1) = uicontextmenu;
    tUserData = struct([]);
    tUserData(1).sChannel = sCh;
    tUserData(1).nUnit = vUnits(u);
    hItem(1) = uimenu(hMenu(end), 'Label', sprintf('&Classify Unit %d', vUnits(u)), 'Callback', @ClassifyUnit, 'UserData', tUserData);
    
    % Waveform indices of this unit
    vIndx = find(tSpikes.hierarchy.assigns == vUnits(u));
    
    % Keep just 1000 waveforms, and discard the rest
    %vIndxKeep = randperm(length(vIndx));
    %vIndxKeep = vIndxKeep(1:min([1000 length(vIndx)]));
    % Calculate number of spikes to display
    nDisplaySpikesCount = round(nDisplaySpikesPercent * length(vIndx));
    vIndxKeep = unique(round(linspace(1, length(vIndx), nDisplaySpikesCount)));
    
    % Always show spikes with extreme amplitudes
    vSums = sum(abs(mWaveforms(vIndx,:)),2);
    vIndxKeep = [vIndxKeep find(vSums > nanmean(vSums)+(3*nanstd(vSums)))'];
    
    %   1) Plot waveforms of unit
    nX = .05; nY = .08;
    nW = .90; nH = .8;
    nSpace = 0;
    hWaveAx = axes('position', [nX+(p-1)*(nW/nCols)+nSpace nY+(nH/3)*2 nW/nCols-nSpace nH/3], 'parent', hCurrPanel);
    set(hWaveAx, 'uicontextmenu', hMenu(end))
    
    vYY = mWaveforms(vIndx(vIndxKeep),:)';
    vYY = detrend(vYY);
    vXX = repmat([([1:size(vYY,1)]/nFs)*1000]', 1, size(vYY, 2));
    if vUnits(u) == 0, mCol = [.5 .5 .5]; % outlier
    else mCol = FV.mColors(u,:); end
    plot(vXX, vYY, '-', 'color', mCol, 'linewidth', .5)
    hold on
    hLin = plot([([1:size(vYY,1)]/nFs)*1000], median(vYY'), '--');
    ThemeObject(hLin)
    axis tight
    ThemeObject(hWaveAx)
    hTit = title(['Unit# ' num2str(vUnits(u)) ' - ' num2str(length(vIndx)) ' spikes']);
    ThemeObject(hTit)
    if p>1, set(hWaveAx, 'xtick', [])
    else xlabel('ms'); ylabel('mV'); end
    
    % If waveforms axis is clicked, change value of checkbox in the Cluster Control window (if open)
    % tag: ClusterControlFig_DAQ_3
    sCmd = ['hf=findobj(''tag'',''ClusterControlFig_' sCh ''');h=findobj(hf,''string'',' num2str(vUnits(u)) ',''style'',''checkbox'');set(h,''value'',~get(h,''value''));sCB=get(h,''callback'');figure(hf);sCB(h)'];
    set(hWaveAx, 'buttondownfcn', sCmd);
    
    %   2) Plot 2D histogram of unit
    hHistAx = axes('position', [nX+(p-1)*(nW/nCols)+nSpace nY+nH/3 nW/nCols-nSpace nH/3], 'parent', hCurrPanel);
    set(hHistAx, 'uicontextmenu', hMenu(end))
    hist2d(detrend(mWaveforms(vIndx, :)')', 200);
    ThemeObject(hHistAx)
    set(hHistAx, 'xtick', [], 'ytick', [])
    hTxt = text(0,0,'Detrended');
    ThemeObject(hTxt)
    set(hTxt, 'units', 'normalized', 'position', [.01 .90 0])
    
    %   3) ISI histogram
    hAx = axes('position', [nX+(p-1)*(nW/nCols)+nSpace nY nW/nCols-nSpace nH/3], 'parent', hCurrPanel);
    set(hAx, 'uicontextmenu', hMenu(end))

    nMaxISI = 50; % ms
    vTheseSpiketimes = DropDeadtimeSpikes((vSpiketimes(vIndx) ./ nFs));
    [vN, vX] = GetISIHistogram(vTheseSpiketimes, nMaxISI / 1000);
    vX = vX .* 1000; % change from s to ms
    hBar = bar(vX, vN);
    ThemeObject(hBar)
    hold on
    
    % Plot deadtime
    plot([FV.nDeadTimeMs FV.nDeadTimeMs], [0 max(vN)+10], 'r:')
    ThemeObject(hAx)
    set(hAx, 'xtick', [0:10:max(vX)], 'xlim', [0 nMaxISI/2], ...
        'ylim', [0 max(vN)+1])
    if p==1 ylabel('N'); end
    xlabel('ISI (ms)')
    
    drawnow
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the ISI histogram of a spike-train
%   GetISIHistogram(S, R) where S is spiketimes in seconds and R is
%   the max range of the histogram, e.g. 0.05 s
% 
%   Returned results are in seconds.
%
function [vN, vT] = GetISIHistogram(vSpiketimes, nMaxRange)

vSpiketimes = vSpiketimes .* 1000;  % change from s to ms
nMaxRange = nMaxRange * 1000;       % change from s to ms

% Decide what width to use for bars in histograms
nSpkLen = length(find(diff(vSpiketimes) <=50 ));
if nSpkLen > 2500, nStep = 0.25;
elseif nSpkLen > 500, nStep = 0.5;
else nStep = 1; end

vRange = 0:nStep:nMaxRange;
[vN, vT] = hist(diff(vSpiketimes), [vRange vRange(end)+1]); % msec, linear scale
vN = vN(1:end-1); % spikes/bin
vT = vT(1:end-1) ./ 1000; % change from ms to s

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowOverlappingSpikeClusters(varargin)
% Plot overlapping, color-coded unit/spike clusters from a single
% electrode/multi-trode
%
if ~IsDataLoaded, return, end
[FV,hWin] = GetStruct();
if ~isempty(varargin)
    if isobject(varargin{1}), sCh = '';
    elseif isnumeric(varargin{1}), sCh = '';
    else sCh = varargin{1}; end
else sCh = ''; end

% Select channel from list of sorted units
csFieldnames = fieldnames(FV.tSpikes)';
vDel = [];
for fn = 1:length(csFieldnames)
    if ~isfield(FV.tSpikes.(csFieldnames{fn}), 'hierarchy'), vDel(end+1) = fn; end
end
csFieldnames(vDel) = [];
if isempty(sCh)
    if ~CheckIfSorted, return, end
    [sCh, bResult] = SelectChannelNumber(csFieldnames);
    if ~bResult, return, end
end

% Check if this is a multitrode
if ~isempty(strfind(sCh,'__')) % is
    vIndx = strfind(sCh,'__');
    nFs = FV.tData.([sCh(1:vIndx(1)-1) '_KHz']) * 1000; % sampling frequency of 1st ch (Hz)
else
    nFs = FV.tData.([sCh '_KHz']) * 1000; % sampling frequency (Hz)
end
tSpikes = FV.tSpikes.(sCh);
vUnits = unique(tSpikes.hierarchy.assigns);
mWaveforms = tSpikes.waveforms;

% Open new figure if a windows is not already open for this channel
sFigTitle = 'Spiky Cluster Control';
sFigTag = sprintf('ClusterControlFig_%s', sCh);
sAxTag = sprintf('ClusterControlAx_%s', sCh);
hMainWin = findobj('Tag', sFigTag);
hAx = findobj('Tag', sAxTag);
if isempty(hMainWin) && isempty(hAx)
    hMainWin = figure;
    hAx = axes('position', [.1 .1 .74 .8], 'Tag', sAxTag);
else
    figure(hMainWin);
    axes(hAx);
end
ThemeObject(hMainWin)
set(hMainWin, 'Tag', sFigTag, 'name', sFigTitle, 'MenuBar', 'figure')
hold on
ThemeObject(hAx)
set(hAx, 'Tag', sAxTag)
drawnow
% Remove all checkboxes from window
delete(findobj(hMainWin, 'Style', 'checkbox'));

sCmd = 'vH=findobj(gcf,''Style'',''checkbox'');for i=vH(1:end-1)'',set(i,''value'',~get(i,''value''));sCB=get(i,''callback'');sCB(i);end';

hCtrl = uicontrol(hMainWin, 'Style', 'checkbox', 'units', 'normalized', 'String', 'All', ...
    'Position', [.85 .88 .15 .05], 'callback', [''], 'HorizontalAlignment', 'left', ...
    'value', 1, 'callback', sCmd);
ThemeObject(hCtrl)

% Calculate number of spikes to show per unit
% The total number of displayed spikes should not exceed 5,000
% If there are less than 5,000 spikes, then show all.
nMaxSpikesTotal = 5000;
nDisplaySpikesPercent = min([1 nMaxSpikesTotal / length(tSpikes.hierarchy.assigns)]); % percent spikes to show from each unit

nMaxDisplayUnits = 20; % max # of units to display
for u = 1:min([length(vUnits) nMaxDisplayUnits]) % one row per unit
    % Line color
    if vUnits(u) == 0, mCol = [.2 .2 .2]; mColTxt = [.4 .4 .4]; % outliers
    else mCol = FV.mColors(u,:); mColTxt = FV.mColors(u,:); end

    if ~ishandle(hMainWin), return, end % exit if window was closed
    
    % Create new line objects for this unit if none already exist
    sLineTag = num2str(vUnits(u));
    hLines = findobj(hAx, 'Tag', sLineTag);
    if isempty(hLines)
        vIndx = find(tSpikes.hierarchy.assigns == vUnits(u));
        vIndxKeep = randperm(length(vIndx));
        
        % Calculate number of spikes to display
        nDisplaySpikesCount = round(nDisplaySpikesPercent * length(vIndx));
        if nDisplaySpikesCount < 500
            nDisplaySpikesCount = length(vIndx);
        end
        vIndxKeep = unique(round(linspace(1, length(vIndx), nDisplaySpikesCount)));
        
        % Always show the largess spikes
        vSums = sum(abs(mWaveforms(vIndx,:)),2);
        vIndxKeep = [vIndxKeep find(vSums > nanmean(vSums)+(3*nanstd(vSums)))'];
        
        vYY = mWaveforms(vIndx(vIndxKeep),:)';
        vXX = repmat([([1:size(vYY,1)]/nFs)*1000]', 1, size(vYY, 2));
        if ~ishandle(hMainWin), return, end % exit if window was closed
        figure(hMainWin);
        hLines = plot(vXX, vYY, '-', 'color', mCol, 'linewidth', .5);
        set(hLines, 'tag', sLineTag);
        
        % Hide outliers by default
        if vUnits(u) == 0, set(hLines, 'visible', 'off'); end
    else
        set(hLines, 'color', mCol)
    end

    if vUnits(u) == 0, sUnitName = 'Outliers';
    else sUnitName = [num2str(vUnits(u))]; end
    nCheckBoxHeight = .0375;
    nMaxBoxesPerColumn = round(.88 / nCheckBoxHeight);
    hCheckbox = uicontrol(hMainWin, 'Style', 'checkbox', 'units', 'normalized', 'String', sUnitName, 'HorizontalAlignment', 'left', 'value', 1);
    ThemeObject(hCheckbox)
    set(hCheckbox, 'foregroundColor', mColTxt)
    if u <= nMaxBoxesPerColumn
        set(hCheckbox, 'Position', [.85 .88-([u]*nCheckBoxHeight) .15 .05], 'Tag', num2str(vUnits(u)), 'callback', @ToggleWaveforms);
    else
        set(hCheckbox, 'Position', [.92 .88-([u-nMaxBoxesPerColumn-1]*nCheckBoxHeight) .15 .05], 'Tag', num2str(vUnits(u)), 'callback', @ToggleWaveforms);
    end
    if ~isempty(hLines)
        switch get(hLines(1), 'visible')
            case 'on'
                set(hCheckbox, 'value', 1)
            case 'off'
                set(hCheckbox, 'value', 0)
        end
    end
    xlabel('ms'); ylabel('mV')
    hold on
    drawnow
end

% Plot threshold line(s)
cFields = fieldnames(FV.tSpikeThresholdsNeg);
cFieldsPos = fieldnames(FV.tSpikeThresholdsPos);

nIndx = find(strcmpi(cFields, sCh));
nIndxPos = find(strcmpi(cFieldsPos, sCh));

if ~ishandle(hAx); return; end
if ~isempty(FV.tSpikeThresholdsNeg) && ~isempty(fieldnames(FV.tSpikeThresholdsNeg)) && ~isempty(nIndx) % negative threshold
    nT_neg = FV.tSpikeThresholdsNeg.(cFields{nIndx});
    hLin = plot(get(hAx, 'xlim'), [nT_neg nT_neg], '--', 'Tag', 'threshold_neg');
    ThemeObject(hLin)
end
if ~isempty(FV.tSpikeThresholdsPos) && ~isempty(fieldnames(FV.tSpikeThresholdsPos))  && ~isempty(nIndxPos) % positive threshold
    nT_pos = FV.tSpikeThresholdsPos.(cFieldsPos{nIndxPos});
    hLin = plot(get(hAx, 'xlim'), [nT_pos nT_pos], '--', 'Tag', 'threshold_pos');
    ThemeObject(hLin)
end

if ~ishandle(hMainWin); return, end

uicontrol(hMainWin, 'Style', 'popupmenu', 'units', 'normalized', ...
    'Position', [.1 .85 .2 .05], 'String', ...
    {'Reassign spike' 'Reassign group' 'Assign as outlier' 'Merge clusters' ...
    'Merge visible' 'Split cluster' 'Undo last action'}, 'callback', @ClusterControl)

uicontrol(hMainWin, 'Style', 'text', 'units', 'normalized', ...
    'Tag', 'StatusField', 'Position', [.31 .855 .53 .04], 'backgroundcolor', [.1 .1 .1], ...
    'foregroundcolor', 'w', 'horizontalalignment', 'left', 'visible', 'off')

hHeader = header(['Clustered spikes: ' GetChannelDescription(sCh)], 12);
ThemeObject(hHeader)
set(hHeader, 'interpreter', 'none')

% Pop-up menu
csDescriptions = cell(size(csFieldnames));
for c = 1:length(csFieldnames)
    csDescriptions{c} = GetChannelDescription(csFieldnames{c});
end

uicontrol(hMainWin, 'Style', 'popupmenu', 'units', 'normalized', ...
    'Position', [0 .95 .2 .05], 'String', ...
    csDescriptions, 'userdata', csFieldnames, 'callback', ...
    'cUD=get(gcbo,''userdata'');nV=get(gcbo,''value'');close(gcf);Spiky.main.ShowOverlappingSpikeClusters('' '''''''' cUD{nV} '''''''' ''))')

% Set the click Callback function
hCursor = datacursormode(hMainWin);
datacursormode on
set(hCursor, 'UpdateFcn', @ShowCurrentLine, 'SnapToDataVertex', 'on');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClusterControl(varargin)
% UI for moving spikes between clusters, assign spikes as outliers manually
%
% Usage:
%   ClusterControl()
% 
[FV, ~] = GetStruct();
persistent FV_this

hMenu = gco;
hFig = get(gco, 'Parent');
cStrings = get(hMenu, 'string');
nSelected = get(hMenu, 'value');
sSelected = cStrings{nSelected};
hStatus = findobj(hFig, 'Tag', 'StatusField');
sChannel = get(hFig, 'tag');
sChannel = sChannel(19:end); % minus ClusterControlFig

switch sSelected
    case 'Assign as outlier'
        FV_this = FV;
        % Reclassify a spike as belonging to the outlier cluster
        hLine = findobj(gcf, 'type', 'line', 'marker', 's');
        if isempty(hLine) % check that a spike is selected
            set(hStatus, 'string', 'Select a spike before reclassifying', 'visible', 'on')
        else
            sTags = get(hLine, 'tag');
            if length(sTags) > 2
                warndlg('Spikes must be assigned as outliers individually.', 'Spiky')
                return
            end
            nSpikeCluster = unique(str2num(sTags)); % cluster selected spike(s)
            vYdata = get(hLine, 'ydata');
            nNewCluster = 0;
            vIndx = find(FV.tSpikes.(sChannel).hierarchy.assigns == nSpikeCluster);
            mWaveforms = FV.tSpikes.(sChannel).waveforms(vIndx, :);
            mYdata = repmat(vYdata, size(mWaveforms,1), 1);
            nIndx = vIndx(all(mYdata == mWaveforms, 2));
            FV.tSpikes.(sChannel).hierarchy.assigns(nIndx) = nNewCluster;
            SetStruct(FV)
            close(hFig)
            ShowOverlappingSpikeClusters(sChannel);
            return
        end
        
    case 'Reassign spike'
        FV_this = FV;
        % Reclassify a spike as belonging to a different cluster
        hLine = findobj(gcf, 'type', 'line', 'marker', 's');
        if isempty(hLine) % check that a spike is selected
            set(hStatus, 'string', 'Select a spike before reclassifying', 'visible', 'on')
        else
            nSpikeCluster = str2num(get(hLine, 'tag')); % selected spike
            vYdata = get(hLine, 'ydata');
            set(hStatus, 'string', 'Select a spike from the new cluster', 'visible', 'on')
            waitfor(figure('visible', 'off', 'Tag', 'UIHoldFigure')) % Create an invisible figure
            hLine = findobj(gcf, 'type', 'line', 'marker', 's');
            nNewCluster = str2num(get(hLine, 'tag'));
            %if strcmp(questdlg(sprintf('Move spike from unit %d to unit %d?', nSpikeCluster, nNewCluster), 'Spiky', 'OK', 'Cancel', 'OK'), 'OK')
            vIndx = find(FV.tSpikes.(sChannel).hierarchy.assigns == nSpikeCluster);
            mWaveforms = FV.tSpikes.(sChannel).waveforms(vIndx, :);
            mYdata = repmat(vYdata, size(mWaveforms,1), 1);
            nIndx = vIndx(all(mYdata == mWaveforms, 2));
            FV.tSpikes.(sChannel).hierarchy.assigns(nIndx) = nNewCluster;
            SetStruct(FV)
            close(hFig)
            ShowOverlappingSpikeClusters(sChannel);
            return
            %end
        end

    case 'Reassign group'
        % Reclassify all spikes that intersect a drawn line
        FV_this = FV;

        % Get handles of all visible waveforms
        vAx = findobj(get(hFig,'children'),'type','axes');
        axes(vAx(1))
        hHandles = findobj(gca, 'visible', 'on', 'type', 'line'); % visible spikes

        set(hStatus, 'string', 'Draw a line through spikes to include', 'visible', 'on')
        
        % Get the drawn line interactively
        hold on; hLin = plot(NaN,NaN);
        datacursormode off
        
        global g_pnt1 g_pnt2 g_lin
        waitforbuttonpress
        g_pnt1 = get(gca,'CurrentPoint'); % button down detected
        g_lin = plot([g_pnt1(3) g_pnt1(3)], [g_pnt1(3) g_pnt1(3)], '--o', 'linewidth', 2);
        ThemeObject(g_lin)
        set(gcf, 'windowButtonMotionFcn', 'global g_pnt1 g_pnt2 g_lin; g_pnt2 = get(gca,''CurrentPoint''); set(g_lin, ''xdata'', [g_pnt1(1) g_pnt2(1)], ''ydata'', [g_pnt1(3) g_pnt2(3)])')
        waitforbuttonpress
        set(gcf, 'windowButtonMotionFcn', '')        
        vX_ln = get(g_lin, 'xdata');
        vY_ln = get(g_lin, 'ydata');
        delete(g_lin)

        datacursormode on
        
        % Iterate over spikes
        SpikyWaitbar(0, length(hHandles));
        hSelectedHandles = [];
        for h = 1:length(hHandles)
            SpikyWaitbar(h, length(hHandles));
            % Get current waveform
            vX_spk = get(hHandles(h), 'xdata');
            vY_spk = get(hHandles(h), 'ydata');
            if any(isnan(vX_spk)), continue, end
            % Iterate over line segments comprising spike waveform
            for i = 1:1:length(vX_spk)-1
                vX_ln2 = vX_spk(i:i+1);
                vY_ln2 = vY_spk(i:i+1);
                % Check if lines intersect
                A = [(vY_ln(2) - vY_ln(1)), -(vX_ln(2) - vX_ln(1));...
                    (vY_ln2(2) - vY_ln2(1)), -(vX_ln2(2) - vX_ln2(1)) ];

                b = [(vY_ln(2) - vY_ln(1))*vX_ln(1) - (vX_ln(2) - vX_ln(1))*vY_ln(1);...
                    (vY_ln2(2) - vY_ln2(1))*vX_ln2(1) - (vX_ln2(2) - vX_ln2(1))*vY_ln2(1) ];
                if cond(A)
                    %PInt = inv(A)*b; % intersection point
                    PInt = A\b; % intersection point
                    if (       PInt(1) >= min(vX_ln(1),vX_ln(2))) ...
                            && (PInt(1) >= min(vX_ln2(1),vX_ln2(2))) ...
                            && (PInt(1) <= max(vX_ln(1),vX_ln(2))) ...
                            && (PInt(1) <= max(vX_ln2(1),vX_ln2(2)))
                        hSelectedHandles(end+1) = hHandles(h);
                        continue
                    end
                end
            end
        end
        SpikyWaitbar(length(hHandles), length(hHandles));

        if ~isempty(hSelectedHandles)
            % De-select threshold lines
            hSelectedHandles(strcmp(get(hSelectedHandles, 'tag'), 'threshold_neg')) = [];
            hSelectedHandles(strcmp(get(hSelectedHandles, 'tag'), 'threshold_pos')) = [];
            
            % Highlight intersecting waveforms
            set(hSelectedHandles, 'marker', 'square', 'markersize', 4, 'MarkerEdgeColor', [.8 .8 .8], 'LineWidth', 2)
            
            vSpikeClusters = str2num(cell2mat(get(hSelectedHandles, 'tag'))); % unit IDs of selected spikes
            set(hStatus, 'string', 'Select a spike from the new cluster', 'visible', 'on')
            waitfor(figure('visible', 'off', 'Tag', 'UIHoldFigure')) % Create an invisible figure
            hLine = findobj(gcf, 'type', 'line', 'marker', 's');
            nNewCluster = str2num(get(hLine, 'tag'));
            for h = 1:length(hSelectedHandles)
                vYdata = get(hSelectedHandles(h), 'ydata');
                vIndx = find(FV.tSpikes.(sChannel).hierarchy.assigns == vSpikeClusters(h));
                mWaveforms = FV.tSpikes.(sChannel).waveforms(vIndx, :);
                mYdata = repmat(vYdata, size(mWaveforms,1), 1);
                nIndx = vIndx(all(mYdata == mWaveforms, 2));
                FV.tSpikes.(sChannel).hierarchy.assigns(nIndx) = nNewCluster;
            end
            SetStruct(FV)
            close(hFig)
            ShowOverlappingSpikeClusters(sChannel);
            return
        end

    case 'Merge clusters'
        FV_this = FV;
        % Merge two clusters
        set(hStatus, 'string', 'Click on spike from cluster to merge from', 'visible', 'on')
        try
            waitfor(figure('visible', 'off', 'Tag', 'UIHoldFigure')) % Create an invisible figure
        catch
            return
        end
        hLine = findobj(gcf, 'type', 'line', 'marker', 's');
        nSpikeClusterFrom = str2double(get(hLine, 'tag'));
        
        set(hStatus, 'string', 'Click on spike from cluster to merge to', 'visible', 'on')
        waitfor(figure('visible', 'off', 'Tag', 'UIHoldFigure')) % Create an invisible figure
        hLine = findobj(gcf, 'type', 'line', 'marker', 's');
        nSpikeClusterTo = str2double(get(hLine, 'tag'));
        if nSpikeClusterFrom == nSpikeClusterTo, warndlg('You cant merge a cluster into itself.', 'Spiky'); return, end

        if nSpikeClusterTo == 0 || nSpikeClusterFrom == 0
            warndlg('You cannot merge outliers with sorted clusters. You can only reassign individual outliers to existing clusters.')
            return
        end

        if strcmp(questdlg(sprintf('Merge cluster %d into %d?', nSpikeClusterFrom, nSpikeClusterTo), 'Spiky', 'OK', 'Cancel', 'OK'), 'OK')
            vFromIndx = find(FV.tSpikes.(sChannel).hierarchy.assigns == nSpikeClusterFrom);
            FV.tSpikes.(sChannel).hierarchy.assigns(vFromIndx) = nSpikeClusterTo;
            FV.tSpikes.(sChannel).hierarchy.tree(end+1,:) = [nSpikeClusterTo nSpikeClusterFrom .5 0];
            % Delete clicked objects and associated checkboxes in open figure
            vObjFrom = findobj(gcf, 'Tag', num2str(nSpikeClusterFrom));
            vObjTo = findobj(gcf, 'Tag', num2str(nSpikeClusterTo));
            delete([vObjFrom; vObjTo])
            SetStruct(FV)
            ShowOverlappingSpikeClusters(sChannel);
            return
        end

    case 'Merge visible'
        % Merge all visible clusters
        FV_this = FV;
        % Get list of visible clusters
        hLine = findobj(gcf, 'type', 'line', 'visible', 'on');
        sVisibleClusters = str2double(unique(get(hLine, 'tag')));
        vNotNaNIndx = find(~isnan(sVisibleClusters));
        nSpikeClusterTo = sVisibleClusters(vNotNaNIndx(1)); % merge into first cluster

        %if strcmp(questdlg(sprintf('Merge all visible clusters into cluster %d?', nSpikeClusterTo), 'Spiky', 'OK', 'Cancel', 'OK'), 'OK')
        for i = 2:length(sVisibleClusters)
            nSpikeClusterFrom = sVisibleClusters(i);
            vFromIndx = find(FV.tSpikes.(sChannel).hierarchy.assigns == nSpikeClusterFrom);
            FV.tSpikes.(sChannel).hierarchy.assigns(vFromIndx) = nSpikeClusterTo;
            FV.tSpikes.(sChannel).hierarchy.tree(end+1,:) = [nSpikeClusterTo nSpikeClusterFrom .5 0];
            % Delete clicked objects and associated checkboxes in open figure
            vObjFrom = findobj(gcf, 'Tag', num2str(nSpikeClusterFrom));
            vObjTo = findobj(gcf, 'Tag', num2str(nSpikeClusterTo));
            delete([vObjFrom; vObjTo])
        end
        SetStruct(FV)
        ShowOverlappingSpikeClusters(sChannel);
        return
        %end
        
    case 'Split cluster'
        FV_this = FV;
        % Split selected cluster into its two parent clusters, if they exist
        set(hStatus, 'string', 'Select cluster to split', 'visible', 'on')
        waitfor(figure('visible', 'off', 'Tag', 'UIHoldFigure')) % Create an invisible figure
        hLine = findobj(gcf, 'type', 'line', 'marker', 's');
        nSpikeClusterSplit = str2num(get(hLine, 'tag'));
        mTree = FV.tSpikes.(sChannel).hierarchy.tree;
        if isempty(mTree)
            sStr = 'This cluster does not have parents clusters and cannot be split.';
            warndlg(sStr)
            sp_disp(sStr)
            return
        end
        if nSpikeClusterSplit == 0 % outliers cant be split
            sStr = 'Outliers have not been sorted, do not belong to a cluster and cannot be split. You can only reassign individual outliers to existing clusters.';
            warndlg(sStr)
            sp_disp(sStr);
            return
        end
        vRowIndx = find(mTree(:,1) == nSpikeClusterSplit);
        if isempty(vRowIndx)
            sStr = 'This cluster does not have parents clusters and cannot be split.';
            warndlg(sStr)
            sp_disp(sStr)
            return
        elseif strcmp(questdlg(sprintf('Split cluster %d into its parent clusters %d and %d?', nSpikeClusterSplit, nSpikeClusterSplit, mTree(vRowIndx(end), 2)), 'Spiky', 'OK', 'Cancel', 'OK'), 'OK')
            nNewCluster = mTree(vRowIndx(end), 2);
            vOverclusterAssigns = FV.tSpikes.(sChannel).overcluster.assigns;
            vIndxChange = find(vOverclusterAssigns == nNewCluster);
            FV.tSpikes.(sChannel).hierarchy.assigns(vIndxChange) = nNewCluster;
            FV.tSpikes.(sChannel).hierarchy.tree(vRowIndx(end),:) = [];
            SetStruct(FV)
            close(hFig)
            ShowOverlappingSpikeClusters(sChannel);
            return
        end

    case 'Undo last action'
        if isempty(FV_this), return, end
        FV = FV_this;
        SetStruct(FV)
        close(hFig)
        ShowOverlappingSpikeClusters(sChannel);
        return
end

if ishandle(hStatus)
    set(hStatus, 'string', '', 'visible', 'off')
end
return

%%%%%%%%%%%%%%%%%%%%%%s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ToggleWaveforms(varargin)
% Toggle visibility of a selected group of waveforms in 'Cluster Control'
% window
%
if ~isempty(varargin) hCallbackHandle = varargin{1};
else hCallbackHandle = gco; end % callback handle
sTag = get(hCallbackHandle, 'Tag'); % tag of all objects to toggle
hToggleHandles = findobj(gcf, 'Tag', sTag); % find handles to toggle
hToggleHandles(find(hToggleHandles == hCallbackHandle)) = []; % remove callback handle from list
if get(hCallbackHandle, 'value'), set(hToggleHandles, 'visible', 'on');
else set(hToggleHandles, 'visible', 'off'); end
% Set YLIM to min/max of all visible lines
% - get all line objects
hVisibleLines = findobj(findobj(gcf, 'Type', 'line'), 'flat', 'Visible', 'on');
if ~isempty(hVisibleLines)
    cYData = get(hVisibleLines, 'ydata');
    if ~iscell(cYData)
        cYData = {cYData};
    end
    vRem = [];
    for c = 1:length(cYData)
        if isnan(cYData{c}), vRem(end+1) = c; end
    end
    cYData(vRem) = [];
    try
        mLines = cell2mat(cYData);
    catch, return, end
    nMin = min(mLines(:));
    nMax = max(mLines(:));
    if nMin ~= nMax
        set(findobj(gcf, 'Type', 'axes'), 'ylim', [nMin nMax])
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SortSpikes(varargin)
% Sort spikes on selected channel
if ~IsDataLoaded, return, end
[FV, hWin] = GetStruct();
if isempty(FV.tSpikes)
    warndlg('You must first run spike detection.')
    return
end
[sCh, bResult] = SelectChannelNumber(fieldnames(FV.tSpikes)');
if ~bResult, return, end
tSpikes = FV.tSpikes.(sCh);
tSpikes.Fs = tSpikes.Fs(1);
if isempty(tSpikes.waveforms) % abort if there are no spikes
    uiwait(warndlg(['There are no spikes associated with the selected channel.']))
    return
end
% Remove fields that will interfere with sorting
if isfield(tSpikes, 'hierarchy'), tSpikes = rmfield(tSpikes, 'hierarchy'); end
if isfield(tSpikes, 'overcluster'), tSpikes = rmfield(tSpikes, 'overcluster'); end
if isfield(tSpikes, 'quality'), tSpikes = rmfield(tSpikes, 'quality'); end
% Over-clustering
tSpikes = ss_kmeans(tSpikes);
% Aggregation
tSpikes = ss_energy(tSpikes);
tSpikes = ss_aggregate(tSpikes);
% Save results
FV.tSpikes.(sCh) = tSpikes;
SetStruct(FV)
% Show sorted units
ShowSpikeClusters([], sCh)
figure(hWin)
ViewTrialData();
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display list of all available keyboard shortcuts
function KeyboardShortcuts(varargin)
[~, hWin] = GetStruct();
hMenuItems = findobj(hWin, 'type', 'uimenu');
cData(:, 1) = get(hMenuItems, 'Label');
cData(:, 2) = get(hMenuItems, 'Accelerator');
vValid = ~cellfun(@isempty, cData(:, 2));
hFig = figure();
set(hFig, 'Name', 'Spiky Keyboard Shortcuts', 'ToolBar', 'none', 'menuBar','none')
ThemeObject(hFig)
cColumnNames = {'Label', 'Accelerator'};
cColumnFormat = {{'char' 'fixed'}, {'char' 'fixed'}};
cColumnEditable =  [false false];
hTable = uitable('Units', 'normalized', 'Position', [0 0 1 1], 'Data', flipud(cData(vValid, :)), ...
    'ColumnName', cColumnNames, 'ColumnEditable', cColumnEditable, 'ColumnWidth', {400, 100});
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SortAllChannels(varargin)
% Run spike sorting on all thresholded channels
if ~IsDataLoaded, return, end
bSortGood = 0;
[FV, ~] = GetStruct();
if isempty(FV.tSpikes)
    uiwait(warndlg('You must first run spike detection.'))
    return
end
% iterate over thresholded channels
csChannels = fieldnames(FV.tSpikes);
if length(csChannels) > 1
    hWait = waitbar(0, 'Sorting spikes on all thresholded channels...');
    set(hWait, 'position', get(hWait, 'position')+[0 150 0 0])
end
for nCh = 1:length(csChannels)
    sCh = csChannels{nCh}; % channel name
    if length(csChannels) > 1, waitbar(nCh/length(csChannels), hWait, ['Sorting spikes on channel ' sCh]); end
    % waveforms and spiketimes
    tSpikes = FV.tSpikes.(sCh);
    if isempty(tSpikes.waveforms); continue; end
    % Skip if there are less than 5 spikes
    if any(size(tSpikes.waveforms) < 5)
        uiwait(warndlg(['There are too few spikes on channel ' csChannels{nCh} ' to sort.']))
        tSpikes.hierarchy = struct('assigns', repmat(0, size(tSpikes.waveforms, 1), 1));
        FV.tSpikes.(sCh) = tSpikes;
        SetStruct(FV)
        continue;
    end
    tSpikes.Fs = tSpikes.Fs(1);
    % Over-clustering
    tSpikes = ss_kmeans(tSpikes);
    % Aggregation
    tSpikes = ss_energy(tSpikes);
    tSpikes = ss_aggregate(tSpikes);
    % Save results
    FV.tSpikes.(sCh) = tSpikes;
    SetStruct(FV)
    bSortGood = 1;
end
if length(csChannels) > 1, close(hWait); end
if bSortGood
    ShowOverlappingSpikeClusters;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sCh, bResult] = SelectChannelNumber(varargin)
% Select channel interactively from drop-down list.
%
% Usage:
%   SelectChannelNumber(CH, TIT, DEF)
%       CH   list of all channels (cell)
%       TIT  title of dialog
%       DEF  default channel (string)
%
if ~IsDataLoaded, return, end
sCh = []; bResult = 0;
[FV,hWin] = GetStruct();

% Strings representing all individual and stereo/tetrode electrodes
sTit = 'Select channel';
if isempty(varargin)
    csChannels = FV.csChannels;
else
    csChannels = varargin{1};
    if length(varargin) > 1, sTit = varargin{2}; end
    if length(varargin) == 3, sDefCh = varargin{3}; else sDefCh = ''; end
end
csChannelNames = csChannels;

% Get channel descriptions
csDescriptions = GetChannelDescription(csChannels);

% Append channel descriptions to channel DAQ names
for c = 1:length(csChannels)
    sDAQ = csChannels{c};
    sDescr = csDescriptions{c};
    if ~isempty(sDescr)
        csChannels{c} = [sDAQ '  (' sDescr ')'];
    end
end

% If a default selection exists, re-order csChannels so default appears
% first in list
if ~isempty(sDefCh)
    % Find match
    nMatch = [];
    for i = 1:length(csChannels)
        if length(csChannels{i}) >= length(sDefCh)
            if strcmp(csChannels{i}(1:length(sDefCh)), sDefCh)
                nMatch = i;
                break
            end
        end
    end
    % Re-order
    if ~isempty(nMatch)
        vNewOrder = [nMatch setdiff(1:length(csChannels), nMatch)];
        csChannels = csChannels(vNewOrder);
        csDescriptions = csDescriptions(vNewOrder);
        csChannelNames = csChannelNames(vNewOrder);
    end
end

bResult = 0;
if isempty(csChannels), return, end
global g_nCh
g_nCh = 1;
if length(csChannels) == 1
    % Return only channel if channel count == 1
    sCh = csChannelNames{1}; bResult = 1;
else
    global g_bClosePressed
    g_bClosePressed = 0;
    hFig = figure;
    set(hFig, 'position', [get(0,'PointerLocation')-50 360 25], 'menu', 'none', ...
        'Name', sTit, 'NumberTitle', 'off', 'visible', 'off', ...
        'CloseRequestFcn', 'delete(gcf); global g_bClosePressed; g_bClosePressed = 1;')
    if exist('centerfig') centerfig(hFig, hWin); end
    set(hFig, 'visible', 'on')
    uicontrol(hFig, 'Style', 'popupmenu', 'Position', [10 5 300 20], ...
        'String', csChannels, 'Callback', 'global g_nCh; g_nCh=get(gco, ''value'');', ...
        'value', 1);
    uicontrol(hFig, 'Style', 'pushbutton', 'Position', [315 5 40 20], 'Callback', 'delete(gcf)', 'String', 'OK'); % OK button
    uiwait % wait for user to close window or click OK button
    % if user pressed the close button, we dont return a channel
    if ~isempty(g_nCh) && ~g_bClosePressed
        sCh = csChannelNames{g_nCh};
        bResult = 1;
    else sCh = []; end
end
clear global g_nCh
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ViewWaveforms(varargin)
% View all waveforms (unsorted)
% 
% Usage:
%   ViewWaveforms()
%
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
nSubplots = 0; % get number of channels for which spikes exist
csFieldnames = fieldnames(FV.tSpikes);
for i = 1:length(csFieldnames)
    if ~isempty(FV.tSpikes.(csFieldnames{i})), nSubplots = nSubplots + 1; end;
end
if isempty(csFieldnames)
    hWarn = warndlg('No waveforms have been extracted on any channel. You must run spike detection before viewing waveforms.');
    uiwait(hWarn); return
end
hFig = figure;
ThemeObject(hFig)
set(hFig, 'name', 'Spiky Waveforms', 'menubar', 'none')
p = 1;
% Find how many channels should be shown
for i = 1:length(csFieldnames) % iterate over individual channels
    if isempty(FV.tSpikes.(csFieldnames{i})), continue, end
    if isempty(FV.tSpikes.(csFieldnames{i}).waveforms), continue, end
    sCh = csFieldnames{i};
    % Check if this is a multitrode
    if ~isempty(strfind(sCh,'__')) % is
        vIndx = strfind(sCh,'__');
        nFs = FV.tData.([sCh(1:vIndx(1)-1) '_KHz']) * 1000; % sampling frequency of 1st ch (Hz)
    else
        nFs = FV.tData.([sCh '_KHz']) * 1000; % sampling frequency (Hz)
    end
    subplot(2,nSubplots,p);
    
    nWf = size(FV.tSpikes.(csFieldnames{i}).waveforms, 1); % number of waveforms
    nMaxWf = 5000; % max number of waveforms to display
    iShow = randperm(nWf);
    if nWf > nMaxWf
        iShow = iShow(1:nMaxWf); % show random waveforms
    end
    
    vXX = repmat([([1:size(FV.tSpikes.(csFieldnames{i}).waveforms,2)]/nFs)*1000]', ...
        1, length(iShow));
    vYY = FV.tSpikes.(csFieldnames{i}).waveforms(iShow, :)';
    
    hLine = plot(vXX, vYY, ':', 'linewidth', .5);
    ThemeObject(hLine)
    
    % Show all threshold lines
    % TODO
    
    hTit = title([sCh ' - ' num2str(nWf) ' spikes'], 'interpreter', 'none');
    ThemeObject(hTit)
    ylabel('Amplitude')
    axis tight;
    ThemeObject(gca)
    set(gca, 'xlim', [0 max(vXX(:,1))], 'xtick', 0:.5:max(vXX(:,1)));
    subplot(2,nSubplots,p+nSubplots); p = p + 1;
    hist2d(FV.tSpikes.(csFieldnames{i}).waveforms);
    ThemeObject(gca)
    set(gca, 'xtick', []);
    ylabel('Amplitude (mV)')
    xlabel('ms')
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ViewWaveformPCs(varargin)
%
if ~IsDataLoaded, return, end
[FV,hWin] = GetStruct();
% Select channel from list of sorted units
csFieldnames = fieldnames(FV.tSpikes)';
vDel = [];
for fn = 1:length(csFieldnames)
    if ~isfield(FV.tSpikes.(csFieldnames{fn}), 'hierarchy'), vDel(end+1) = fn; end
end
csFieldnames(vDel) = [];
if ~CheckIfSorted, return, end
[sCh, bResult] = SelectChannelNumber(csFieldnames);
if ~bResult, return, end

% assign colors
vUnique = unique(FV.tSpikes.(sCh).hierarchy.assigns);
vIndx = find(vUnique > 0);
nIndx = find(vUnique == 0);

mColsTemp = ones(max(FV.tSpikes.(sCh).hierarchy.assigns), 3);
mColsTemp(vUnique(vIndx), :) = FV.mColors(vIndx, :);

FV.tSpikes.(sCh).overcluster.colors = mColsTemp;
ssg_databrowse3d(FV.tSpikes.(sCh), FV.tSpikes.(sCh).hierarchy.assigns, unique(FV.tSpikes.(sCh).hierarchy.assigns))
hFig = gcf;
ThemeObject(hFig)
set(hFig, 'Name', 'Spiky Principal Components', 'ToolBar', 'figure')
hAx = get(hFig, 'Children');
ThemeObject(hAx(strcmp(get(hAx, 'type'), 'axes')))

sDescr = GetChannelDescription(sCh);
hHeader = header([sDescr '  ' sCh], 10);
ThemeObject(hHeader)
set(hHeader, 'interpreter', 'none')

% set colors on legend
hLeg = findobj(hFig, 'type', 'axes', 'Tag', 'legend');
set(hLeg, 'color', 'none', 'textcolor', 'w', 'location', 'northeast')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function csDescr = GetChannelDescription(varargin)
% Get the descriptive string of a channel name.
% 
% Usage:
%   GetChannelDescription(S) where S is a string
%   GetChannelDescription(C) where C is a cell of strings
% 
% The function returns a cell or string as a function of C or S,
% respectively.
%

[FV, ~] = GetStruct();
if ischar(varargin{1})
    cCh = {varargin{1}};
else
    cCh = varargin{1};
end

csDescr = cell(1, length(cCh));
for iCs = 1:length(cCh)
    sCh = cCh{iCs};
    if isfield(FV, 'tChannelDescriptions')
        if isfield(FV.tChannelDescriptions, 'sChannel') && isfield(FV.tChannelDescriptions, 'sDescription')
            nIndx = find(strcmp({FV.tChannelDescriptions.sChannel}, sCh));
            if ~isempty(nIndx)
                csDescr{iCs} = FV.tChannelDescriptions(nIndx).sDescription;
            else
                csDescr{iCs} = '';
            end
        end
    end
end

% If a single string was passed, return a string
if ischar(varargin{1})
    csDescr = csDescr{1};
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function KeyboardMode(varargin)
% Enter keyboard mode inside this function.
%
[FV, hWin] = GetStruct();
clc
sp_disp(sprintf('%s\n%s\n%s\n\t%s\n\t%s', ...
    '** YOU ARE NOW IN SPIKY''S KEYBOARD MODE **', ...
    'The FV structure containing all spike data and settings has been loaded into', ...
    'this workspace. You may modify its contents directly. To save changes, run:', ...
    'SetStruct(FV)', ...
    'return' ) );
keyboard
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AddTetrode(varargin)
% Add a multi-trode (2, 3 or 4 wires)
% AddTetrode(sCh1__2) will use the channel names supplied
% AddTetrode() will ask user which channels to combine
%
if ~IsDataLoaded, return, end
global g_vElectrodes
g_vElectrodes = [1 1 1 1];
[FV, hWin] = GetStruct();

% Get electrode names
csFieldnames = fieldnames(FV.tSpikes)';
cFields = [{''} csFieldnames];
if ~isempty(varargin{2})
    % Get channel names from function input
    sNewName = varargin{2};
    % Split multitrode name into its constituents
    vStartIndx = [1 strfind(sNewName, '__')+2];
    vEndIndx = [vStartIndx(2:end)-3 length(sNewName)];
    cSelChannels = {};
    vIndx = [];
    for i = 1:length(vStartIndx)
        sCh = sNewName(vStartIndx(i):vEndIndx(i)); % channel name
        % get channel index
        nIndx = find(strcmp(csFieldnames, sNewName(vStartIndx(i):vEndIndx(i))));
        if isempty(nIndx) % check if channel exists
            hWarn = warndlg(['Channel ' sCh ' was attempted combined in a multitrode but this channel doest not exist in this datafile. Close this window to continue combining the remaining channels.']);
            uiwait(hWarn)
        else
            vIndx(end+1) = nIndx + 1;
            % assign channel
            cSelChannels{end+1} = sCh;
        end
    end
    % check that we remain with at least 2 channels
    if length(vIndx) > 2
        hWarn = warndlg('Too few channels are used to define the new multitrode. At least 2 are required.');
        uiwait(hWarn)
        return
    end
    % re-create sNewName since there's a chance it changed durin the input error check above
    sNewName = cSelChannels{1};
    for nCh = 2:length(cSelChannels)
        sNewName = [sNewName '__' cSelChannels{nCh}];
    end
else
    % Get channel names from user input
    if length(csFieldnames) < 2 % abort if spikes have been detected on < 2 ch
        warndlg(['You must first run spike detection on at least two channels.'])
        return
    end
    hFig = figure; % create window
    vPos = get(0,'PointerLocation');
    set(hFig, 'position', [vPos(1) vPos(2)-150 250 160], 'menu', 'none', 'Name', 'Define stereo/tetrode', 'NumberTitle', 'off')
    nPos = 125;
    for i = 1:4
        uicontrol(hFig, 'Style', 'text', 'Position', [10 nPos 80 20], 'String', ['Electrode ' num2str(i)]);
        uicontrol(hFig, 'Style', 'popupmenu', 'Position', [90 nPos 150 20], 'String', cFields, 'Callback', ['global g_vElectrodes; g_vElectrodes(' num2str(i) ') = get(gco, ''value'');']);
        nPos = nPos - 25;
    end
    uicontrol(hFig, 'Style', 'pushbutton', 'Position', [85 nPos-15 80 30], 'Callback', 'close(gcf)', 'String', 'OK' ); % OK button
    uiwait % wait for user to close window or click OK button

    % Check if selection is valid
    if length(find(g_vElectrodes>1)) < 2
        warndlg('You must select 2 or more electrodes')
        return
    end

    % Selected combination
    cSelChannels = cFields(g_vElectrodes(g_vElectrodes>1));
    sNewName = sprintf('%s__',cSelChannels{:});
    sNewName = sNewName(1:end-2);

    % Check that combination was not selected previously
    bComboExists = 0;
    if any(strcmp(fieldnames(FV.tSpikes), sNewName))
        bComboExists = 1;
        warndlg(['This particular electrode combination has already been selected.'])
        return
    end

    % Channel indices
    vIndx = g_vElectrodes(g_vElectrodes>1);
end

% Get number of spikes on each channel
vNumSpikes = [];
for i = vIndx
    vNumSpikes(end+1) = length(FV.tSpikes.(cFields{i}).spiketimes); % number of spikes
end
nMaxNumSpikes = max(vNumSpikes); % channels with largest number of spikes (used to condition data in next loop)

% Create matrix containing spiketimes of all channels
mSpiketimes = [];
mCont = [];
vFs = []; vThreshT = []; vThreshV = [];
for i = vIndx
    % sampling rate (KHz)
    vFs(end+1) = FV.tData.([cFields{i} '_KHz']);
    % spiketimes of this channel
    vSpiketimes = FV.tSpikes.(cFields{i}).spiketimes;
    % subtract absolute begin time of 1st sample
    nBeginTime = FV.tData.([cFields{i} '_TimeBegin']); % sec
    vSpiketimes = vSpiketimes - (nBeginTime*vFs(end)*1000);
    % join spiketimes
    mSpiketimes = [mSpiketimes [vSpiketimes; repmat(-1, nMaxNumSpikes-length(FV.tSpikes.(cFields{i}).spiketimes), 1)]]; % pad with -1 to maxlen
    % continuous data (used to concatenate waveforms when no spike was detected on a channel)
    mCont = [mCont FV.tData.(cFields{i})'];
    % Threshold crossing time
    vThreshT(end+1) = FV.tSpikes.(cFields{i}).threshT;
    % Threshold crossing voltage
    vThreshV(end+1,1:2) = FV.tSpikes.(cFields{i}).threshV;
end
% Check that sampling rate is the same across channels
if length(unique(vFs)) > 1
    warndlg(['Sampling rates on selected channels are different. Concatenation of spikes assumes the sampling rate is the same across all channels.'])
    return
end
% Concatenate spikes across selected channels
[mConcatSpikes, vCombSpiketimes] = concatenate_spikes(mSpiketimes, mCont, FV.nMaxCrossChannelSpikeJitter, FV.nSpikeTrigPreMs, FV.nSpikeTrigPostMs, vFs(1)*1000);
% Add back the begin time of 1st sample
vCombSpiketimes = vCombSpiketimes + (nBeginTime*vFs(end)*1000);
FV.tSpikes.(sNewName).waveforms = mConcatSpikes;
FV.tSpikes.(sNewName).spiketimes = vCombSpiketimes';
FV.tSpikes.(sNewName).Fs = vFs(1)*1000;
FV.tSpikes.(sNewName).threshV = vThreshV(1,1:2); % TODO: NOT VALID!
FV.tSpikes.(sNewName).threshT = vThreshT(1);
SetStruct(FV)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RemoveTetrode(varargin)
%
[FV, ~] = GetStruct();
% Find which are stereo/tetrodes
csFieldnames = fieldnames(FV.tSpikes);
vIndx = [];
for fn = 1:length(csFieldnames)
    if ~isempty(strfind(csFieldnames{fn}, '__')), vIndx(end+1) = fn; end
end
if isempty(vIndx)
    warndlg('No defined electrodes to remove!');
    return
end
csFieldnames = {'' csFieldnames{vIndx}};
hFig = figure;
set(hFig, 'position', [get(0,'PointerLocation')-50 360 25], 'menu', 'none', 'Name', 'Select channel', 'NumberTitle', 'off')
global g_nCh
g_nCh = 1;
uicontrol(hFig, 'Style', 'popupmenu', 'Position', [10 5 300 20], 'String', csFieldnames, 'Callback', 'global g_nCh; g_nCh=get(gco, ''value'');', 'value', 1);
uicontrol(hFig, 'Style', 'pushbutton', 'Position', [315 5 40 20], 'Callback', 'close(gcf)', 'String', 'OK' ); % OK button
uiwait % wait for user to close window or click OK button
if g_nCh == 1, return, end
FV.tSpikes = rmfield(FV.tSpikes, csFieldnames{g_nCh});
clear global g_nCh
SetStruct(FV)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FV = SetFVDefaults()
% Sets parameters in the FV structure to default values.
%
FV.sLoadedTrial = '';
FV.csDisplayChannels = {};      % channels to be displayed (analogue, continuous)
FV.csChannels = {};             % names of all channels in datafile
FV.csDigitalChannels = {};      % list of digital channels
FV.CurrentTheme = 'spiky';      % default theme
FV.mColors = distinguishable_colors(100, [.2 .2 .2]);
FV.tSpikeThresholdsNeg = struct([]);
FV.tSpikeThresholdsPos = struct([]);
FV.tGain = struct([]);
FV.tYlim = struct([]);
FV.tSpikes = struct([]);
FV.nSpikeTrigPreMs = 0.75; % ms
FV.nSpikeTrigPostMs = 1.75; % ms
FV.vXlim = [];
FV.nDeadTimeMs = 0.01; % ms
FV.bPlotRasters = 0;
FV.bPanOn = 0;
FV.nMaxCrossChannelSpikeJitter = 0.5; % maximal allowed spiketimmer jitter across electrodes (ms)
FV.tAmplitudeUnit = struct('sUnit', 'uV', 'nFactor', 10^6);
% Experiment variables
FV.tExperimentVariables = struct([]);
g_bBatchMode = false;
FV.sBatchLock = 'off';
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DistributeSettings(varargin)
% Save the current settings as settings for all other files in current directory.
% Overwrite settings of other files if they exist by default.
%
[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end
global g_hSpike_Viewer

% Get current file list
hFileList = findobj(g_hSpike_Viewer, 'Tag', 'MenuFileList');
csPaths = flipud(get(get(hFileList, 'Children'), 'UserData'));
if isempty(csPaths), return; end
if ischar(csPaths), return; end % single file loaded

sThisTrial = FV.sLoadedTrial;
FV = rmfield(FV, {'tSpikes' 'tData'}); % remove data fields
FV.tSpikes = struct([]);
FV.tData = struct([]);
FV.vXlim = []; % remove the x-limits since they might not apply to other files

% Iterate over all files in directory
SpikyWaitbar(0, length(csPaths));
for f = 1:length(csPaths)
    SpikyWaitbar(f, length(csPaths));
    drawnow

    % set sLoadedTrial
    sNewFile = csPaths{f};
    if strcmp(sNewFile, sThisTrial), continue, end % dont overwrite already loaded trial

    % load current data and settings (if they exist)
    sSPBFile = [sNewFile(1:end-4) '.spb'];
    if exist(sSPBFile, 'file') == 2
        tFileData = load(sSPBFile, 'FV', '-MAT');
        tFV_file = tFileData.FV;
    else
        tFV_file = [];
    end

    % Keep ALL fields existing in associated .spb file
    if ~isempty(tFV_file)
        FV.tData = tFV_file.tData;
    end
    FV.sLoadedTrial = sNewFile;
    % save settings (and imported fields)
    save([FV.sLoadedTrial(1:end-4) '.spb'], 'FV', '-v7.3')
end
sp_disp(sprintf('Done distributing current settings to %d files in current directory', length(csPaths)))
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RemoveSettings(varargin)
% Save the current settings as settings for all other files in current directory.
% Overwrite settings of other files if they exist by default.
%
[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end
global g_hSpike_Viewer

% Warn that this action is irreversible
switch questdlg(sprintf('DANGER ZONE:\nAll *.spb files associated with current files will be deleted from disk. This action is irreversible. Are you sure you want to continue?'), 'Spiky', 'Yes', 'No', 'Cancel', 'Cancel')
    case 'No', return
    case 'Cancel', return
end

% Get current file list
hFileList = findobj(g_hSpike_Viewer, 'Tag', 'MenuFileList');
csPaths = flipud(get(get(hFileList, 'Children'), 'UserData'));
if isempty(csPaths), return; end
if ischar(csPaths), return; end % single file loaded

% Iterate over all files in directory
SpikyWaitbar(0, length(csPaths));
for f = 1:length(csPaths)
    SpikyWaitbar(f, length(csPaths));
    delete([csPaths{f}(1:end-4) '.spb'])
end
sp_disp(sprintf('Done deleting settings on %d loaded files.', length(csPaths)))
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CreateMergeFile(varargin)
% Merge all open files into a common merge file.
% 
% User can select which types of data to merge.
%
global g_hSpike_Viewer g_vSelCh
[FV, ~] = GetStruct();

if ~IsDataLoaded(), return, end

% Ask user what data to merge
hFig = figure;
% TODO Remove Imported Channels? This appears to no longer be used.
cFields = {'Spikes' 'Events' 'Imported channels' 'Filtered channels' 'Displayed channels'};
vVals = [1 1 1 1 0];
nCL = length(cFields)*25+15;
vHW = get(g_hSpike_Viewer, 'position');
set(hFig, 'visible', 'off', 'position', [vHW(1:2) 250 nCL+25], 'menu', 'none', 'Name', 'Merge Files', 'NumberTitle', 'off')
drawnow
if exist('centerfig') centerfig(hFig, g_hSpike_Viewer); end
set(hFig, 'visible', 'on')
for i = 1:length(cFields)
    uicontrol(hFig, 'Style', 'checkbox', 'Position', [10 nCL 150 20], 'String', [' ' cFields{i}], ...
        'HorizontalAlignment', 'left', 'backgroundcolor', [.8 .8 .8], 'value', vVals(i));
    nCL = nCL - 25;
end
uicontrol(hFig, 'Style', 'pushbutton', 'Position', [75 nCL 50 20], ...
    'Callback', 'global g_vSelCh; g_vSelCh=flipud(get(get(gcf,''children''),''value'')); g_vSelCh=[g_vSelCh{:}];close(gcf)', ...
    'String', 'OK' ); % OK button

uiwait % wait for user to close window or click OK button
g_vSelCh = g_vSelCh(1:end-1);
cMergeData = cFields(logical(g_vSelCh)); % selected data
if isempty(cMergeData), return, end % no options were selected
if any(strcmpi(cMergeData, 'spikes')), bSpikes = 1;
else bSpikes = 0; end

if any(strcmpi(cMergeData, 'events')), bEvents = 1;
else bEvents = 0; end

if any(strcmpi(cMergeData, 'imported channels')), bImported = 1;
else bImported = 0; end

if any(strcmpi(cMergeData, 'filtered channels')), bFiltered = 1;
else bFiltered = 0; end

if any(strcmpi(cMergeData, 'displayed channels')), bDisplayed = 1;
else bDisplayed = 0; end

% Ask whether to save over existing file if it exists
sPath = [FV.sDirectory filesep 'MergeFile.spb'];
sAppendReplace = '';
if exist(sPath, 'file')
    sAppendReplace = questdlg('Merge file already exists on disk. Do you want to replace this file or write/append the new merge data to it?', ...
        'Merge file exists', 'Replace', 'Append', 'Append');
    if isempty(sAppendReplace), return, end
end

% Reset merge file contents
FV_merge.tSpikes = struct([]);
FV_merge.tData = struct([]);
FV_merge.tChannelDescriptions = struct([]);
FV_merge.csDisplayChannels = {};
FV_merge.csDigitalChannels = {};
tFiltCh = struct([]);      % temporary variable that holds filtered channels
tImported = struct([]); % temporary variable that holds imported data
tDisplayed = struct([]); % temporary variable that holds displayed channels
FV_merge.sDirectory = FV.sDirectory;
sDiffSampleRateAction = 'none';

if strcmpi(sAppendReplace, 'Append')
    % Get existing merge data
    % Note that we don't append to existing spike data, since spikes may
    % have been dejittered (in which case matrices wont merge readily).
    % Thus, whenever spikes are selected to be merged, they will always
    % replace ALL spikes in the merge structure, regardless.
    FV_old = load([FV_merge.sDirectory filesep 'MergeFile.spb'], '-MAT');
    if bImported || bFiltered || bEvents
        FV_merge.tData = FV_old.FV.tData;
        FV_merge.csDisplayChannels = FV_old.FV.csDisplayChannels;
        FV_merge.csDigitalChannels = FV_old.FV.csDigitalChannels;
    end
end

StartMergeMode(); % start merge mode

% Iterate over files in current directory
hFileList = findobj(g_hSpike_Viewer, 'Tag', 'MenuFileList');
csPaths = flipud(get(get(hFileList, 'Children'), 'UserData'));
if isempty(csPaths)
    warndlg('No additional files have been loaded into your workspace (see File->Files). You first need to load a directory, a directory tree or a custom list of files into your workspace before distributing settings and merging results.', 'Spiky');
    return
end
bResult = SpikyWaitbar(0, length(csPaths));
if iscell(csPaths); nLen = length(csPaths);
else nLen = 1; end
for f = 1:nLen
    bWaitResult = SpikyWaitbar(f, length(csPaths));
    if ~bWaitResult % waitbar was closed
        StopMergeMode();
        SpikyWaitbar(length(csPaths), length(csPaths));
        uiwait(warndlg('Merge mode has been cancelled.'))
        return
    end
    if iscell(csPaths)
        bResult = LoadTrial(csPaths{f}); % load trial data and settings
    else
        bResult = LoadTrial(csPaths);
    end
    drawnow
    [FV, ~] = GetStruct();
    
    % Get spike mergedata
    if bSpikes
        DetectSpikes(); % detect spikes on thresholded channels
        [tSpikes, ~] = GetStruct('tSpikes'); % get spikes from the updated FV
        FV.tSpikes = tSpikes;

        % concatenate multitrodes
        csMultiTrodes = {}; % get names of multitrodes
        csCh = fieldnames(FV.tSpikes);
        for nCh = 1:length(csCh)
            if ~isempty(strfind(csCh{nCh}, '__')) % multitrode found
                csMultiTrodes{end+1} = csCh{nCh};
            end
        end
        for nCh = 1:length(csMultiTrodes) % concatenate
            AddTetrode([], csMultiTrodes{nCh})
        end
        SaveResults(); % save results to this file's .spb file
        
        % Get spikes from all thresholded channels
        csChannels = fieldnames(FV.tSpikes);
        
        % Collect spiketime and associated information
        for nCh = 1:length(csChannels)
            % Skip channels that have no spikes
            if isempty(FV.tSpikes.(csChannels{nCh}).waveforms)
                continue
            end
            
            % Check that sampling rate is equal to other file
            nFs1 = 1; nFs2 = 1;
            if isfield(FV_merge, 'tData')
                % Check that channel exists in current FV file and merge file
                if isfield(FV_merge.tData, csChannels{nCh}) && isfield(FV.tData, csChannels{nCh})
                    % Compare sampling rates
                    if FV.tData.([csChannels{nCh} '_KHz']) ~= FV_merge.tData.([csChannels{nCh} '_KHz'])
                        nFs1 = FV_merge.tData.([csChannels{nCh} '_KHz']);
                        nFs2 = FV.tData.([csChannels{nCh} '_KHz']);
                        if strcmp(sDiffSampleRateAction, 'none') % re-use last selected action
                            % Warn that sampling rates are different
                            sMsg = sprintf('Sampling rate of channel %s is not equal across all files. You can choose to either ignore more spikes from this channel, or resample waveforms to the sampling rate of the first file (%d KHz). Your selection here will be repeated for all applicable channels.', csChannels{nCh}, FV_merge.tData.([csChannels{nCh} '_KHz']));
                            sAns = questdlg(sMsg, 'Spiky', 'Ignore Spikes', 'Resample Spikes', 'Cancel Merging', 'Resample Spikes');
                            switch sAns
                                case 'Ignore Spikes'
                                    sDiffSampleRateAction = 'ignore';
                                case 'Resample Spikes'
                                    sDiffSampleRateAction = 'resample';
                                otherwise
                                    return;
                            end
                        end
                    end
                end
            end

            % Get waveforms
            if nFs1 ~= nFs2 % resample or ignore if sample rates are different in this file
                switch sDiffSampleRateAction
                    case 'ignore'
                        continue;
                    case 'resample'
                        % Resample all spikes to samplerate of first file
                        mWaveforms = FV.tSpikes.(csChannels{nCh}).waveforms;
                        nNewSpikeLen = NaN;
                        if isfield(FV_merge.tSpikes.(csChannels{nCh}), 'waveforms')
                           nNewSpikeLen = size(FV_merge.tSpikes.(csChannels{nCh}).waveforms, 2);
                        end
                        if isnan(nNewSpikeLen)
                            nNewSpikeLen = round(size(mWaveforms, 2) * nFs1/nFs2);
                        end
                        mNewWaveforms = zeros(size(mWaveforms, 1), nNewSpikeLen);
                        for wf = 1:size(mWaveforms, 1)
                            vWaveform = mWaveforms(wf, :);
                            vNewWaveform = interp1(1:length(vWaveform), vWaveform, ...
                                linspace(1, length(vWaveform), nNewSpikeLen), 'spline');
                            mNewWaveforms(wf, :) = vNewWaveform;
                        end
                        
                        % Make sure all spikes cross threshold at time threshT
                        threshT = FV.tSpikes.(csChannels{nCh}).threshT; % samples
                        threshT = round(threshT * nFs1/nFs2); % samples
                        nThreshV = FV.tSpikes.(csChannels{nCh}).threshV;
                        vIndx = find(mNewWaveforms(:,threshT) >= nThreshV(1));                        
                        vRandShift = rand(1,length(vIndx)) ./ 10;
                        mNewWaveforms(vIndx,threshT) = (repmat(nThreshV(1), length(vIndx), 1) - abs((nThreshV(1).*vRandShift')))';
                        
                        mWaveforms = mNewWaveforms;
                        Fs = nFs1 * 1000; % Hz
                        
                        % Spiketimes must be change with new sample rate in mind
                        vSpiketimes = FV.tSpikes.(csChannels{nCh}).spiketimes; % samples
                        vSpiketimes = round(vSpiketimes .* (nFs1/nFs2)); % samples
                        
                        threshT = FV.tSpikes.(csChannels{nCh}).threshT; % samples
                        threshT = round(threshT * nFs1/nFs2); % samples
                end
            else
                mWaveforms = FV.tSpikes.(csChannels{nCh}).waveforms;
                Fs = FV.tSpikes.(csChannels{nCh}).Fs; % Hz
                vSpiketimes = FV.tSpikes.(csChannels{nCh}).spiketimes; % samples
                threshT = FV.tSpikes.(csChannels{nCh}).threshT;
            end

            FV_merge.csDisplayChannels = unique([csChannels{nCh} FV_merge.csDisplayChannels]);
            FV_merge.tData(1).(csChannels{nCh}) = [];
            FV_merge.tData(1).([csChannels{nCh} '_TimeBegin']) = 0;
            threshV = FV.tSpikes.(csChannels{nCh}).threshV;
            
            % Assign waveforms and spiketimes
            if isfield(FV_merge.tSpikes, csChannels{nCh}) % channel already exists
                % If new waveforms have different dimensions, check Fs:
                %   if Fs is same, then resample new waveforms to same length
                if size(mWaveforms, 2) ~= size(FV_merge.tSpikes.(csChannels{nCh}).waveforms, 2) && ~isempty(mWaveforms)
                    nNewLen = size(FV_merge.tSpikes.(csChannels{nCh}).waveforms, 2);
                    mWaveforms2 = interp1( 1:size(mWaveforms,2), ... % x
                        mWaveforms', ... % y
                        linspace(1, size(mWaveforms,2), nNewLen) ); % method
                    mWaveforms = mWaveforms2';
                end
                FV_merge.tSpikes.(csChannels{nCh}).waveforms = ... % waveforms
                    [FV_merge.tSpikes.(csChannels{nCh}).waveforms; mWaveforms];
                FV_merge.tSpikes.(csChannels{nCh}).spiketimes = ... % spiketimes
                    [FV_merge.tSpikes.(csChannels{nCh}).spiketimes; vSpiketimes];
                FV_merge.tSpikes.(csChannels{nCh}).threshV = ... % threshold
                    [FV_merge.tSpikes.(csChannels{nCh}).threshV; threshV];
                FV_merge.tSpikes.(csChannels{nCh}).threshT = ... % threshold
                    [FV_merge.tSpikes.(csChannels{nCh}).threshT; threshT];
                FV_merge.tSpikes.(csChannels{nCh}).Fs = ... % samling frequency
                    [FV_merge.tSpikes.(csChannels{nCh}).Fs; Fs];
            else % channel doesnt already exists
                FV_merge.tSpikes(1).(csChannels{nCh})(1).waveforms = mWaveforms;
                FV_merge.tSpikes.(csChannels{nCh}).spiketimes = vSpiketimes;
                FV_merge.tSpikes.(csChannels{nCh}).threshV = threshV;
                FV_merge.tSpikes.(csChannels{nCh}).threshT = threshT;
                FV_merge.tSpikes.(csChannels{nCh}).Fs = Fs;
                % assign sampling frequency
                vIndx = strfind(csChannels{nCh},'__');
                if ~isempty(vIndx) % is multitrode
                    FV_merge.tData(1).([csChannels{nCh}(1:vIndx(1)-1) '_KHz']) = FV.tData.([csChannels{nCh}(1:vIndx(1)-1) '_KHz']);
                else % is unitrode
                    FV_merge.tData(1).([csChannels{nCh} '_KHz']) = FV.tData.([csChannels{nCh} '_KHz']);
                end
            end
        end
    end % end of spike merging

    bWaitResult = SpikyWaitbar(f, length(csPaths));
    
    % Merge digital lines
    if bEvents
        FV_merge.csDigitalChannels = unique([FV_merge.csDigitalChannels FV.csDigitalChannels]);
        for c = 1:length(FV.csDigitalChannels)
            vUpTimes = FV.tData.([FV.csDigitalChannels{c} '_Up']);
            vDownTimes = FV.tData.([FV.csDigitalChannels{c} '_Down']);
            
            % Remove downtimes that occur before the 1st uptime
            vDownTimes(vDownTimes <= min(vUpTimes)) = [];
            
            % Remove uptimes that occur after the last downtime
            vUpTimes(vUpTimes >= max(vDownTimes)) = [];
            
            if isfield(FV_merge.tData, [FV.csDigitalChannels{c} '_Up'])
                FV_merge.tData.([FV.csDigitalChannels{c} '_Up']) = ...
                    [FV_merge.tData.([FV.csDigitalChannels{c} '_Up']) vUpTimes];
                FV_merge.tData.([FV.csDigitalChannels{c} '_Down']) = ...
                    [FV_merge.tData.([FV.csDigitalChannels{c} '_Down']) vDownTimes];
            else
                FV_merge.tData(1).([FV.csDigitalChannels{c} '_Up']) = vUpTimes;
                FV_merge.tData.([FV.csDigitalChannels{c} '_Down']) = vDownTimes;
            end
            FV_merge.tData.([FV.csDigitalChannels{c} '_KHz']) = FV.tData.([FV.csDigitalChannels{c} '_KHz']);
            FV_merge.tData.([FV.csDigitalChannels{c} '_TimeBegin']) = FV.tData.([FV.csDigitalChannels{c} '_TimeBegin']);
            FV_merge.tData.([FV.csDigitalChannels{c} '_TimeEnd']) = FV.tData.([FV.csDigitalChannels{c} '_TimeEnd']);
        end
    end % end of events merging
    
    % Add fields that indicate start and end of file
    if ~isfield(FV_merge.tData, 'FileStart'), FV_merge.tData(1).FileStart = {}; end
    if ~isfield(FV_merge.tData, 'FileEnd'), FV_merge.tData(1).FileEnd = {}; end
    FV_merge.tData(1).FileStart(end+1).Timestamp = FV.tData.([FV.csDigitalChannels{c} '_TimeBegin']); % file start marker
    FV_merge.tData(1).FileStart(end).File = FV.sLoadedTrial;
    FV_merge.tData(1).FileEnd(end+1).Timestamp = FV.tData.([FV.csDigitalChannels{c} '_TimeEnd']); % file end marker
    FV_merge.tData(1).FileEnd(end).File = FV.sLoadedTrial;
    
    % Merge imported fields
    if bImported
        csFields = fieldnames(FV.tData);
        csChannels = {};
        for nCh = 1:length(csFields) % find imported fields
            if ~isempty(strfind(csFields{nCh}, '_Imported')), csChannels{end+1} = csFields{nCh}(1:end-9); end
        end
        csChannels = unique(csChannels);
        for nCh = 1:length(csChannels) % iterate over imported channels
            vCont = ChannelCalculator(FV.tData.(csChannels{nCh}), csChannels{nCh});
            [vCont, FV] = AdjustChannelGain(FV, vCont, csChannels{nCh});
            
            % Time vector
            nBegin = FV.tData.([csChannels{nCh} '_TimeBegin']); % sampling start, sec
            nEnd = FV.tData.([csChannels{nCh} '_TimeEnd']); % sampling end, sec
            nFs = FV.tData.([csChannels{nCh} '_KHz']) * 1000; % sampling frequency Hz
            vTime = (nBegin+1/nFs):(1/nFs):(nBegin+length(vCont)/nFs); % absolute time, sec

            % Apply channel filter
            [vCont, vTimeOut, nFsOut] = GetFilteredChannel(csChannels{nCh}, vCont, 'decimate');
            if ~isempty(vTimeOut), vTime = vTimeOut; end
            if ~isempty(nFsOut), nFs = nFsOut; end
            
            % Copy data to temporary variable
            if ~isfield(tImported, csChannels{nCh}), nIndx = 1;
            else nIndx = length(tImported.(csChannels{nCh})) + 1; end
            tImported(1).(csChannels{nCh})(nIndx).cont = vCont;
            tImported.(csChannels{nCh})(nIndx).timebegin = nBegin; % store begin time
            tImported.(csChannels{nCh})(nIndx).timeend = nEnd; % store end time
            tImported.(csChannels{nCh})(nIndx).time = vTime; % store continuous time vector, sec
            tImported.(csChannels{nCh})(nIndx).fs = nFs; % sampling rate
        end
        clear vTime vCont
    end
    
    % Merge filtered signals (EMG, LFP, ECoG etc)
    % Note: Any signal that has been Filtered through the GUI will be exported to
    % the Merge File (store in temporary var first and join with FV_merge after outer loop)
    if bFiltered
        if isfield(FV, 'tFilteredChannels')
            csChannels = {FV.tFilteredChannels.sChannel};
        else
            csChannels = {};
        end
        for nCh = 1:length(csChannels)
            % Check if channel has already been merged
            if isfield(tImported, csChannels(nCh)), continue; end
            
            % Check first that channel exists
            if ~isfield(FV.tData, csChannels{nCh}) continue, end
            
            % Raw signal
            sCh = csChannels{nCh};
            vCont = ChannelCalculator(FV.tData.(csChannels{nCh}), csChannels{nCh});
            [vCont, FV] = AdjustChannelGain(FV, vCont, sCh);
            
            % Time vector
            nBegin = FV.tData.([csChannels{nCh} '_TimeBegin']); % sampling start, sec
            nEnd = FV.tData.([csChannels{nCh} '_TimeEnd']); % sampling end, sec
            nFs = FV.tData.([csChannels{nCh} '_KHz']) * 1000; % sampling frequency Hz
            vTime = (nBegin+1/nFs):(1/nFs):(nBegin+length(vCont)/nFs); % absolute time, sec
            
            % Apply channel filter
            [vCont, vTimeOut, nFsOut] = GetFilteredChannel(csChannels{nCh}, vCont, 'decimate');
            if ~isempty(vTimeOut), vTime = vTimeOut; end
            if ~isempty(nFsOut), nFs = nFsOut; end
            
            % Copy data to temporary variable
            if ~isfield(tFiltCh, csChannels{nCh}), nIndx = 1;
            else nIndx = length(tFiltCh.(csChannels{nCh})) + 1; end
            tFiltCh(1).(csChannels{nCh})(nIndx).cont = vCont;
            tFiltCh.(csChannels{nCh})(nIndx).timebegin = nBegin; % store begin time
            tFiltCh.(csChannels{nCh})(nIndx).timeend = nEnd; % store end time
            tFiltCh.(csChannels{nCh})(nIndx).time = vTime; % store continuous time vector, sec
            tFiltCh.(csChannels{nCh})(nIndx).fs = nFs; % new sampling rate
        end
        clear vTime vCont
    end
    
    % Merge displayed fields
    % Displayed fields are merged at the original sampling rate, unless
    % filter have been configured for these channels. A displayed field
    % will not be merged again if it was already merged due it being an
    % Imported or Filtered channel.
    if bDisplayed
        % Get displayed channels
        csChannels = unique(FV.csDisplayChannels);
        
        % Verify that all displayed channels exist
        for cs = 1:length(csChannels)
            if ~isfield(FV.tData, csChannels(cs))
                csChannels(cs) = [];
            end
        end
        
        % Iterate over displayed channels
        for nCh = 1:length(csChannels)
            % Check if channel has already been merged
            if isfield(tFiltCh, csChannels(nCh)), continue; end
            if isfield(tImported, csChannels(nCh)), continue; end
            
            % Apply custom math
            vCont = ChannelCalculator(FV.tData.(csChannels{nCh}), csChannels{nCh});
            
            % Adjust gain
            [vCont, FV] = AdjustChannelGain(FV, vCont, csChannels{nCh});
            
            % Time vector
            nBegin = FV.tData.([csChannels{nCh} '_TimeBegin']); % sampling start, sec
            nEnd = FV.tData.([csChannels{nCh} '_TimeEnd']); % sampling end, sec
            nFs = FV.tData.([csChannels{nCh} '_KHz']) * 1000; % sampling frequency Hz
            vTime = (nBegin+1/nFs):(1/nFs):(nBegin+length(vCont)/nFs); % absolute time, sec

            % Apply channel filter
            [vCont, vTimeOut, nFsOut] = GetFilteredChannel(csChannels{nCh}, vCont, 'decimate');
            if ~isempty(vTimeOut), vTime = vTimeOut; end
            if ~isempty(nFsOut), nFs = nFsOut; end
            
            % Copy data to temporary variable
            if ~isfield(tDisplayed, csChannels{nCh}), nIndx = 1;
            else nIndx = length(tDisplayed.(csChannels{nCh})) + 1; end
            
            tDisplayed(1).(csChannels{nCh})(nIndx).cont = vCont;
            tDisplayed.(csChannels{nCh})(nIndx).timebegin = nBegin; % store begin time
            tDisplayed.(csChannels{nCh})(nIndx).timeend = nEnd; % store end time
            tDisplayed.(csChannels{nCh})(nIndx).time = vTime; % store continuous time vector, sec
            tDisplayed.(csChannels{nCh})(nIndx).fs = nFs; % sampling rate
        end
        clear vTime vCont
    end
    
    % Keep ALL channel descriptions (regardless of whether they were included in the merge file or not)
    % As first step, check for and ignore duplicates
    for i = 1:length(FV.tChannelDescriptions)
        sCh = FV.tChannelDescriptions(i).sChannel;
        if isempty(FV_merge.tChannelDescriptions)
            FV_merge.tChannelDescriptions = FV.tChannelDescriptions;
        else
            j = strcmp({FV_merge.tChannelDescriptions.sChannel}, sCh);
            if ~any(j)
                FV_merge.tChannelDescriptions(end+1).sChannel = FV.tChannelDescriptions(i).sChannel;
                FV_merge.tChannelDescriptions(end).sDescription = FV.tChannelDescriptions(i).sDescription;
            end
        end
    end
end

% Move imported, filtered and displayed channels into FV_merge. For each category, a
% continuous vector is created containing the concatenated signal from all files.
% Missing data segments are filled with NaNs
tAll.tImported = tImported;
tAll.tFiltCh = tFiltCh;
tAll.tDisplayed = tDisplayed;
csOrder = fieldnames(tAll);
for c = find([bImported bFiltered bDisplayed])
    cFields = fieldnames(tAll.(csOrder{c}));
    for nF = 1:length(cFields)
        
        % Initialize continuous vector (min begintime to max endtime)
        vTimeBegin = [tAll.(csOrder{c}).(cFields{nF}).timebegin]; % file endtimes, sec
        vTimeEnd = [tAll.(csOrder{c}).(cFields{nF}).timeend]; % file endtimes, sec
        vFs = [tAll.(csOrder{c}).(cFields{nF}).fs];
        if length(unique(vFs)) > 1
            uiwait(warndlg(sprintf('Failed importing channel %s. Sampling rate must be same across all files.', cFields{nF})))
            continue
        end
        
        nD = 1/vFs(1); % duration of one sampling points
        vContAll = [];
        nAbsBeginTime = round(min(vTimeBegin) / nD - 1); % begin time of first trial, samples
        vIndxAll = [];
        
        for i = 1:length(vTimeEnd)
            nBegin = tAll.(csOrder{c}).(cFields{nF})(i).timebegin; % sec
            nBegin = round(nBegin / nD); % samples
            vCont = tAll.(csOrder{c}).(cFields{nF})(i).cont;
            vIndx = (nBegin:[nBegin+length(vCont)-1])-nAbsBeginTime;
            vIndxAll = [vIndxAll vIndx];
            vContAll(vIndx) = vCont; % samples
        end
        
        % Replace empty segments with NaNs
        vContAll(setdiff(1:length(vContAll), vIndxAll)) = NaN;
        
        % store continuous vector in FV_merge
        FV_merge.tData(1).(cFields{nF}) = vContAll;
        FV_merge.tData.([cFields{nF} '_KHz']) = vFs(1)/1000;
        FV_merge.tData.([cFields{nF} '_TimeBegin']) = min(vTimeBegin);
        FV_merge.tData.([cFields{nF} '_TimeEnd']) = max(vTimeEnd);
    end
end

% 
% % Move imported data into FV_merge after creating one continuous vector
% % containing the concatenated signal from all files. Missing data segments
% % are filled with NaNs
% if bImported
%     %cFields = fieldnames(tImported);
%     for nF = 1:length(cFields)
%         % initialize continuous vector (min begintime to max endtime)
%         vTimeBegin = [tImported.(cFields{nF}).timebegin]; % file endtimes, sec
%         vTimeEnd = [tImported.(cFields{nF}).timeend]; % file endtimes, sec
%         vFs = [tImported.(cFields{nF}).fs];
%         if length(unique(vFs)) > 1
%             uiwait(warndlg('Sorry, cannot export all imported channels. Sampling rate must be same across all files. Use of different sampling rates for the same channel is currently not supported. Will continue to try to export remaining imported channels.'))
%             continue
%         end
%         nD = 1/vFs(1); % duration of one sampling points
%         vContAll = [];
%         nAbsBeginTime = round(min(vTimeBegin) / nD - 1); % begin time of first trial, samples
%         vIndxAll = [];
%         for i = 1:length(vTimeEnd)
%             nBegin = tImported.(cFields{nF})(i).timebegin; % sec
%             nBegin = round(nBegin / nD); % samples
%             vCont = tImported.(cFields{nF})(i).cont;
%             vIndx = (nBegin:[nBegin+length(vCont)-1])-nAbsBeginTime;
%             vIndxAll = [vIndxAll vIndx];
%             vContAll(vIndx) = vCont; % samples
%         end
%         % Replace empty segments with NaNs
%         vContAll(setdiff(1:length(vContAll), vIndxAll)) = NaN;
%         
%         % store continuous vector in FV_merge
%         FV_merge.tData(1).(cFields{nF}) = vContAll;
%         FV_merge.tData.([cFields{nF} '_KHz']) = vFs(1)/1000;
%         FV_merge.tData.([cFields{nF} '_TimeBegin']) = min(vTimeBegin);
%         FV_merge.tData.([cFields{nF} '_TimeEnd']) = max(vTimeEnd);
%     end
% end
% 
% % Move FILTERED data into FV_merge after creating one continuous vector
% % containing the concatenated signal from all files. Missing data segments
% % are filled with NaNs
% if bFiltered
%     %cFields = fieldnames(tFiltCh);
%     for nF = 1:length(cFields)
%         % initialize continuous vector (min begintime to max endtime)
%         %vTimeBegin = [tFiltCh.(cFields{nF}).timebegin]; % file endtimes, sec
%         %vTimeEnd = [tFiltCh.(cFields{nF}).timeend]; % file endtimes, sec
%         %vFs = [tFiltCh.(cFields{nF}).fs];
%         %if length(unique(vFs)) > 1
%         %    uiwait(warndlg('Sorry, cannot export all imported channels. Sampling rate must be same across all files. Use of different sampling rates for the same channel is currently not supported. Will continue to try to export remaining imported channels.'))
%         %    continue
%         %end
%         %nD = 1/vFs(1); % duration of one sampling points
%         %vContAll = [];
%         %nAbsBeginTime = round(min(vTimeBegin) / nD - 1);
%         %vIndxAll = [];
%         %for i = 1:length(vTimeEnd)
%             %nBegin = tFiltCh.(cFields{nF})(i).timebegin; % sec
%             %nBegin = round(nBegin / nD); % samples
%             %vCont = tFiltCh.(cFields{nF})(i).cont;
%             %vIndx = (nBegin:[nBegin+length(vCont)-1])-nAbsBeginTime;
%             %vIndxAll = [vIndxAll vIndx];
%             %vContAll(vIndx) = vCont; % samples
%         %end
%         % Replace empty segments with NaNs
%         %vContAll(setdiff(1:length(vContAll), vIndxAll)) = NaN;
%         
%         % store continuous vector in FV_merge
%         %FV_merge.tData(1).(cFields{nF}) = vContAll;
%         %FV_merge.tData.([cFields{nF} '_KHz']) = vFs(1)/1000;
%         %FV_merge.tData.([cFields{nF} '_TimeBegin']) = min(vTimeBegin);
%         %FV_merge.tData.([cFields{nF} '_TimeEnd']) = max(vTimeEnd);
%     end
% end

% Set FV as FV_merge
if strcmpi(sAppendReplace, 'Replace')
    % In Replace mode, reset FV to defaults
    FV = SetFVDefaults();
elseif strcmpi(sAppendReplace, 'Append')
    % In Append mode, load and use as basis existing mergefile
    FV_old = load([FV_merge.sDirectory filesep 'MergeFile.spb'], '-MAT');
    FV = FV_old.FV;
end
if bSpikes, FV.tSpikes = FV_merge.tSpikes; end
FV.tData = FV_merge.tData;
FV.sDirectory = FV_merge.sDirectory;
FV.csDisplayChannels = FV_merge.csDisplayChannels;
FV.csDigitalChannels = FV_merge.csDigitalChannels;
FV.tChannelDescriptions = FV_merge.tChannelDescriptions;
FV.bPlotRasters = 1;
FV.sLoadedTrial = [FV_merge.sDirectory filesep 'MergeFile.spb'];
SetStruct(FV)
save(FV.sLoadedTrial, 'FV', '-v7.3') % Save merged results to disk
set(g_hSpike_Viewer, 'Name', sprintf('%s - Spiky', GetCurrentFile()))

ViewTrialData();

% Check merge file for potential problems
%
% TODO
%   If start times repeat, allow user to force trials to be sequential?
vStart = FV.tData.DAQ_Start_Up; % start time duplicates
if length(unique(vStart)) < length(vStart)
    uiwait(warndlg(sprintf('There is a potential problem with the merged files:\nThere are duplicates of file Start times')))
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StartMergeMode(varargin)
% Initialize workspace for being in MERGE MODE
% print message that we are in a MERGE MODE.
%
[FV, hWin] = GetStruct();
figure(hWin)
global g_bMergeMode
g_bMergeMode = 1;

% Update menu items
set(findobj(hWin, 'Label', 'E&xit Merge Mode'), 'Enable', 'on')
set(findobj(hWin, 'Label', 'Set &Threshold...'), 'Enable', 'off')
set(findobj(hWin, 'Label', '&Auto-Detect Thresholds...'), 'Enable', 'off')
set(findobj(hWin, 'Label', '&Run Spike Detection'), 'Enable', 'off')
set(findobj(hWin, 'Label', 'Add &Multitrode...'), 'Enable', 'off')
set(findobj(hWin, 'Label', 'Remove Multitrode...'), 'Enable', 'off')
set(findobj(hWin, 'Label', 'Set Max &Jitter...'), 'Enable', 'off')
set(findobj(hWin, 'Label', 'Set Pre-Trigger...'), 'Enable', 'off')
set(findobj(hWin, 'Label', 'Set Post-Trigger...'), 'Enable', 'off')
%set(findobj(hWin, 'Label', 'Set Spike Deadtime...'), 'Enable', 'off')
set(findobj(hWin, 'Label', '&Create Merge File'), 'Enable', 'off')
set(findobj(hWin, 'Label', '&Distribute Settings'), 'Enable', 'off')
set(findobj(hWin, 'Label', 'Show P&ower Spectral Densities'), 'Enable', 'off')
set(findobj(hWin, 'Label', 'Digitize Channel'), 'Enable', 'off')
set(findobj(hWin, 'Label', 'Digitize Manually'), 'Enable', 'off')
hRaster = findobj(hWin, 'Label', 'Plot &Rasters');
set(hRaster, 'Enable', 'off', 'Checked', 'on')
FV.bPlotRasters = 1;
SetStruct(FV)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StopMergeMode(varargin)
% Clean up and exit from merge mode.
%
[~, hWin] = GetStruct();
global g_bMergeMode
g_bMergeMode = 0;
% Update menu items
set(findobj(hWin, 'Label', 'Set &Threshold...'), 'Enable', 'on')
set(findobj(hWin, 'Label', '&Auto-Detect Thresholds...'), 'Enable', 'on')
set(findobj(hWin, 'Label', '&Run Spike Detection'), 'Enable', 'on')
set(findobj(hWin, 'Label', 'Define &Multitrode...'), 'Enable', 'on')
set(findobj(hWin, 'Label', 'Set Max &Jitter...'), 'Enable', 'on')
set(findobj(hWin, 'Label', 'Set Pre-Trigger...'), 'Enable', 'on')
set(findobj(hWin, 'Label', 'Set Post-Trigger...'), 'Enable', 'on')
%set(findobj(hWin, 'Label', 'Set Spike Deadtime...'), 'Enable', 'on')
set(findobj(hWin, 'Label', '&Create Merge File'), 'Enable', 'on')
set(findobj(hWin, 'Label', '&Distribute Settings'), 'Enable', 'on')
set(findobj(hWin, 'Label', 'Show P&ower Spectral Densities'), 'Enable', 'on')
set(findobj(hWin, 'Label', 'Digitize Channel'), 'Enable', 'on')
set(findobj(hWin, 'Label', 'Digitize Manually'), 'Enable', 'on')
set(findobj(hWin, 'Label', 'Plot &Rasters'), 'Enable', 'on')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadMergeFile(varargin)
% Load a file containing data of already merged files.
%
[FV,hWin] = GetStruct();

% Check if file path is given in function input
if ~isempty(varargin{1}) && ischar(varargin{1}) && length(varargin) == 2
    sPath = varargin{1};
    sFile = varargin{2};
else
    % Select merge file
    [sFile, sPath] = uigetfile( {'*.spb','Spiky files (*.spb)'}, 'Select merged file');
    if sFile == 0, return, end
end
FV.sDirectory = sPath;
% Load merge file
StartMergeMode();
sFilePath = [sPath sFile];
if exist(sFilePath, 'file') > 0
    load(sFilePath, 'FV', '-MAT');
    FV.sDirectory = sPath;
    FV.sLoadedTrial = sFile;
else
    uiwait(warndlg('No merge file exists in the current directory.'))
    StopMergeMode();
    return
end
SetStruct(FV)
ViewTrialData()
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FigureExport(varargin)
% Export the Spiky window to a graphics file
global g_hSpike_Viewer % handle to Spiky window
% Select output path and format
[sFile, sPath] = uiputfile( ...
    {'*.fig' , 'MATLAB Figure (*.fig)'; ...
    '*.ai','Adobe Illustrator File (*.ai)'; ...
    '*.bmp','Bitmap file (*.bmp)'; ...
    '*.eps','EPS file (*.eps)'; ...
    '*.emf','Enhanced metafile (*.emf)'; ...
    '*.jpg','JPEG image (*.jpg)'; ...
    '*.pdf','Portable Document Format (*.pdf)'; ...
    '*.png','Portable Network Graphics File (*.png)'; ...
    '*.tif','TIFF image (*.tif)'}, ...
    'Export Spiky window as');
if all(sFile == 0), return, end % user cancelled save dialog
saveas(g_hSpike_Viewer, [sPath sFile])
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InvertColors(varargin)
% Toggle the InvertHardcopy parameter of the Spiky window. This feature
% decides whether the background color of the window should be included in
% exports.
global g_hSpike_Viewer % handle to Spiky window
hMenuItem = findobj(g_hSpike_Viewer, 'Label', '&Invert Colors'); % handle of 'Invert Color' option in meny
switch get(g_hSpike_Viewer, 'InvertHardcopy')
    case 'on' % turn off if on
        set(g_hSpike_Viewer, 'InvertHardcopy', 'off')
        set(hMenuItem, 'Checked', 'off');
    case 'off' % turn on if off
        set(g_hSpike_Viewer, 'InvertHardcopy', 'on')
        set(hMenuItem, 'Checked', 'on');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChannelStatistics(varargin)
% Show basic channel statistics, such as number of detected spikes
% 
% Format:
%   Channel_name  Fs  #spikes
%
[FV,hWin] = GetStruct();
% iterate over thresholded channels
csChannels = fieldnames(FV.tSpikes);
cAboutText = {};
for nCh = 1:length(csChannels)
    cAboutText{end+1} = [csChannels{nCh} '      ' ...
        num2str(size(FV.tSpikes.(csChannels{nCh}).waveforms,1)) ' spikes' ...
        '   ' num2str(FV.tSpikes.(csChannels{nCh}).Fs(1)/1000) 'kHz'];
end
if isempty(cAboutText)
    cAboutText = {'No spiking statistics available.'};
end
hBox = msgbox(cAboutText, 'Channel Statistics');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomAmplitude(varargin)
% Click on an axes and zoom vertically
%
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
global g_hSpike_Viewer % handle to Spiky window

sPointer = get(g_hSpike_Viewer, 'Pointer');
set(g_hSpike_Viewer, 'Pointer', 'crosshair')
zoom off
set(g_hSpike_Viewer, 'WindowButtonDownFcn', 'set(gcbo, ''WindowButtonDownFcn'', [])')
waitfor(g_hSpike_Viewer, 'WindowButtonDownFcn') % wait for user to click in figure
hAx = gca;
ThemeObject(hAx)
zoom yon % zoom axis
try % wait for user to click in figure
    waitfor(hAx, 'ylim')
catch
    return
end
zoom off
set(g_hSpike_Viewer,'Pointer','arrow')
% store y-limits
if ishandle(hAx)
    sCh = get(hAx, 'Tag');
    FV.tYlim.(sCh) = get(hAx, 'ylim');
    ThemeObject(hAx)
end
SetStruct(FV)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotRasters(varargin)
% 
%
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();

% Update menu item
global g_hSpike_Viewer % handle to Spiky window
hMenuItem = findobj(g_hSpike_Viewer, 'Label', 'Plot &Rasters'); % handle of 'Plot Rasters'
switch get(hMenuItem, 'Checked')
    case 'on' % turn off if on
        set(hMenuItem, 'Checked', 'off'); FV.bPlotRasters = 0;
    case 'off' % turn on if off
        if ~CheckIfSorted
            set(hMenuItem, 'Checked', 'off'); FV.bPlotRasters = 0;
        else set(hMenuItem, 'Checked', 'on'); FV.bPlotRasters = 1; end
end
SetStruct(FV)
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotInstSpikeRate(varargin)
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
% update menu item
global g_hSpike_Viewer % handle to Spiky window
hMenuItem = findobj(g_hSpike_Viewer, 'Label', 'Plot &Instantaneous Rate');
switch get(hMenuItem, 'Checked')
    case 'on' % turn off if on
        set(hMenuItem, 'Checked', 'off'); FV.bPlotInstSpikeRate = 0;
    case 'off' % turn on if off
        if ~CheckIfSorted
            set(hMenuItem, 'Checked', 'off'); FV.bPlotInstSpikeRate = 0;
        else set(hMenuItem, 'Checked', 'on'); FV.bPlotInstSpikeRate = 1; end
end
SetStruct(FV)
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExitSpiky(varargin)
% Quit and exit Spiky
[~, hWin] = GetStruct();

% Check whether there are changes to be saved
if ~isempty(strfind(get(hWin, 'Name'), '*'))
    switch questdlg('Save changes to current file?', 'Spiky', 'Yes', 'No', 'Cancel', 'Yes')
        case 'Yes', SaveResults;
        case 'Cancel', return
    end
end
% Close Spiky windows
delete(hWin)
CloseSpikyWindows
% Clear global variables
clear g_hSpike_Viewer g_bMergeMode
clear OpenSettings % clear persistent variables in OpenSettings()
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function QuickSort(varargin)
% QuickSort! is a function for automatically detectingspikes, dejittering, removing
% outliers and sorting all currently displayed channels. Its a tool for quickly
% giving an impression of the data. Some steps are skipped in Merge Mode. It assumes
% thresholds has already been define.
[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end
if ~isfield(FV, 'tData'), return, end % check data has been loaded
if ~IsMergeMode() % skip spike detection if we're in Merge Mode
    DetectSpikes()
    if isempty(FV.tSpikeThresholdsNeg) % check thresholds have been set
        uiwait(warndlg('Threshold have not been defined. QuickSort! was aborted.', 'Spiky'))
        return
    end
end
DejitterSpikes()
RemoveOutlierSpikes()
SortAllChannels()
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sTxt] = ShowCurrentLine(obj, event_obj)
% Highlight a line when clicked on and bring to front
sTxt = {['']};

hObj = gco;
sTypeOfCurrObj = get(hObj, 'Type');
if ~strcmpi(sTypeOfCurrObj, 'line'), return, end

% Reset color on unclicked lines
hToggleHandles = findobj(gcf, 'Tag', get(hObj, 'tag'), 'Type', 'line');
cCol = get(hToggleHandles(2:end), 'color');
if iscell(cCol)
    mCol = reshape([cCol{:}], 3, length(hToggleHandles)-1)';
else
    mCol = reshape(cCol, 3, length(hToggleHandles)-1)';
end
set(hToggleHandles(2:end), 'color', median(mCol,1), 'linewidth', 1, 'marker', 'none');

% Reset marker on all line objects
hHandles = findobj(gcf, 'Type', 'line');
set(hHandles, 'marker', 'none')

% Change thickness and add markers to clicked line
set(hObj, 'color', median(mCol,1), 'linewidth', 3, 'marker', 's', 'markeredgecolor', 'w', 'markersize', 4)

% Move clicked line to top of stack
hAx = get(gco, 'parent');
hHandles = get(hAx, 'children');
set(hAx, 'children', [hObj; setdiff(hHandles, hObj)]);

% If an UIHoldFigure figure is open, close it/them
close(findobj('Tag', 'UIHoldFigure'));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IdentifySpike(varargin)
[FV, ~] = GetStruct();
sAxTag = get(gca, 'Tag');
vSpikeTimes = eval(['FV.tSpikes.' sAxTag '.spiketimes']); % datapoints
nFs = eval(['FV.tData.' sAxTag '_KHz']); % Fs, kHz
vXY = get(gca, 'currentPoint');
nSpikeTime = round(vXY(1,1)*nFs*1000); % time of clicked spike (dp's)
[nMin, nIndx] = min(abs(vSpikeTimes - nSpikeTime));
try
    vAssigns = eval(['FV.tSpikes.' sAxTag '.hierarchy.assigns']); % datapoints
    sUnit = num2str(vAssigns(nIndx));
    vUniqeAssigns = unique(unique(vAssigns));
    nColIndx = find(vUniqeAssigns == vAssigns(nIndx));
    if nColIndx == 1, vCol = [.5 .5 .5];
    else vCol = FV.mColors(nColIndx, :); end
catch
    sUnit = 'not sorted';
    vCol = [1 1 1];
end
hDot = plot(vXY(1,1), vXY(1,2), 'o', 'markersize', 10, 'Color', vCol, 'linewidth', 2);
sTxt = sprintf('   Spike %d\n   Unit %s', nIndx, sUnit);
text(vXY(1,1), vXY(1,2), sTxt, 'color', vCol)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CloseSpikyWindows(varargin)
hWin = get(0, 'children');
for i = 1:length(hWin)
    if ~isempty(strfind(lower(get(hWin(i), 'name')), 'spiky')) && ~strcmp(get(hWin(i), 'Tag'), 'Spiky')
        close(hWin(i))
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BringSpikyFiguresToFront(varargin)
hWin = get(0, 'children');
for i = 1:length(hWin)
    if ~isempty(strfind(lower(get(hWin(i), 'name')), 'spiky')) && ~strcmp(get(hWin(i), 'Tag'), 'Spiky')
        figure(hWin(i))
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CascadeSpikyWindows(varargin)
hWin = get(0, 'children');
% Get handles of windows to cascade
hTileWin = [];
for i = 1:length(hWin)
    if ~isempty(strfind(lower(get(hWin(i), 'name')), 'spiky')) && ~strcmp(get(hWin(i), 'Tag'), 'Spiky')
        hTileWin(end+1) = hWin(i);
    end
end
if isempty(hTileWin), return, end
vScreenSize = get(0, 'screenSize'); % screen size
nWidth = vScreenSize(3);
nHeight = vScreenSize(4);
% Tile windows
nStep = 25;
for h = 1:length(hTileWin)
    set(hTileWin(h), 'position', [h*nStep (nHeight-420)-(h+2)*nStep 560 420], 'toolbar', 'figure')
    figure(hTileWin(h))
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cSelectedFields, nManualValue, nButton] = SelectImportVariable(cFields, sTitle, ...
    csButton, bManualOpt, nDefManVal, sDefField, bMultiSel, varargin)
% This is a support function for ImportData. It prompts the user to select one
% or more variable from a listbox.
%
% Inputs:
%   cFields     cell of fieldnames
%   sTitle      Window title
%   csButton    button text
%   bManualOpt  enable option to set value manually (0/1)
%   nDefManVal  default manual value
%   sDefField   default selected field (needs to be empty if not used)
%   bMultiSel   option to enable/diable multi-select
%   cDescr      cell of fieldnames descriptions
%
nFields = length(cFields);

if bManualOpt, cFields{end+1} = 'Set Manually'; end
if bManualOpt && ~isempty(nDefManVal) && ~isnan(nDefManVal)
    cFields{end+1} = sprintf('Default (%d)', nDefManVal);
end

% Modify fieldnames to include description strings
if ~isempty(varargin)
    cDescr = varargin{1};
    cFieldsMod = {};
    for i = 1:nFields
        if isempty(cDescr{i}) || strcmp(cDescr{i}, ' ')
            cFieldsMod{i} = cFields{i};
        else
            cFieldsMod{i} = sprintf('%s (%s)', cFields{i}, cDescr{i});
        end
    end
    cFieldsMod{end+1} = cFields{end};
else cFieldsMod = cFields; end

% Default selected field
vDefField = [];
if isempty(sDefField), vDefField = 1;
else
    if iscell(sDefField)
        for c = 1:length(sDefField)
            vIndx = find(strcmpi(cFields, sDefField{c}));
            if ~isempty(vIndx)
                vDefField(end+1) = vIndx;
            end
        end
    else
            vIndx = find(strcmpi(cFields, sDefField));
            if ~isempty(vIndx)
                vDefField(end+1) = vIndx;
            end
    end
end

if bMultiSel, nMax = 10;
else nMax = 1; end

hFig = figure;
nElHeight = 16;
vPos = get(hFig, 'position');
set(hFig, 'Menu', 'none', 'Name', sTitle, 'NumberTitle', 'off', 'Position', [vPos(1) vPos(2) 300 length(cFields)*nElHeight+25]);

hList = uicontrol(hFig, 'Position', [0 26 300 length(cFields)*nElHeight], 'Style', 'listbox', ...
    'String', cFieldsMod, 'Tag', 'ImportListBox', 'Max', nMax);
if ~isempty(vDefField), set(hList, 'Value', vDefField); end

for b = 1:length(csButton)
    nW = 300 / length(csButton);
    nX = nW * (b - 1);
    uicontrol(hFig, 'Position', [nX 0 nW 30], 'Style', 'pushbutton', 'String', csButton{b}, ...
        'Callback', 'global nSelectedIndx sButton; sButton = get(gcbo, ''string''); nSelectedIndx = get(findobj(''Tag'', ''ImportListBox''), ''value''); close');
end

uiwait(hFig)
global nSelectedIndx sButton
nButton = find(strcmp(csButton, sButton)); % get ID of button that was pressed

if isempty(nSelectedIndx)
    cSelectedFields = [];
    nManualValue = [];
    nButton = [];
    return
end

cSelectedFields = cFields(nSelectedIndx);
nManualValue = NaN;
if length(cSelectedFields) == 1
    if strcmp(cSelectedFields(1), 'Set Manually')
        nManualValue = GetInputNum(sTitle, cSelectedFields{1}, nDefManVal);
    elseif ~isempty(cell2mat(strfind(cSelectedFields, 'Default')))
        nManualValue = nDefManVal;
    else
        nManualValue = NaN;
    end
end

clear global nSelectedIndx
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = ChannelProperties(varargin)
if ~IsDataLoaded, return, end
[FV, ~] = GetStruct();
sCh = get(varargin{1},'tag'); % channel name

% Sampling rate
if isfield(FV.tData, [sCh '_KHz'])
    nFs = FV.tData.([sCh '_KHz']); % Hz
else nFs = NaN; end

% Samples
nLen = length(FV.tData.(sCh)); % samples

% Gain
if isfield(FV.tGain, sCh)
    nGain = FV.tGain.(sCh);
else nGain = NaN; end

% Begin time
if isfield(FV.tData, [sCh '_TimeBegin'])
    nTimeBegin = FV.tData.([sCh '_TimeBegin']); % sec
else
    nTimeBegin = NaN;
end

sDescr = GetChannelDescription(sCh);

cData(1, :) = {'Name', sprintf('%s', sCh)};
cData(end+1, :) = {'Description', sprintf('%s', sDescr)};
cData(end+1, :) = {'Sampling rate (Hz)', sprintf('%.2f', nFs*1000)};
cData(end+1, :) = {'Length (samples)', sprintf('%d', nLen)};
cData(end+1, :) = {'Length (s)', sprintf('%.2f', nLen/(nFs*1000))};
cData(end+1, :) = {'Gain', sprintf('%d', nGain)};
cData(end+1, :) = {'Begin time', sprintf('%.4f s', nTimeBegin)};

% Add custom properties
if isfield(FV.tData, [sCh '_Properties'])
    tProps = FV.tData.([sCh '_Properties']);
    for i = 1:length(tProps)
        if isinteger(tProps(i).Value),  sFormat = '%d';
        elseif ischar(tProps(i).Value), sFormat = '%s';
        else                            sFormat = '%.2f'; end
        cData(end+1, :) = {tProps(i).Descr, sprintf(sFormat, tProps(i).Value)};
    end
end

ShowTable(cData, {'Field' 'Value'}, {'' ''}, [0 0], [200 200], ...
    'Channel Properties', [435 min([600 size(cData,1)*30])]);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = DigitizeChannelAuto(hObject, varargin)
DigitizeChannel(hObject, 'auto')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = DigitizeChannelCrossing(hObject, varargin)
DigitizeChannel(hObject, 'crossing')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = DigitizeChannel(varargin)
% Digitize a channel interactively or automatically.
% 
% Usage:
%  DigitizeChannel(H)
%   or
%  DigitizeChannel(H, 'auto/crossing')
%   where   H is a handle to any object that references the channel name in its Tag property
%           'auto'      automatically finds a threshold
%           'crossing'  returns only threshold crossing events
% 
%  [UP DOWN] = DigitizeChannel(S, Y, FS, B, E)
%   where   S is the signal
%           Y the threshold
%           FS the sampling rate (Hz)
%           B the first sample to use (s)
%           E the last sample to use (s)
%

[FV, hWin] = GetStruct();

% Check if data and threshold was supplied with function inputs
if length(varargin) == 5 && isnumeric(varargin{1})
    bInteractive = 0;
elseif length(varargin) == 2 && strcmp('auto', varargin{2})
    bInteractive = 1;
    nY = .5; % threshold for automatic digitization set at 0.5 V
elseif length(varargin) == 2 && strcmp('crossing', varargin{2})
    bInteractive = 1;
else
    bInteractive = 1; % interactive mode; user clicks to set threshold
end

if bInteractive
    % Get data and threshold interactively
    if ~IsDataLoaded, return, end

    % Select channel
    if ~exist('nY')
        zoom off;
        [~, nY] = ginput(1);
        zoom xon
    end
    sTag = get(gca, 'Tag');
    vContData = FV.tData.(sTag);
    nFs = FV.tData.([sTag '_KHz']) * 1000; % sampling frequency (Hz)
    nBeginTime = FV.tData.([sTag '_TimeBegin']); % start of sampling (sec)
    nEndTime = FV.tData.([sTag '_TimeEnd']); % start of sampling (sec)

    % Arbitrary channel operation
    vContData = ChannelCalculator(vContData, sTag);
else
    % Get data and threshold from function input
    vContData = varargin{1};
    nY = varargin{2};
    nFs = varargin{3}; % Hz
    nBeginTime = varargin{4};
    nEndTime = varargin{5};
end

% Threshold signal and find UP and DOWN times
vIndxHIGH = find(vContData >= nY);
vIndxLOW = find(vContData < nY);

% Check UP and DOWN times if any were found
if ~isempty(vIndxHIGH) && ~isempty(vIndxLOW)
    vI = (find(diff(vIndxHIGH) > 1) + 1);
    vIndxUP = vIndxHIGH(unique([1; vI(:)]));
    vI = (find(diff(vIndxLOW) > 1) + 1);
    vIndxDOWN = vIndxLOW(unique([1; vI(:)]));
    
    % Drop index==1 (because no event could have occurred in first frame)
    vIndxUP(vIndxUP == 1) = [];
    vIndxDOWN(vIndxDOWN == 1) = [];
    
    % Compute corresponding time vector
    vTime = linspace(nBeginTime, nEndTime, length(vContData));
    
    % Convert indices to temporal values
    vUpTimes = vTime(vIndxUP);      % sec
    vDownTimes = vTime(vIndxDOWN);  % sec
    
    % Subtract 1/2 timesample from indices; the mid-point is the best
    % estimate for when the signal actually changed state
    vUpTimes = sort(vUpTimes - ((1/nFs)/2));
    vDownTimes = sort(vDownTimes - ((1/nFs)/2));
    
    % Force signal to go DOWN at end of each file (if this is a merged file)
    if isfield(FV, 'tData')
        if isfield(FV.tData, 'FileEnd')
            vDownTimes = sort([vDownTimes [FV.tData.FileEnd.Timestamp]]);
        end
    end
    
    % For each UP event, find first DOWN event that follows
    nDownMatches = zeros(size(vUpTimes));
    for u = 1:length(vUpTimes)
        vDiff = vUpTimes(u) - vDownTimes;
        vDiff(vDiff > 0) = nan;
        
        [~, nMaxI] = max(vDiff);
        if ~isempty(nMaxI)
            nDownMatches(u) = vDownTimes(nMaxI);
        end
    end
    vDownTimes = nDownMatches;
else
    vUpTimes = [];
    vDownTimes = [];
end

if strcmp('crossing', varargin{2})
    vDownTimes = vUpTimes + 0.005; % 5 ms events
end

if bInteractive
    % Save UP/DOWN times in data structure
    FV.tData.([sTag '_Up']) = vUpTimes;         % UP sample times
    FV.tData.([sTag '_Down']) = vDownTimes;     % DOWN sample times
    FV.csDisplayChannels = {FV.csDisplayChannels{~strcmp(FV.csDisplayChannels, sTag)}};
    FV.csDigitalChannels = unique([FV.csDigitalChannels sTag]);

    % Save threshold value used for digitization
    if ~isfield(FV, 'tEventThresholds')
        FV.tEventThresholds = struct([]);
        nIndx = 1;
    else
        % Remove old thresholds
        FV.tEventThresholds(strcmp(sTag, {FV.tEventThresholds.sChannel})) = [];
        nIndx = length(FV.tEventThresholds) + 1;
    end
    
    FV.tEventThresholds(nIndx).sChannel = sTag;
    FV.tEventThresholds(nIndx).nThreshold = nY;

    SetStruct(FV)
    ViewTrialData()
else
    % Store results in output variable
    varargout{1} = vUpTimes;
    varargout{2} = vDownTimes;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DuplicateChannel(varargin)
% Duplicate a channel with options for resampling and filtering
% 
% Usage:
%   DuplicateChannel()
% 
% Batch = Yes
% 
[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end

% Select channel and get resampling parameters interactively
global g_bBatchMode
persistent p_sCh p_nHiPass p_nLoPass p_bRectify p_sDescr
if ~g_bBatchMode || isempty(p_sCh)
    [p_sCh, ~] = SelectChannelNumber(FV.csChannels, 'Select channel', p_sCh);
    if isempty(p_sCh), return, end
    if isempty(p_nHiPass); p_nHiPass = NaN; end
    if isempty(p_nLoPass); p_nLoPass = NaN; end
    if isempty(p_bRectify); p_bRectify = 0; end
    if isempty(p_sDescr); p_sDescr = ''; end
    
    cPrompt = {'Description','High-pass (Hz)','Low-pass (Hz)','Rectify (1=yes, 0=no)'};
    cAnswer = inputdlg(cPrompt,'Resampling options', 1, ...
        {p_sDescr, num2str(p_nHiPass), num2str(p_nLoPass), num2str(p_bRectify)});
    if isempty(cAnswer), return, end
    p_sDescr = cAnswer{1};
    p_nHiPass = str2double(cAnswer{2});
    p_nLoPass = str2double(cAnswer{3});
    p_bRectify = str2double(cAnswer{4});
end

% Get channel data
vCont = FV.tData.(p_sCh);
nFs = FV.tData.([p_sCh '_KHz']) * 1000;

% Filter channel
if ~(p_nHiPass == 0 && p_nLoPass == nFs)
    [vContNew, ~, nFsNew] = FilterChannel(vCont, [], nFs, p_nLoPass, p_nHiPass, p_bRectify, 'decimate');
else
    vContNew = vCont;
end

FV.tData.([p_sCh '_Copy_Imported']) = 1;
FV.tData.([p_sCh '_Copy']) = vContNew;
FV.tData.([p_sCh '_Copy_KHz']) = nFsNew / 1000;
FV.tData.([p_sCh '_Copy_KHz_Orig']) = nFsNew / 1000;
FV.tData.([p_sCh '_Copy_TimeBegin']) = FV.tData.([p_sCh '_TimeBegin']);
FV.tData.([p_sCh '_Copy_TimeEnd']) = FV.tData.([p_sCh '_TimeEnd']);

FV.csDisplayChannels = unique([FV.csDisplayChannels [p_sCh '_Copy']]);
sNewCh = [p_sCh '_Copy'];
FV.csChannels = unique([FV.csChannels sNewCh]);

% Save channel description
iCh = find(strcmp(sNewCh, {FV.tChannelDescriptions.sChannel}));
if isempty(iCh)
    FV.tChannelDescriptions(end+1).sChannel = sNewCh;
    FV.tChannelDescriptions(end).sDescription = p_sDescr;
else
    FV.tChannelDescriptions(iCh(1)).sChannel = sNewCh;
    FV.tChannelDescriptions(iCh(1)).sDescription = p_sDescr;
    % Remove duplicates
    if length(iCh) > 1
        FV.tChannelDescriptions(iCh(2:end)) = [];
    end
end

% Copy channel gain
FV.tGain.(sNewCh) = FV.tGain.(p_sCh);

if ~g_bBatchMode; BatchRedo([], 'DuplicateChannel'); end
SetStruct(FV)
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DeleteChannel(varargin)
% Delete a channel. If input arguments are provided, the first string
% encountered will be the channel to delete.
[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end

sCh = '';
if nargin > 0
    for i = 1:nargin
        if ischar(varargin{i}) sCh = varargin{i}; end
    end
end
if ~any(strcmp(sCh, [FV.csDisplayChannels FV.csChannels]))
    csAll = [FV.csDisplayChannels FV.csChannels];
    [sCh, ~] = SelectChannelNumber(csAll);
    if isempty(sCh), return, end % no channel selected
end

sAns = questdlg(['Are you sure you want to delete channel ' sCh '?'], 'Spiky', 'Yes', 'No', 'No');
if strcmp(sAns, 'No'), return; end

% Delete fields
if isfield(FV.tData, sCh), FV.tData = rmfield(FV.tData, sCh); end
if isfield(FV.tData, [sCh '_KHz']), FV.tData = rmfield(FV.tData, [sCh '_KHz']); end
if isfield(FV.tData, [sCh '_KHz_Orig']), FV.tData = rmfield(FV.tData, [sCh '_KHz_Orig']); end
if isfield(FV.tData, [sCh '_TimeBegin']), FV.tData = rmfield(FV.tData, [sCh '_TimeBegin']); end
if isfield(FV.tData, [sCh '_TimeEnd']), FV.tData = rmfield(FV.tData, [sCh '_TimeEnd']); end
if isfield(FV.tData, [sCh '_Imported']), FV.tData = rmfield(FV.tData, [sCh '_Imported']); end
if isfield(FV.tData, [sCh '_Up']), FV.tData = rmfield(FV.tData, [sCh '_Up']); end
if isfield(FV.tData, [sCh '_Down']), FV.tData = rmfield(FV.tData, [sCh '_Down']); end

FV.csDisplayChannels = unique({FV.csDisplayChannels{~strcmp(FV.csDisplayChannels, sCh)}});
FV.csDigitalChannels = unique({FV.csDigitalChannels{~strcmp(FV.csDigitalChannels, sCh)}});
FV.csChannels = unique({FV.csChannels{~strcmp(FV.csChannels, sCh)}});

FV.tChannelDescriptions = FV.tChannelDescriptions(~strcmp({FV.tChannelDescriptions(:).sChannel}, sCh));
FV.tChannelDescriptions(strcmp({FV.tChannelDescriptions(:).sChannel}, sCh)) = [];

SetStruct(FV);
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UndigitizeChannel(varargin)
% Undigitizes channels that were previously digitized through the GUI
[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end

[sCh, ~] = SelectChannelNumber(FV.csDigitalChannels);
if isempty(sCh), return, end % no channel was selected (i.e. dialog was closed)

if ~isfield(FV.tData, sCh)
    % Check whether channel has an analog equivalent. If not, we cannot 'undigitize'it
    warndlg(['Channel ' sCh ' cannot be un-digitized as it has no analog equivalent.'], 'Spiky')
    return
end

% 'Un-digitize' channel
FV.csDigitalChannels(strcmpi(FV.csDigitalChannels, sCh)) = [];
SetStruct(FV)
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ModifyEventDuration(varargin)
% Changes the duration of all Up-Down events on selected channel to set value
[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end
persistent p_nNewDur p_sCh
global g_bBatchMode

if ~g_bBatchMode || isempty(p_sCh)
    [p_sCh, ~] = SelectChannelNumber(FV.csDigitalChannels);
end
if isempty(p_sCh), return, end % no channel was selected (i.e. dialog was closed)
if ~isfield(FV.tData, [p_sCh '_Up']); return; end

% Current event duration (average)
vUp = FV.tData.([p_sCh '_Up']);
vDown = FV.tData.([p_sCh '_Down']);
if isempty(vUp) || isempty(vDown), return; end
if length(vUp) ~= length(vDown)
    sp_disp(['Channel ' p_sCh ' has an unequal number of UP and DOWN events. Nothing done...'])
    return;
end

% Get all event durations
vDurs = round((vDown - vUp).*1000)./1000; % s (rounded to nearest ms)
if length(unique(vDurs)) > 1
    sp_disp(['Channel ' p_sCh ' events have different durations.'])
end

% Ask for new event duration
if ~g_bBatchMode || isempty(p_nNewDur)
    cAns = inputdlg('New event duration (ms):', 'Spiky', 1, {num2str(p_nNewDur)});
    if isempty(cAns), return, end
    p_nNewDur = str2double(cAns{1}); % ms
end
if isempty(p_nNewDur) return; end

% Update events durations
vNewDown = vUp + (p_nNewDur / 1000);
FV.tData.([p_sCh '_Down']) = vNewDown;

if ~g_bBatchMode; BatchRedo([], 'ModifyEventDuration'); end
SetStruct(FV)
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CreateVirtualEvent(varargin)
% Create a virtual event
% 
% Virtual events are copies of existing events or combinations of two
% events computed using the logical operators AND, OR, XOR and NOT.
%
% Usage:
%   CreateVirtualEvent()    Display GUI to create event
%

[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end

hFig = figure('visible', 'off', 'position', [0 0 300 100], 'WindowStyle', 'modal');
centerfig(hFig, GetGUIHandle())
ThemeObject(hFig, 'visible', 'on', 'toolbar', 'none', 'menubar', 'none', ...
    'units', 'normalized', 'name', 'Create Virtual Event', ...
    'closeRequestFcn', 'set(gcbf,''userdata'',0)');

% Text
csTxt = fliplr({'Channel 1' 'Operator' 'Channel 2' 'New name'});
for i = 1:length(csTxt)
    hTxt(i) = uicontrol(hFig, 'Style', 'text', 'units', 'normalized', ...
        'Position', [0 .25*(i-1) .25 .25], 'String', csTxt{i});
end
ThemeObject(hTxt)

% Channel selectors
vY = [.25 .75];
for i = 1:2
    hUI(i) = uicontrol(hFig, 'Style', 'popup', 'units', 'normalized', ...
        'Position', [.25 vY(i) .75 .25], ...
        'String', GetChannelDescription(FV.csDigitalChannels));
end

% Operator and NOT selectors
csOps = {'AND', 'OR' 'XOR'};
hUI(3) = uicontrol(hFig, 'Style', 'popup', 'units', 'normalized', 'Position', [.25 .5 .25 .25], ...
    'String', csOps );
hUI(4) = uicontrol(hFig, 'Style', 'popup', 'units', 'normalized', 'Position', [.5 .5 .25 .25], ...
    'String', {'' 'NOT'} );

% New name input
hUI(5) = uicontrol(hFig, 'Style', 'edit', 'units', 'normalized', 'Position', [.25 0 .5 .25]);

% Create button
uicontrol(hFig, 'Style', 'pushbutton', 'units', 'normalized', 'Position', [.75 0 .25 .25], ...
    'String', 'Create', 'Callback', 'set(gcf,''userdata'', 1)');

waitfor(hFig, 'UserData') % wait for Create button to be pressed
if get(hFig, 'UserData') == 0
    delete(hFig);
    return
end

% Get channels and operator
csCh1 = FV.csDigitalChannels{get(hUI(2), 'value')};
csCh2 = FV.csDigitalChannels{get(hUI(1), 'value')};
sOp = csOps{get(hUI(3), 'value')};
if get(hUI(4), 'value') == 2, sNot = '~';
else sNot = ''; end
sName = get(hUI(5), 'string');

delete(hFig)
SpikyWaitbar(0, 10)

% Get event data
[mEvents1(:,1), mEvents1(:,2)] = GetEventPairs(csCh1, 'all'); % s
[mEvents2(:,1), mEvents2(:,2)] = GetEventPairs(csCh2, 'all'); % s

% Convert to samples
nFs = max([FV.tData.([csCh1 '_KHz']) FV.tData.([csCh2 '_KHz'])]) * 2500; % Hz
mEvents1 = round(mEvents1 * nFs); % samples
mEvents2 = round(mEvents2 * nFs); % samples

nStart = min([mEvents1(:,1); mEvents2(:,1)]) - 1; % samples
nEnd = max([mEvents1(:,2); mEvents2(:,2)]); % samples
mEvents1 = mEvents1 - nStart;
mEvents2 = mEvents2 - nStart;
SpikyWaitbar(2, 10)

% Expand UP/DOWN indices to a continuous digital vector
mCh = sparse(nEnd - nStart, 2);

mCh(mEvents1(:,1), 1) = mEvents1(:,1); % channel 1
mCh(mEvents1(:,2), 1) = -mEvents1(:,1);
mCh(cumsum(mCh(:,1)) > 0, 1) = 1;
SpikyWaitbar(3, 10)

mCh(mEvents2(:,1), 2) = mEvents2(:,1); % channel 2
mCh(mEvents2(:,2), 2) = -mEvents2(:,1);
mCh(cumsum(mCh(:,2)) > 0, 2) = 1;
SpikyWaitbar(4, 10)

mCh(mCh < 0) = 0;

% Create a sparse matrix where events are represented by ones
%for u = 1:size(mEvents1, 1)
%    mCh(mEvents1(u,1):mEvents1(u,2), 1) = 1;
%end
%for u = 1:size(mEvents2, 1)
%    mCh(mEvents2(u,1):mEvents2(u,2), 2) = 1;
%end
%%

eval(sprintf('vNew = %s(mCh(:,1), %smCh(:,2));', lower(sOp), sNot))
SpikyWaitbar(5, 10)

% Digitize
vUp = find(diff(vNew) == 1) + nStart + 1; % samples
vDown = find(diff(vNew) == -1) + nStart + 1; % samples

% Convert to seconds
vUp = vUp ./ nFs; % s
vDown = vDown ./ nFs; % s

FV.tData.([sName '_KHz']) = nFs / 1000;
FV.tData.([sName '_TimeBegin']) = vUp(end);
FV.tData.([sName '_TimeEnd']) = vDown(end);
FV.tData.([sName '_Up']) = vUp';
FV.tData.([sName '_Down']) = vDown';

FV.csDigitalChannels = unique([FV.csDigitalChannels sName]);
if isfield(FV, 'tChannelDescriptions')
    if ~any(strcmp([FV.tChannelDescriptions(:).sChannel], sName))
        FV.tChannelDescriptions(end+1).sChannel = sName;
        FV.tChannelDescriptions(end).sDescription = sName;
    end
end
SpikyWaitbar(8, 10)

SetStruct(FV)
ViewTrialData()

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ModifyEventInterval(varargin)
% Set a minimum interval between events on selected channel. Events that
% occur after a preceding event at less than the minimum interval are deleted.
[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end
persistent p_nMinInt p_sCh
global g_bBatchMode

if ~g_bBatchMode || isempty(p_sCh)
    [p_sCh, ~] = SelectChannelNumber(FV.csDigitalChannels);
end
if isempty(p_sCh), return, end % no channel was selected (i.e. dialog was closed)
if ~isfield(FV.tData, [p_sCh '_Up']); return; end

% Current event duration (average)
vUp = FV.tData.([p_sCh '_Up']);
vDown = FV.tData.([p_sCh '_Down']);
if isempty(vUp) || isempty(vDown), return; end
if length(vUp) ~= length(vDown)
    sp_disp(['Channel ' p_sCh ' has an unequal number of UP and DOWN events. Nothing done...'])
    return;
end

% Get event intervals
vInt = diff(vUp) .* 1000; % ms

% Get minimum interval interactively
if ~g_bBatchMode || isempty(p_nMinInt)
    cAns = inputdlg('Minimum event interval (ms):', 'Spiky', 1, {num2str(p_nMinInt)});
    if isempty(cAns), return, end
    p_nMinInt = str2double(cAns{1}); % ms
end
if isempty(p_nMinInt); return; end

% Remove events that violate the minimum interval
vRem = [false vInt < p_nMinInt];
vUp(vRem) = [];
vDown(vRem) = [];

FV.tData.([p_sCh '_Up']) = vUp;
FV.tData.([p_sCh '_Down']) = vDown;

if ~g_bBatchMode; BatchRedo([], 'ModifyEventInterval'); end
SetStruct(FV)
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bResult = SpikyWaitbar(nStatus, nLen)
% Create and display a custom waitbar in the toolbar of the Spiky GUI.
% 
% Usage: SpikyWaitbar(0, 20)    Initialize waitbar with length 20
%        SpikyWaitbar(10, 20)   Waitbar reaches 10 of 20
%        SpikyWaitbar(20, 20)   Waitbar reaches end and is destroyed
%        SpikyWaitbar(0, 0)     Function returns current position (0 if no bar exists)
% 
% Function returns false if the Cancel button is pressed.
%
global g_hSpike_Viewer
persistent p_nStatus

if nStatus == 0 && nLen == 0
    bResult = p_nStatus;
end

bResult = true;
if nLen > 30
    nStatus = ceil(nStatus / ceil(nLen/30));
    nLen = ceil(nLen / ceil(nLen/30));
end

hToolbar = findobj(g_hSpike_Viewer, 'Tag', 'Spiky_Waitbar');
if nStatus == 0 % Initialize
    % Hide Action buttons in toolbar
    set(findobj(hToolbar, '-regexp', 'Tag', 'Spiky_WaitbarAction_'), 'visible', 'off')
    % Dont re-initialize if a progress bar already exists
    if isempty(findobj(hToolbar, 'Tag', 'Spiky_Waitbar_Cancel'))
        % Cancel button
        mCData = im2double(imread('cancel.png'));
        mCData(mCData == 1) = NaN;
        uitoggletool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_Waitbar_Cancel');
        cCData(:,:,1) = repmat(1,[16 16]); % red
        cCData(:,:,2) = repmat(0,[16 16]); % green
        cCData(:,:,3) = repmat(0,[16 16]); % blue
        for i = 1:nLen
            uipushtool('Parent', hToolbar, 'cdata', cCData, 'Enable', 'off', 'Tag', sprintf('Spiky_Waitbar_%d', i));
        end
    end
    p_nStatus = nStatus;
else % Update
    hToolbar = findobj(g_hSpike_Viewer, 'Tag', 'Spiky_Waitbar');
    set(get(hToolbar, 'children'), 'separator', 'off')
    for i = 1:nStatus
        set(findobj(hToolbar, 'Tag', sprintf('Spiky_Waitbar_%d', i)), 'Enable', 'On');
        p_nStatus = 0;
    end
    for j = (nStatus+1):nLen
        set(findobj(hToolbar, 'Tag', sprintf('Spiky_Waitbar_%d', j)), 'Enable', 'Off');
    end
    if nStatus >= nLen % Terminate
        hButtons = findobj(hToolbar, '-regexp', 'Tag', 'Spiky_Waitbar_');
        delete(hButtons)
        % Show Action buttons in toolbar
        set(findobj(hToolbar, '-regexp', 'Tag', 'Spiky_WaitbarAction_'), 'visible', 'on')
        p_nStatus = 0;
    end
    p_nStatus = nStatus;
end
hCancelButton = findobj(g_hSpike_Viewer, 'Tag', 'Spiky_Waitbar_Cancel');
if ~isempty(hCancelButton)
    if strcmp(get(hCancelButton(1), 'State'), 'on')
        bResult = false; % process running waitbar was cancelled
    end
end
drawnow % make sure toolbar is updated
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChannelDescriptions(varargin)
% Set descriptive strings and units for all channels interactively.
% 
% Usage:
%   ChannelDescriptions()
%
if ~IsDataLoaded, return, end
[FV, hWin] = GetStruct();

hFig = figure('closeRequestFcn', 'set(gcbf,''userdata'',1)', 'visible', 'off');
centerfig(hFig, hWin)
set(hFig, 'visible', 'on')
cColumnNames = {'Channel', 'Description', 'Unit'};
cColumnFormat = {{'numeric' 'Adjustable'}, {'numeric' 'Adjustable'}};
cColumnEditable =  [false true true];

% Create list of selectable channels
cFields = fieldnames(FV.tData);
cData = {};
for i = 1:length(cFields)
    if ~isempty(strfind(cFields{i}, '_TimeBegin'))
        cData{end+1, 1} = cFields{i}(1:end-10); % channel name
        
        % Get channel description
        sDescr = '';
        if isfield(FV, 'tChannelDescriptions')
            if ~isfield(FV.tChannelDescriptions, 'sChannel'), nIndx = [];
            else
                nIndx = find(strcmp({FV.tChannelDescriptions.sChannel}, cData{end, 1}));
            end
            % If current channel does not appear in the tChannelDescriptions structure, create it
            if isempty(nIndx)
                FV.tChannelDescriptions(end+1).sChannel = cData(end, 1);
                FV.tChannelDescriptions(end).sDescription = '';
                nIndx = length(FV.tChannelDescriptions);
            end
            sDescr = FV.tChannelDescriptions(nIndx).sDescription;
        end
        cData{end, 2} = sDescr;
        
        % Get channel unit
        if isfield(FV.tData, [cData{end, 1} '_Unit'])
            cData{end, 3} = FV.tData.([cData{end, 1} '_Unit']);
        else
            cData{end, 3} = '';
        end
    end
end

ThemeObject(hFig)
set(hFig, 'Name', 'Spiky Channel Descriptions', 'ToolBar', 'none', 'menuBar','none')
hTable = uitable('Units', 'normalized','Position', [0 0 1 1], 'Data', cData, 'ColumnName', ...
    cColumnNames, 'ColumnEditable', cColumnEditable, 'ColumnWidth',{150, 277, 100});
waitfor(hFig, 'userdata') % wait for figure to be closed
cData = get(hTable, 'data'); % modified data

% Assign new descriptions
if ~isfield(FV, 'tChannelDescriptions')
    FV.tChannelDescriptions = struct([]);
end
for c = 1:size(cData, 1)
    FV.tChannelDescriptions(c).sChannel = cData{c, 1};
    FV.tChannelDescriptions(c).sDescription = cData{c, 2};
end

% Assign new units
for c = 1:size(cData, 1)
    if isempty(cData{c, 3}), continue; end
    FV.tData.([cData{c, 1} '_Unit']) = cData{c, 3};
end

delete(hFig)
SetStruct(FV)
ViewTrialData()
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExperimentDescriptions(varargin)
% Add, remove and edit custom file description variables
%
[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end

hFig = figure('closeRequestFcn', 'set(gcbf,''userdata'',''closed'')');
cColumnNames = {'Variable', 'Value'};
cColumnFormat = {{'char' 'Adjustable'}, {'char' 'Adjustable'}};
cColumnEditable =  [true true];

% create list of selectable channels
if ~isfield(FV, 'tExperimentVariables')
    FV_temp = SetFVDefaults();
    FV.tExperimentVariables = FV_temp.tExperimentVariables;
end
cData = repmat({''},200,2);
for i = 1:length(FV.tExperimentVariables)
    cData{i, 1} = FV.tExperimentVariables(i).sVariable;
    cData{i, 2} = FV.tExperimentVariables(i).sValue;
end

ThemeObject(hFig)
set(hFig, 'Name', 'Spiky Experiment Variables', 'ToolBar', 'none', 'menuBar','none')
hTable = uitable('Units', 'normalized','Position', [0 0 1 1], 'Data', cData, 'ColumnName', ...
    cColumnNames, 'ColumnEditable', cColumnEditable, 'ColumnWidth',{150, 352});
waitfor(hFig, 'userdata', 'closed') % wait for figure to be closed
cData = get(hTable, 'data'); % modified data

FV.tExperimentVariables = struct([]); % clear variables
for c = 1:size(cData, 1)
    if isempty(cData{c,1}) || isempty(cData{c,2}), continue, end
    FV.tExperimentVariables(end+1).sVariable = cData{c, 1};
    FV.tExperimentVariables(end).sValue = cData{c, 2};
end

delete(hFig)
SetStruct(FV)
ViewTrialData()
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DeleteSortingData(varargin)
% Remove all spike sorting results. All assignments of spikes into clusters
% or outliers will be removed. Raw waveforms, or dejittered waveforms, are
% not affected.
%
[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end

sAns = questdlg('All sorting data, including assignment of spikes as outliers, will be removed. Waveforms and dejittering is not affected. Continue?', ...
    'Spiky', 'Yes', 'No', 'No');
if strcmp(sAns, 'No'), return; end

% Iterate over channels containing spikes
cFields = fieldnames(FV.tSpikes);
for c = 1:length(cFields)
    % Remove fields pertaining to sorting results
    if isfield(FV.tSpikes.(cFields{c}), 'overcluster')
        FV.tSpikes.(cFields{c}) = rmfield(FV.tSpikes.(cFields{c}), 'overcluster');
    end
    if isfield(FV.tSpikes.(cFields{c}), 'hierarchy')
        FV.tSpikes.(cFields{c}) = rmfield(FV.tSpikes.(cFields{c}), 'hierarchy');
    end
    if isfield(FV.tSpikes.(cFields{c}), 'outliers')
        FV.tSpikes.(cFields{c}) = rmfield(FV.tSpikes.(cFields{c}), 'outliers');
    end
end
SetStruct(FV)
ViewTrialData()
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vCont, FV] = AdjustChannelGain(FV, vCont, sCh)
% Adjust channel amplitude by specified gain. The returned signal is by
% default in microvolts, unless a different range (V or mV) has been
% selected in View -> Amplitude Unit or via SetAmplitudeUnit()
%

% Adjust gain of signal
bIsGain = 0;
if ~isempty(FV.tGain)
    if isfield(FV.tGain, sCh)
        if isnumeric(FV.tGain.(sCh))
            nGain = FV.tGain.(sCh);
            bIsGain = 1;
        end
    end
end

if ~bIsGain
    % If no gain has been specified for a channel, open 'Set Channel Gains' dialog
    bResetGain = 0;
    if isempty(FV.tGain), bResetGain = 1;
    elseif ~isfield(FV.tGain, sCh), bResetGain = 1; end
    if bResetGain, FV.tGain(1).(sCh) = 1; end % default gain of 1
    SetStruct(FV)
    SetGain('noplot')
    [FV, ~] = GetStruct(); % get updated values
    nGain = FV.tGain.(sCh);
end

vCont = vCont ./ nGain; % divide by gain

% Scale to microvolts (or optionally to V or mV)
if nGain ~= 1
    vCont = vCont .* FV.tAmplitudeUnit.nFactor;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sHash = GetGitHash
% Get current GIT SHA1 hash (unique identifier of currently used commit).
% The local SHA1 hash can be compared to other versions of Spiky to
% determine if versions are the same or not.
%

% Set paths
sWT_dir = which('spiky');
sGITPath = [sWT_dir(1:strfind(sWT_dir, mfilename)-1) '.git' filesep 'refs' filesep 'heads' filesep 'master'];
sVERSIONPath = [sWT_dir(1:strfind(sWT_dir, mfilename)-1) 'VERSION'];

% Get build number from SVN entries files
sHash = 'Unknown';
if exist(sGITPath, 'file')
    fid = fopen(sGITPath);
    sHash = fgetl(fid);
    fclose(fid);
    % Copy hash into VERSION file in top directory
    hFID = fopen(sVERSIONPath, 'w');
    fprintf(hFID, '%s', sHash);
    fclose(hFID);
    sHash = sHash(1:10);
elseif exist(sVERSIONPath, 'file')
    % If the GIT has does not exist in .git, check if VERSION file exists
    hFID = fopen(sVERSIONPath, 'r');
    sHash = char(fread(hFID, 10, 'schar'))';
    fclose(hFID);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MergeChannels(varargin)
% Merge selected channels
%
[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end
persistent p_sName
if isempty(p_sName), p_sName = 'Merged'; end
csCurrentCh = FV.csDisplayChannels;

% Select channels
csChannels = SelectChannels();

% None or one channel selected
if length(csChannels) < 2
    SetStruct(FV)
    ViewTrialData()
    return
end

% Get data of channels
cCont = cell(1, length(csChannels));
cTime = cell(1, length(csChannels));
vFs = nan(1, length(csChannels));
vBegin = nan(1, length(csChannels));
vEnd = nan(1, length(csChannels));
for nCh = 1:length(csChannels) % iterate over imported channels
    vCont = ChannelCalculator(FV.tData.(csChannels{nCh}), csChannels{nCh});
    [vCont, FV] = AdjustChannelGain(FV, vCont, csChannels{nCh});
    
    % Time vector
    nBegin = FV.tData.([csChannels{nCh} '_TimeBegin']); % sampling start, sec
    nEnd = FV.tData.([csChannels{nCh} '_TimeEnd']); % sampling end, sec
    nFs = FV.tData.([csChannels{nCh} '_KHz']) * 1000; % sampling frequency Hz
    vTime = (nBegin+1/nFs):(1/nFs):(nBegin+length(vCont)/nFs); % absolute time, sec
    
    % Apply channel filter
    [vCont, vTimeOut, nFsOut] = GetFilteredChannel(csChannels{nCh}, vCont, 'decimate');
    if ~isempty(vTimeOut), vTime = vTimeOut; end
    if ~isempty(nFsOut), nFs = nFsOut; end
    cCont{nCh} = vCont;
    cTime{nCh} = vTime;
    vFs(nCh) = nFs;
    vBegin(nCh) = nBegin;
    vEnd(nCh) = nEnd;
end

% Check that sampling rate is same for all merging channels
if length(unique(vFs)) > 1
    uiwait(warndlg('Sampling rate of all merging channels must be the same. Aborting merge.', 'Spiky'));
    return
end

% Initialize vector that will contain all data
nDur = max(vEnd) - min(vBegin); % duration, s
nLen = ceil(nFs * nDur); % new length, samples
nFs = nLen / nDur;
vContAll = nan(1, nLen);

% Insert mergin channels
for nCh = 1:length(vFs)
    nBegin = round((vBegin(nCh) - min(vBegin)) * nFs) + 1; % relative start time, samples
    vContAll(nBegin:(nBegin-1+length(cCont{nCh}))) = cCont{nCh};
end

% Get new channel name interactively
cAns = inputdlg('New channel name:', 'Channel Name', 1, {p_sName});
if ~isempty(cAns), p_sName = cAns{1}; end

% Insert new, merged channel into FV
FV.tData.(p_sName) = vContAll;
FV.tData.([p_sName '_TimeBegin']) = min(vBegin);
FV.tData.([p_sName '_TimeEnd']) = min(vBegin) + (nLen / nFs);
FV.tData.([p_sName '_KHz']) = nFs / 1000;

% Display previously viewed channels, minus the merged ones
FV.csDisplayChannels = unique([setdiff(csCurrentCh, csChannels) p_sName]);
SetStruct(FV)
ViewTrialData
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InvertChannel(varargin)
% Select a continuous channel manually and invert it.
%

[FV, ~] = GetStruct();
if ~IsDataLoaded, return, end

% Select channel
[sCh, ~] = SelectChannelNumber(FV.csChannels);

% Invert channel data
FV.tData.(sCh) = FV.tData.(sCh) .* -1;

% Force inverted channel to be displayed
FV.csDisplayChannels = unique([FV.csDisplayChannels sCh]);

SetStruct(FV)
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Undo
% TODO
function Undo(varargin)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Redo(varargin)
% Redo previous action
%
BatchRedo([], 'redo');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ToggleStatus(varargin)
% Toggle status of the current object (eg. Grid, Normal Time etc)
%
[FV, ~] = GetStruct();
switch lower(get(gcbo,'Checked'))
    case 'on'
        set(gcbo, 'Checked', 'off');
    case 'off'
        set(gcbo, 'Checked', 'on');
end
SetStruct(FV)
if ~isempty(FV.csDisplayChannels)
    ViewTrialData
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShiftBeginTime(varargin)
% Shift the BeginTime time of selected channel
%
[FV, ~] = GetStruct();

% Check that tag of handle references a known channel
sCh = get(varargin{1}, 'Tag');
if ~(isempty(sCh) || ~any(strcmp(fieldnames(FV.tData), sCh)))
    % Request timeshift in msec
    nNum = GetInputNum('Milliseconds to shift begin time (can be positive or negative):', 'Shift begin time', 0);

    nBeginTime = FV.tData.([sCh '_TimeBegin']);
    nEndTime = FV.tData.([sCh '_TimeEnd']);

    nBeginTime = nBeginTime + (nNum / 1000); % sec
    nEndTime = nEndTime + (nNum / 1000); % sec

    FV.tData.([sCh '_TimeBegin']) = nBeginTime;
    FV.tData.([sCh '_TimeEnd']) = nEndTime;
end
SetStruct(FV)
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ReplaceValues(varargin)
% Replace a particular value in selected channel. Note that the FROM value
% is searched for after applying define (if any) math operations. The new
% value is, however, inserted into the original data vector.
%
[FV, ~] = GetStruct();
persistent p_cAns
if isempty(p_cAns); p_cAns = {'' ''}; end

% Check that tag of handle references a known channel
sCh = get(varargin{1}, 'Tag');
if ~(isempty(sCh) || ~any(strcmp(fieldnames(FV.tData), sCh)))
    % Get value to replace and new value
    cAns = inputdlg({'Value to replace:' 'New value:'}, 'Replace value', 1, p_cAns);
    if isempty(cAns); return; end
    if isempty(cAns{1}) || isempty(cAns{2}); return; end
    p_cAns = cAns;
    nFrom = str2double(p_cAns{1});
    nTo = str2double(p_cAns{2});
    vDataOrig = FV.tData.(sCh);
    vData = ChannelCalculator(vDataOrig, sCh);
    if isnan(nFrom)
        iChange = isnan(vData);
    else
        iChange = vData == nFrom;
    end
    vDataOrig(iChange) = nTo;
    FV.tData.(sCh) = vDataOrig;
end
SetStruct(FV)
ViewTrialData
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowDAQInfo(varargin)
% Show info stored in DAQ file
%
[FV, ~] = GetStruct();
if isfield(FV, 'sLoadedTrial') % a trial is loaded
    sFilename = FV.sLoadedTrial;
    if ~isempty(sFilename) % we have a filename
        if exist(sFilename, 'file') % file exists
            global tInfo
            tInfo = daqread(sFilename, 'info');
            assignin('base', 'tInfo', tInfo); % copy tInfo to base workspace
            openvar('tInfo')
            return
        end
    end
end
warndlg('No DAQ file loaded or loaded file is not a DAQ (.daq) file.', 'Spiky')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MeasureLine(varargin)
% Measure time between two positions along displayed channel.
%
hObj = varargin{1}; % get calling object
[FV, ~] = GetStruct();

% Find axis with same tag and make that current
sCh = get(hObj, 'Tag');
hAx = findobj('Type', 'axes', 'Tag', sCh);
if isempty(hAx), return
else axes(hAx(1)), end

vCont = ChannelCalculator(FV.tData.(sCh), sCh); % continuous trace (V)
nBeginTime = FV.tData.([sCh '_TimeBegin']); % start of sampling (sec)
nEndTime = FV.tData.([sCh '_TimeEnd']); % start of sampling (sec)
nFs = FV.tData.([sCh '_KHz']) * 1000; % Hz
vTime = GetTime(nBeginTime, nEndTime, vCont, nFs);

hold on

% Point A
[nXa, ~] = ginput(1);
% Snap to nearest data point on time axis
[~, nMinIndx] = min(abs(vTime - nXa));
nXa_samples = nMinIndx;
nXa = vTime(nMinIndx);
nYa = vCont(nMinIndx);
hPnta = plot(nXa, nYa, 'ro');

% Point B
[nXb nYb] = ginput(1);
% Snap to nearest data point on time axis
[~, nMinIndx] = min(abs(vTime - nXb));
nXb_samples = nMinIndx;
nXb = vTime(nMinIndx);
nYb = vCont(nMinIndx);
hPntb = plot(nXb, nYb, 'ro');

% Plot connecting line
hLin = plot([nXa nXb], [nYa nYb], 'r-');

% Compute distances
nXd = nXb - nXa;
nYd = nYb - nYa;

waitfor(msgbox(sprintf('Time: %.2f msec\nLength: %d samples\nAmplitude: %.2f mV', ...
    nXd*1000, nXb_samples-nXa_samples, nYd*1000), 'Measure'))

if ishandle(hPnta), delete(hPnta); end
if ishandle(hPntb), delete(hPntb); end
if ishandle(hLin), delete(hLin); end
hold off
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowEventStatistics(varargin)
% Show the number of UP/DOWN events on each channel.
%
[FV, hWin] = GetStruct();

% Iterate over digital channels
csChannels = FV.csDigitalChannels;

cNames = {'Name', 'kHz', '#Up', '#Down', 'Event Durations (ms)', 'Description'};
cFormat = {{'char' 'Fixed'}, {'numeric' 'Fixed'}, {'numeric' 'Fixed'}, ...
    {'numeric' 'Fixed'}, {'char' 'Fixed'}, {'char' 'Fixed'}};
cEditable =  [1 0 0 0 0];

for nCh = 1:length(csChannels)
    cData{nCh, 1} = csChannels{nCh};
    cData{nCh, 2} = FV.tData.([csChannels{nCh} '_KHz']);
    cData{nCh, 3} = length(FV.tData.([csChannels{nCh} '_Up']));
    cData{nCh, 4} = length(FV.tData.([csChannels{nCh} '_Down']));
    vDown = FV.tData.([csChannels{nCh} '_Down']);
    vUp = FV.tData.([csChannels{nCh} '_Up']);
    iMax = min([length(vDown) length(vUp)]);
    vDur = vDown(1:iMax) - vUp(1:iMax);
    vDur = unique(round(vDur.*1000));
    cData{nCh, 5} = sprintf('%.0f,', vDur);
    cData{nCh, 5} = cData{nCh, 5}(1:end-1);
    if strcmp(cData{nCh, 5}, '0'), cData{nCh, 5} = ''; end
    
    % Get average event duration
    if cData{nCh, 3} == cData{nCh, 4}
        FV.tData.([csChannels{nCh} '_Down']) = FV.tData.([csChannels{nCh} '_Up']);
    end
    nIndx = strcmp(csChannels{nCh}, {FV.tChannelDescriptions.sChannel});
    if ~isempty(FV.tChannelDescriptions(nIndx))
        cData{nCh, 6} = FV.tChannelDescriptions(nIndx).sDescription;
    end
end

vWidths = [150, 100, 70, 70, 160, 220];
ShowTable(cData, cNames, cFormat, cEditable, vWidths, 'Event Statistics', [800 200]);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowTable(cData, cNames, cFormat, vEditable, vWidths, sTitle, vSize)
% Display data in a table
% Tabular data can be exported to a spreadsheet by selecting data and
% copy/paste with Ctrl+C / Ctrl-D
% 
% Usage:
%   ShowTable(C, N, F, E, S, [W H]), where
%       C is a cell containing cell data
%       N
%       F
%       E is a logical vector denoting which columns are editable
%       W is a vector of column widths
%       S is the figure title
%       [W H] is the Width and Height of the table window
% 
%

% Initialize figure
hFig = figure('visible', 'off');
vPos = get(hFig, 'position');
ThemeObject(hFig, 'Name', sTitle, 'ToolBar', 'none', 'menuBar','none', 'position', [vPos(1:2) vSize]);

if exist('centerfig', 'file')
    [~, hWin] = GetStruct();
    centerfig(hFig, hWin);
end

% Create table
uitable('Units', 'normalized', 'Position', [0 0 1 1], ...
    'Data', cData, ...
    'ColumnName', cNames, ...
    'ColumnWidth', num2cell(vWidths), ...
    'ColumnEditable', logical(vEditable) );

set(hFig, 'visible', 'on')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DeleteSection(varargin)
% Delete section of selected channel. This function is only called by the context
% menu which appears when right-clicking channel traces.
%
[FV, ~] = GetStruct();

sCh = get(gcbo, 'tag'); % channel name
hAx = findobj(gcf, 'type', 'axes', 'tag', sCh); % find corresponding axes
hAx = hAx(1);
if isempty(hAx)
    warndlg('Cannot find corresponding axes for selected channel. Aborting...', 'Spiky')
    return
end
axes(hAx) % make axes current
hold on

% Initialize lines
vYLim = get(hAx, 'ylim');
global hLine1 hLine2;
hLine1 = plot([10^16 10^16], vYLim, 'r:'); % first vertical line
hLine2 = plot([10^16 10^16], vYLim, 'r:'); % second vertical line

% Get line 1 click
set(gcf, 'WindowButtonMotionFcn', 'global hLine1;mXY=get(gca,''CurrentPoint'');set(hLine1,''xdata'',mXY(:,1))')
set(hLine1, 'ButtonDownFcn', 'set(gcf,''WindowButtonMotionFcn'',[])')
waitfor(gcf, 'WindowButtonMotionFcn')

% Get line 2 click
set(gcf, 'WindowButtonMotionFcn', 'global hLine2;mXY=get(gca,''CurrentPoint'');set(hLine2,''xdata'',mXY(:,1))')
set(hLine2, 'ButtonDownFcn', 'set(gcf,''WindowButtonMotionFcn'',[])')
waitfor(gcf, 'WindowButtonMotionFcn')

% Ask whether to delete
switch questdlg('Should the selected region and channel be deleted? This action can be undone by reloading original data files from disk.', ...
        'Spiky', 'Delete', 'Cancel', 'Cancel')
    case 'Delete'
        % Get line X coordinates
        vX1 = get(hLine1, 'xdata');
        vX2 = get(hLine2, 'xdata');
        nTT = sort([vX1(1) vX2(1)]); % sort

        % Convert time to datapoints
        nTimeBegin = FV.tData.([sCh '_TimeBegin']);
        nFs = FV.tData.([sCh '_KHz']) * 1000; % Hz
        vData = FV.tData.(sCh);
        nStart = max([1 round((nTT(1) - nTimeBegin) * nFs)]); % start sample
        nEnd = min([length(vData) round((nTT(2) - nTimeBegin) * nFs)]); % end sample

        % Replace deleted values with NaN's
        vData(nStart:nEnd) = NaN;
        FV.tData.(sCh) = vData;
end

% Delete line objects
delete([hLine1 hLine2])

% Release hold on axes
hold off

SetStruct(FV)
ViewTrialData
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClassifyUnit(varargin)
% Classify a selected unit. This function assumes that the calling object
% (gcbo) identifies the channel and unit number in the UserData property.
%
tUserData = get(gcbo, 'UserData');
[FV, ~] = GetStruct();

% Create quality structure if it doesnt already exist
if ~isfield(FV.tSpikes.(tUserData.sChannel), 'quality')
    FV.tSpikes.(tUserData.sChannel).quality = struct([]);
end

% Get index of this unit in .quality sub-structure
nIndx = [];
nScore = 0; % default answer
if ~isempty(FV.tSpikes.(tUserData.sChannel).quality)
    nIndx = find([FV.tSpikes.(tUserData.sChannel).quality.unit] == tUserData.nUnit);
end
if isempty(nIndx), nIndx = length(FV.tSpikes.(tUserData.sChannel).quality) + 1;
else nScore = FV.tSpikes.(tUserData.sChannel).quality(nIndx).score; end

% Get unit quality from user
cPrompt = {sprintf('Quality (0 - 1), where:\n\n  0 = Outliers / noise\n  1 = Multiple units \n  2 = Putative single unit  \n  3 = Single unit with dropped spikes \n  4 = Single unit\n')};
cAns = inputdlg(cPrompt,'Classify Unit', 1, {num2str(nScore)});
if isempty(cAns), return
else
    FV.tSpikes.(tUserData.sChannel).quality(nIndx).unit = tUserData.nUnit;
    FV.tSpikes.(tUserData.sChannel).quality(nIndx).score = str2num(cAns{1});
end
SetStruct(FV)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetElectrodePosition(varargin)
% Set and store the electrode position used during the experiment.
%

sCh = get(gcbo, 'tag');
[FV, ~] = GetStruct();

% Get input from user
nElecIndx = find(strcmp({FV.tChannelDescriptions.sChannel}, sCh));
if isfield(FV.tChannelDescriptions(nElecIndx), 'nAPPos')
    nAPPos = FV.tChannelDescriptions(nElecIndx).nAPPos;
    nMLPos = FV.tChannelDescriptions(nElecIndx).nMLPos;
    nDVPos = FV.tChannelDescriptions(nElecIndx).nDVPos;
    if isfield(FV.tChannelDescriptions(nElecIndx), 'sPosDescr')
        sPosDescr = FV.tChannelDescriptions(nElecIndx).sPosDescr;
    else sPosDescr = ''; end
else nAPPos = []; nMLPos = []; nDVPos = []; sPosDescr = ''; end
if isempty(sPosDescr), sPosDescr = ''; end
cPrompt = {'Anterior-posterior (mm):','Medial-lateral (mm):', 'Dorsal-ventral (mm):', 'Description (e.g. brain region):'};
cAnswer = inputdlg(cPrompt, 'Electrode Position', 1, {num2str(nAPPos), num2str(nMLPos), num2str(nDVPos), sPosDescr});
if isempty(cAnswer), return, end
FV.tChannelDescriptions(nElecIndx).nAPPos = str2num(cAnswer{1});
FV.tChannelDescriptions(nElecIndx).nMLPos = str2num(cAnswer{2});
FV.tChannelDescriptions(nElecIndx).nDVPos = str2num(cAnswer{3});
FV.tChannelDescriptions(nElecIndx).sPosDescr = cAnswer{4};
SetStruct(FV)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetReceptiveField(varargin)
% Set receptive field properties (as a string) of a selected unit.
%
[FV, ~] = GetStruct();
sCh = get(gcbo, 'tag');

% Get unit number
sChUnit = get(get(gcbo, 'parent'), 'tag');
nIndx = findstr(sChUnit, '-');
nUnit = str2double(sChUnit(nIndx+1:end));

% THIS PART CAN BE SAFELY REMOVED
% First, delete old ReceptiveField field if this file contains one
if isfield(FV.tChannelDescriptions, 'sRF')
    sRF = FV.tChannelDescriptions(strcmpi({FV.tChannelDescriptions.sChannel},sCh)).sRF;
    FV.tChannelDescriptions = rmfield(FV.tChannelDescriptions, 'sRF');
    warndlg(sprintf('Removed sRF field in FV.tChannelDescriptions! (RF=%s)', sRF));
end
% END OF REMOVE PART

% Get existing receptive field info
if isfield(FV.tSpikes.(sCh), 'receptivefield')
    nIndx = find([FV.tSpikes.(sCh).receptivefield.unit] == nUnit);
    if isempty(nIndx) sReceptiveField = '';
    else
        sReceptiveField = FV.tSpikes.(sCh).receptivefield(nIndx).rf;
    end
else
    FV.tSpikes.(sCh).receptivefield = struct([]);
    sReceptiveField = '';
end

cAnswer = inputdlg(sprintf('Receptive field for unit %d on electrode %s:\n\nSeparate multiple fields with commas. If there is a\nprimary receptive field, list it first.\n', nUnit, sCh), ...
    'Receptive Field', 1, {sReceptiveField});
if isempty(cAnswer), return, end

% Save new receptive field info
if isfield(FV.tSpikes.(sCh).receptivefield, 'unit')
    nIndx = find([FV.tSpikes.(sCh).receptivefield.unit] == nUnit);
    if isempty(nIndx) nIndx = length(FV.tSpikes.(sCh).receptivefield) + 1; end
else nIndx = 1; end
FV.tSpikes.(sCh).receptivefield(nIndx).unit = nUnit;
FV.tSpikes.(sCh).receptivefield(nIndx).rf = cAnswer{1};
SetStruct(FV)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BatchRedo(hObject, varargin)
% Repeat last recorded action as a batch job on all selected files
% 
% Usage:    BatchRedo(H, S)
%           where H is the calling object handle and S is an action string:
%               'redo'  Repeat last recorded action
%               'xyz'   Record command xyz as action to repeat
%
global g_hSpike_Viewer g_bBatchMode
persistent sLastAction
sAction = varargin{1};
[FV, hWin] = GetStruct();

% List of files
hFileList = findobj(g_hSpike_Viewer, 'Tag', 'MenuFileList');
csFiles = flipud(get(get(hFileList, 'Children'), 'UserData'));
if iscell(csFiles)
    nFiles = length(csFiles);
else
    nFiles = 1;
    csFiles = {csFiles};
end

% If a file selection is not available, or a merge file is currently
% loaded, then repeat last action only on the current loaded file.
[sDir, sFile] = fileparts(FV.sLoadedTrial);
if nFiles == 0 || strcmpi('MergeFile', sFile)
    if ~isempty(FV.sLoadedTrial)
        if isempty(sDir)
            csFiles = {fullfile(FV.sDirectory, FV.sLoadedTrial)};
        else
            csFiles = {FV.sLoadedTrial};
        end
        nFiles = 1;
    end
end

if strcmp(sAction, 'redo')
    % Redo last recorded BATCH action (B)
    if isempty(sLastAction) % Do nothing
        warndlg('There is no available action to run in batch mode. Note that only actions with the B suffix (e.g. in menus or buttons) can run as batch jobs. To start a batch job, you must first run one of the supported actions and then re-select Batch Redo from the menu.')
    else % Repeat last supported redo action
        if nFiles > 1
            sAns2 = questdlg('Should the batch operation be repeated on the currently loaded file?', 'Spiky', 'Yes', 'No', 'No');
        else
            sAns2 = 'Yes';
        end
        
        bWaitResult = SpikyWaitbar(0,nFiles+1);
        g_bBatchMode = true;

        % Get index of current file
        [FV, hWin] = GetStruct();
        sLoadedTrial = FV.sLoadedTrial;
        for m = 1:nFiles
            % Update progress bar
            bWaitResult = SpikyWaitbar(m, nFiles+1);
            if ~bWaitResult; break, end % cancel button pressed
        
            % If applicable, skip current file
            if strcmp(csFiles{m}, sLoadedTrial) && strcmp(sAns2, 'No')
                continue;
            end
            
            if nFiles > 1, OpenFile([], m); end % Load next movie
            eval(sLastAction)   % Redo action
            SaveResults();      % Save results
            drawnow
        end
        if bWaitResult % Batch complete notification
            sStr = sprintf('Batch job %s completed (%d files processed)', sLastAction, nFiles+1);
        else
            sStr = sprintf('Batch job %s cancelled after file %d/%d', sLastAction, m, nFiles+1);
        end
        ViewTrialData(); % refresh
        bWaitResult = SpikyWaitbar(nFiles+1,nFiles+1);
        sp_disp(sStr);
    end
else % Record which action was last performed
    sLastAction = sAction;
    % change toolbar button status
    hBut = findobj(g_hSpike_Viewer, 'Tag', 'Spiky_WaitbarAction_BatchRedo');
    set(hBut, 'TooltipString', sprintf('Repeat %s on all open files', sLastAction), 'enable', 'on');
end
g_bBatchMode = false;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sp_disp(sStr)
% Display a string on the command line.
%
disp(sprintf('%s Spiky says: %s', datestr(now, 'HH:mm:ss'), sStr))
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RunBatchScript(varargin)
% Select a script interactively and run it as a batch job.
%
RunScript('batch')
SaveResults();
BatchRedo([], 'redo')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RunScript(varargin)
% Select and run a Spiky script interactively.
%
%  Usage:
%   RunScript()        select and run a script interactively
%   RunScript('batch') select a script interactively and run as a batch job
%                      on all files
%   RunScript(S)       run the script S (on the Matlab path)
%
[FV, ~] = GetStruct();
global g_bBatchMode

% If the function to be run was not specified as in input, select file manually
if length(varargin) ~= 1
    [sFile, ~] = uigetfile('*.m', 'Pick an M-file');
elseif strcmpi(varargin{1}, 'batch')
    % Run script as batch job

    % Get script file
    [sFile, ~] = uigetfile('*.m', 'Pick an M-file');
    
    % Save script action in Batch Redo function
    if ~g_bBatchMode
        BatchRedo([], sprintf('RunScript(''%s'')', sFile));
    end
    return
else
    sFile = varargin{1};
end
if sFile == 0, return, end

% The remainder of this function runs only if a script is run for one file

% Save action for later global Redo
if ~g_bBatchMode
    BatchRedo([], sprintf('RunScript(''%s'')', sFile));
end

FV.ScriptError = '';
eval(['FV = ' sFile(1:end-2) '(FV);'])
if ~isempty(FV.ScriptError)
    sp_disp(sprintf('In %s: %s', sFile(1:end-2), FV.ScriptError))
end

SetStruct(FV)
ViewTrialData()

% Execute termination commands from script
if isfield(FV, 'ScriptExitCommand')
    eval(FV.ScriptExitCommand)
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowVersion(varargin)
% Display the current version of Spiky.
%
sWT_dir = fileparts(mfilename('fullpath'));
sVERSIONPath = [sWT_dir filesep 'VERSION'];
if ~exist(sVERSIONPath, 'file')
    sHash = 'Unknown';
else
    hFID = fopen(sVERSIONPath, 'r');
    sHash = char(fread(hFID, 10, 'schar'))';
    fclose(hFID);
end
inputdlg('GitHub SHA-1 Checksum:', 'Spiky Version', 1, {sHash});
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowLicense(varargin)
% Display the Spiky license.
%
[~, hWin] = GetStruct();
sWT_dir = fileparts(mfilename('fullpath'));
sLICENSEPath = [sWT_dir filesep 'LICENSE'];
if ~exist(sLICENSEPath, 'file'), return; end
sLicense = textread(sLICENSEPath, '%s', 'whitespace', '', 'bufsize', 2^16);
hFig = figure;
set(hFig, 'name', 'Spiky License', 'menubar', 'none', 'numbertitle', 'off', 'color', 'w')
uicontrol('parent', hFig, 'style', 'edit', 'backgroundcolor', 'w', 'max', 2, ...
    'units', 'normalized', 'position', [0 0 1 1], 'string', sLicense, ...
    'horizontalalignment', 'left')
if exist('centerfig'); centerfig(hFig, hWin); end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = header(x, fontsize)
%
t = findobj(gcf, 'Tag', 'header');
if isempty(t)
	ax = gca;
	axes('Position', [0 0 1 1], 'Visible', 'off');
	t = text(0.5, 0.98, x, ...
				'Units', 'normalized', ...
				'VerticalAlignment', 'top', ...
				'HorizontalAlignment', 'center', ...
				'Tag', 'header', ...
				'FontSize', fontsize);
	axes(ax);
else
	set(t, 'String', x, 'FontSize', fontsize);
end
if nargout > 0; y = t; end;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hPlot,hFill] = mean_error_plot(vMean, vError, vColor, vX)
% THIS FUNCTION IS DEPRECATED. USE PlotMeanError
[hPlot,hFill] = PlotMeanError(vMean, vError, vColor, vX);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hPlot,hFill] = PlotMeanError(vMean, vError, vColor, varargin)
% Plot a mean and error as a solid line with a filled error region in the
% current axes.
% 
% Usage:
%   [HP, HF] = PlotMeanError(M, E, COL, XX)
%   where   HP and HF are handles of the mean and error objects.
%           M is the mean vector
%           E is the error vector
%           COL is the color of the error fill
%           XX is the x-value vector
%
vMean = vMean(:);
vError = vError(:);

if nargin > 3
    vXt = varargin{1}(:);
else
    vXt = (1:length(vError))';
end

vXb = flipud(vXt);
vYt = vMean + vError;
vYb = flipud(vMean - vError);

% Plot error fill
hFill = fill([vXt;vXb], [vYt;vYb]', vColor);
hold on;
ThemeObject(hFill);
set(hFill, 'facecolor', vColor./2, 'EdgeColor', vColor./2)

% Plot average line
hPlot = plot(vXt, vMean);
ThemeObject(hPlot)
set(hPlot, 'color', vColor, 'LineWidth', 1)
hold off

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = RunAnalysis(hObject, eventdata, sFun, sDir)
% Run selected analysis function from GUI. This function ensures that any
% outputs are captured from the analysis and stored as needed.
% sDir is the analysis function directory, e.g. 'discrete' or 'continuous'
% sFun is the analysis function name (.m not necessary).
%
if ~IsDataLoaded() return, end
[FV, ~] = GetStruct();

% Save batch job
sPath = [sDir filesep() sFun];
[sFullPath, ~] = fileparts(which(sPath));
sEval = [sFullPath filesep() sFun '(FV)'];
sLocalEval = [sFun '(FV)'];
BatchRedo([], sEval)

% Move to directory of analysis
sCurrDir = pwd;
cd(sFullPath)

% Run analysis function
switch nargout(sFun)
    case 0
        eval(sLocalEval);
    case 1
        tSig = eval(sLocalEval);
    otherwise
        warndlg(['Analysis function ' sFun ' returns an invalid number of output arguments.'], 'Spiky')
end
cd(sCurrDir) % return to directory

if ~exist('tSig', 'var'), return; end
if isempty(tSig), return; end

% Validate new signal
cFields = sort(fieldnames(tSig)); % shortest field listed first
if isempty(cell2mat([strfind(cFields, 'KHz')])) ...
        || isempty(cell2mat([strfind(cFields, 'TimeBegin')])) ...
        || isempty(cell2mat([strfind(cFields, 'TimeEnd')]))
    warndlg(['Analysis function ' sFun ' did not return a valid data structure.'], 'Spiky')
    return
end

% Insert new signal into FV.tData and update GUI
for i = 1:length(cFields)
    FV.tData.(cFields{i}) = tSig.(cFields{i});
end
FV.tGain.(cFields{1}) = 1;
FV.csChannels = unique([FV.csChannels cFields{1}]);
FV.csDisplayChannels = unique([FV.csDisplayChannels cFields{1}]);
SetStruct(FV)
ViewTrialData()

return
