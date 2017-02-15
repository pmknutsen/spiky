function FV = spiky_Experiment_Coordinator(FV)
% Experiment coordinator
%
% Coordinates an experimental rig using Open Ephys and Arduino components.
% Assumes the Arduino is fully programmed operating in autonomous mode.
%
% TODO:
%   Start/stop experiment
%   Set Open Ephys filename
%   Monitor serial device in real-time and log to file
%   Send commands to Arduino
%   Automatically load saved Open Ephys files into Spiky
%
%   Button for sending START/STOP commands to arduino
%   Input field to set filename and directory
%
%   Set up alarms/triggers on specific COM events (e.g. turn on a little
%   indicator when COM matches 'REWARD')
%
%   Count count of unique COM command, e.g. how many times has STIMULUS_ON
%   appeared etc.
%
%   See also TODOs below.
%

global Spiky;

%%
%% Initialize GUI
hWin = figure('tag', 'spiky_experimentcoordinator');
set(hWin, 'menuBar', 'none', 'resize', 'on', ...
    'position', [1000 800 800 500], ...
    'name', 'Experiment Coordinator (Open EPhys + Arduino)');
centerfig(hWin);
Spiky.main.ThemeObject(hWin)

hMonPanel = uipanel('position', [0 0 0.4 1]);
hCtrlPanel = uipanel('position', [.4 0 0.6 1]);

% GUI elements
hObj = [];
hText = [];
hFields = [];

% List of serial ports
tCOMInfo = instrhwinfo('serial');
hObj(end+1) = uicontrol(hWin, 'Style', 'popupmenu', 'units', 'normalized', ...
    'Position', [.5 .95 .5 .05], 'String', tCOMInfo.SerialPorts, ...
    'callback', @SelectCOM, 'parent', hMonPanel);
hText(end+1) = uicontrol(hWin, 'Style', 'edit', 'units', 'normalized', ...
    'Position', [0 0.95 0.5 .05], 'string', 'Arduino COM', ...
    'foregroundcolor', 'r', 'tag', 'expcord_COMstatus', 'parent', hMonPanel);

% Serial monitor
hObj(end+1) = uicontrol(hWin, 'Style', 'edit', 'units', 'normalized', ...
    'enable', 'inactive', 'max', 100, 'parent', hMonPanel, ...
    'Position', [0 .35 1 .6], 'tag', 'expcord_COMmonitor', ...
    'horizontalalignment', 'left', 'fontsize', 8, 'fontname', 'courier');

% Serial monitor counters
hText(end+1) = uicontrol(hWin, 'Style', 'edit', 'parent', hMonPanel, ...
    'units', 'normalized', 'Position', [0 .3 1 .05], ...
    'string', 'Counters');
hObj(end+2) = uicontrol(hWin, 'Style', 'edit', 'units', 'normalized', ...
    'enable', 'inactive', 'max', 100, 'parent', hMonPanel, ...
    'Position', [0 .05 1 .25], 'tag', 'expcord_COMcounter', ...
    'horizontalalignment', 'left', 'fontsize', 8, 'fontname', 'courier');

% Trigger
hText(end+1) = uicontrol(hWin, 'Style', 'edit', 'parent', hMonPanel, ...
    'units', 'normalized', 'Position', [0 0 .5 .05], ...
    'string', 'SW Trigger 1');
hObj(end+1) = uicontrol(hWin, 'Style', 'popupmenu', 'units', 'normalized', ...
    'Position', [.5 0 .5 .05], 'tag', 'expcord_SWtrigger', 'string', 'None', ...
    'callback', @SetGetSWTrigger, 'parent', hMonPanel);

% Input fields
uicontrol(hWin, 'Style', 'edit', 'units', 'normalized', ...
    'Position', [.4 0.95 .6 .05], 'tag', '', 'parent', hCtrlPanel);

% Input fields
tDefs = InputDefaults();
for d = 1:length(tDefs)
    % Description
    hText(end+1) = uicontrol(hWin, 'Style', 'edit', 'parent', hCtrlPanel, ...
        'units', 'normalized', 'Position', [0 1-(d*0.05) .4 .05], ...
        'string', tDefs(d).name);
    % Input field
    hFields(end+1) = uicontrol(hWin, 'Style', 'edit', 'parent', hCtrlPanel, ...
        'units', 'normalized', 'Position', [.4 1-(d*0.05) .6 .05], ...
        'string', tDefs(d).default, ...
        'tag', sprintf('expcord_%s', regexprep(lower(tDefs(d).name), '[^\w'']', '')));
end

% Start/Stop button
uicontrol(hWin, 'Style', 'togglebutton', 'units', 'normalized', ...
    'backgroundcolor', 'r', 'string', 'START', ...
    'Position', [.4 0 .1 .1], 'callback', @StartStopExp);

Spiky.main.ThemeObject(hObj)
Spiky.main.ThemeObject(hText, 'fontsize', 8, 'enable', 'inactive', 'fontweight', 'bold')

%%
return

%% Start/Stop experiment
function StartStopExp(hBut, varargin)
% TODO
%   
% Get current status
hCOM = SetGetCOM();
if get(hBut, 'Value')
    % Start experiment
    % Clear COM log
    set(findobj('tag', 'expcord_COMmonitor'), 'string', '');
    set(findobj('tag', 'expcord_COMcounter'), 'string', '');

    % Send SESSION_INIT command Arduino
    if isobject(hCOM)
        fwrite(hCOM, 4)
    else
        SetCOMStatus(0);
    end
    
    % Update button status
    set(hBut, 'backgroundcolor', 'g', 'string', 'STOP')
else
    % Stop experiment

    % Send SESSION_WAIT command Arduino
    if isobject(hCOM)
        fwrite(hCOM, 5)
    else
        SetCOMStatus(0);
    end
    
    % Update button status
    set(hBut, 'backgroundcolor', 'r', 'string', 'START')
end

return

%% Get figure handle
function hWin = GetWin(varargin)
hWin = findobj('tag', 'spiky_experimentcoordinator');
hWin = hWin(end);
return

%% Open connection to COM port
function SelectCOM(varargin)
global Spiky;
hCOM = SetGetCOM();
% Close existing COM connections
if ~isempty(hCOM)
    delete(instrfindall) % if all fails, this cmd closes ALL connections
    Spiky.main.sp_disp('Close all COM connections')
    SetCOMStatus(0);
end
hCOMui = varargin{1};
hCOM = SetGetCOM(serial(hCOMui.String(hCOMui.Value), 'BaudRate', 9600));
try
    fopen(hCOM)
    Spiky.main.sp_disp(sprintf('Connected to %s', hCOMui.String{hCOMui.Value}))
    set(hCOM, 'BytesAvailableFcn', {@UpdateCOMMonitor, hCOM})
    SetCOMStatus(1);
catch
    Spiky.main.sp_disp(sprintf('Failed connecting to %s', hCOMui.String{hCOMui.Value}))
    SetCOMStatus(0);
end
return

%% Set the COM connection status indicator
function SetCOMStatus(bStatus)
hInd = findobj('tag', 'expcord_COMstatus');
switch bStatus
    case 1
        set(hInd, 'foregroundcolor', 'g')
    otherwise
        set(hInd, 'foregroundcolor', 'r')
end
return

%% Store and retrieve COM handle as persistent variable
function hCOM = SetGetCOM(varargin)
persistent p_hCOM
if nargin == 1
    p_hCOM = varargin{1};
end
hCOM = p_hCOM;
return

%% Update the serial port monitor
function UpdateCOMMonitor(varargin)
% TODO  Print to log
%       Limit display to 100 lines
hCOM = varargin{end};
sNewStr = regexprep(strtrim(fscanf(hCOM)), '\t', ' ');
hMon = findobj('tag', 'expcord_COMmonitor');
csCurStr = flipud(get(hMon, 'string'));
if isempty(csCurStr)
    csUpdStr{1} = sNewStr;
else
    csUpdStr = csCurStr;
    csUpdStr{end+1} = sNewStr;
end
set(hMon, 'string', flipud(csUpdStr))

% Update counters
UpdateCOMCounters(sNewStr);
return

%% Set or get software trigger string
function sSWTrigger = SetGetSWTrigger(varargin)
persistent p_sSWTrigger
if nargin > 0
    if isobject(varargin{1})
        nVal = get(varargin{1}, 'value');
        csStr = get(varargin{1}, 'string');
        p_sSWTrigger = csStr{nVal};
    else
        p_sSWTrigger = varargin{1};
    end
end
sSWTrigger = p_sSWTrigger;
return

%% Update COM counters
function UpdateCOMCounters(sStr)
persistent tCounters
if isempty(tCounters)
    tCounters = struct();
end

% Parse the new string
[cTokens, ~] = regexp(sStr, '(\d+).*\[(.*)\]', 'tokens', 'match');
if isempty(cTokens), return; end
nDms = cTokens{1}{1};
sCmd = cTokens{1}{2};

% Replace whitespaces with underscore
sCmd = regexprep(sCmd, '[^\w'']', '_');
if length(sCmd) < 2, return; end
nDms = str2num(nDms);

if isfield(tCounters, sCmd)
    tCounters(1).(sCmd)(1).nCount = tCounters(1).(sCmd).nCount + 1;
    tCounters(1).(sCmd)(1).vDms = [tCounters(1).(sCmd)(1).vDms nDms];
else
    tCounters(1).(sCmd).nCount = 1;
    tCounters(1).(sCmd).vDms = nDms;
end

% Construct text that goes into the counter monitor
csFields = fieldnames(tCounters);
csStr = {};
csStrSW = {};
for f = 1:length(csFields)
    csStr{end+1} = sprintf('%.3d  %s', tCounters.(csFields{f}).nCount, csFields{f});
    csStrSW{end+1} = csFields{f};
end
set(findobj('tag', 'expcord_COMcounter'), 'string', csStr)

% Update SW trigger selector
csStrSW = ['None' csStrSW];
hMenu = findobj('tag', 'expcord_SWtrigger');
nVal = get(hMenu, 'value');
if nVal > length(csStrSW)
    nVal = length(csStrSW);
end
set(hMenu, 'value', nVal, 'string', csStrSW);

return

%% Default input fields and values
function tDefs = InputDefaults

tDefs = struct();
tDefs(1).name = 'Directory';
tDefs(1).default = 'C:\Data\';

tDefs(2).name = 'Filename';
tDefs(2).default = 'testfile';

tDefs(3).name = 'File prefix';
tDefs(3).default = 'exp_';

tDefs(4).name = 'Open Ephys IP';
tDefs(4).default = 'localhost';

return
