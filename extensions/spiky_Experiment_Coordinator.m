function FV = spiky_Experiment_Coordinator(FV)
% Experiment coordinator
%
% Coordinates an experimental rig using Open Ephys and Arduino components.
% Assumes the Arduino is fully programmed operating in autonomous mode.
%
% TODO:
%   Set Open Ephys filename
%   Monitor serial device in real-time and log to file
%   Send commands to Arduino
%   Automatically load saved Open Ephys files into Spiky
%   Connect to Open Ephys
%

InitGUI();
return
%%

% Initialize GUI
function InitGUI()
global Spiky

%% Bring existing GUI to front
hWin = findobj('tag', 'spiky_experimentcoordinator');
if ~isempty(hWin)
    figure(hWin)
    return;
end

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
    'tag', 'expcord_COMselection', 'parent', hMonPanel);
hText(end+1) = uicontrol(hWin, 'Style', 'edit', 'units', 'normalized', ...
    'Position', [0 0.95 0.5 .05], 'string', 'Arduino COM', ...
    'foregroundcolor', 'r', 'parent', hMonPanel);

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
    'string', 'SW Trigger');
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
    if isempty(tDefs(d).tooltip), tDefs(d).tooltip = ''; end
    hFields(end+1) = uicontrol(hWin, 'Style', 'edit', 'parent', hCtrlPanel, ...
        'units', 'normalized', 'Position', [.4 1-(d*0.05) .6 .05], ...
        'string', tDefs(d).default, 'tooltipString', tDefs(d).tooltip, ...
        'tag', sprintf('expcord_%s', TagifyString(tDefs(d).name)));
end

% Status buttons
tDefs = struct();
csButs = {'Arduino COM', 'OE Acquiring', 'OE Recording', 'SW Trigger', 'Running'};
hButs = [];
for b = 1:length(csButs)
    hButs(end+1) = uicontrol(hWin, 'Style', 'togglebutton', 'units', 'normalized', ...
        'backgroundcolor', 'r', 'string', csButs{b}, ...
        'tag', sprintf('expcord_%s', TagifyString(csButs{b})), ...
        'horizontalalignment', 'center', 'callback', @ButtonAction, ...
        'parent', hCtrlPanel, 'Position', [(b*.2)-.2 0 .2 .1]);
end
Spiky.main.ThemeObject(hObj)
Spiky.main.ThemeObject(hText, 'fontsize', 8, 'enable', 'inactive', 'fontweight', 'bold')
return

% Run an action on UI button press
function ButtonAction(hObj, varargin)
if ~isobject(hObj), return; end
sStr = get(hObj, 'string');
nVal = get(hObj, 'value');
switch sStr
    case 'Arduino COM'
        if nVal == 1
            ConnectCOM();
        else
            DisconnectCOM();
        end
    case 'OE Acquiring'
        ToggleOEAcquisition();
    case 'OE Recording'
        ToggleOERecording();
    case 'SW Trigger'
        RunSWTrigger();
    case 'Running'
        StartStopExp(hObj);
    otherwise
        % undefined button
end
return

% Set status of a toggle button
function SetButtonStatus(hBtn, bStatus)
switch bStatus
    case 1
        set(hBtn, 'backgroundcolor', 'g', 'value', 1)
    otherwise
        set(hBtn, 'backgroundcolor', 'r', 'value', 0)
end
return

% Send TCP command to Open Ephys GUI
function [sResp, tInfo] = SendOE(sCmd)
sURL = ['tcp://' get(findobj('tag', 'expcord_openephysipport'), 'string')];
[sResp, tInfo] = zeroMQrr('Send', sURL, sCmd, 1);
% Check if connection is good
if strcmp(sResp, 'failed waiting for a reply')
    sResp = 'nan';
end
return

% Check if Open Ephys is acquiring data
function IsOEAcquiring()
[sResp, ~] = SendOE('IsAcquiring');
SetButtonStatus(findobj('tag', 'expcord_oeacquiring'), str2num(sResp))
return

% Check if Open Ephys is recording data
function IsOERecording()
[sResp, ~] = SendOE('IsRecording');
SetButtonStatus(findobj('tag', 'expcord_oerecording'), str2num(sResp))
return

% Start Open Ephys data acquisition
function ToggleOEAcquisition(varargin)
if nargin == 1
    nVal = varargin{1};
else
    nVal = get(findobj('tag', 'expcord_oeacquiring'), 'value');
end
if nVal == 1
    SendOE('StartAcquisition');
else
    SendOE('StopAcquisition');
end
UpdateToggleButtons();
return

% Start Open Ephys data acquisition
function ToggleOERecording(varargin)
if nargin == 1
    nVal = varargin{1};
else
    nVal = get(findobj('tag', 'expcord_oerecording'), 'value');
end
if nVal == 1
    % Build filename
    sRecCmd = 'StartRecord';
    
    sRecDir = get(findobj('tag', 'expcord_recordingdirectory'), 'string');
    if ~isempty(sRecDir)
        sRecCmd = [sRecCmd ' RecDir=' sRecDir];
    end
    
    sPrependText = get(findobj('tag', 'expcord_filenameprefix'), 'string');
    if ~isempty(sPrependText)
        sRecCmd = [sRecCmd ' PrependText=' sPrependText];
    end
    
    sAppendText = get(findobj('tag', 'expcord_filenamepostfix'), 'string');
    if ~isempty(sAppendText)
        sRecCmd = [sRecCmd ' PrependText=' sAppendText];
    end
    SendOE(sRecCmd);
else
    SendOE('StopRecord');
end
UpdateToggleButtons();
return

% Update the 
function UpdateToggleButtons()
IsOEAcquiring();
IsOERecording();
return

% Start/Stop experiment
function StartStopExp(hBut)
hCOM = SetGetCOM();
if get(hBut, 'Value')
    % Start experiments
    ToggleOEAcquisition(1); % start OE acquisition

    % Clear COM log
    set(findobj('tag', 'expcord_COMmonitor'), 'string', '');
    set(findobj('tag', 'expcord_COMcounter'), 'string', '');
    if isobject(hCOM)
        fwrite(hCOM, 4) % SESSION_INIT
    else
        SetCOMStatus(0);
    end
    SetButtonStatus(hBut, 1)
else
    % Stop experiment
    if isobject(hCOM)
        fwrite(hCOM, 5) % SESSION_WAIT
    else
        SetCOMStatus(0);
    end
    SetButtonStatus(hBut, 0)
    ToggleOEAcquisition(0); % stop OE acquisition
end
return

% Close all COM connection
function DisconnectCOM()
global Spiky
delete(instrfindall) % if all fails, this cmd closes ALL connections
Spiky.main.sp_disp('Close all COM connections')
SetCOMStatus(0);
return

% Open connection to COM port
function ConnectCOM()
global Spiky;
hCOM = SetGetCOM();
% Close existing COM connections
if ~isempty(hCOM)
    DisconnectCOM()
end
hCOMui = findobj('tag', 'expcord_COMselection');
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

% Set the COM connection status indicator
function SetCOMStatus(bStatus)
hInd = findobj('tag', 'expcord_arduinocom');
SetButtonStatus(hInd, bStatus)
return

% Store and retrieve COM handle as persistent variable
function hCOM = SetGetCOM(varargin)
persistent p_hCOM
if nargin == 1
    p_hCOM = varargin{1};
end
hCOM = p_hCOM;
return

% Update the serial port monitor
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

% Set or get software trigger string
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

% Check software trigger
function CheckSWTrigger(sCmd)
sTrig = SetGetSWTrigger();
if strcmp(sTrig, sCmd)
    RunSWTrigger();
end
return

% Execute events that should run when a software trigger occurs
function RunSWTrigger()
hTrig = findobj('tag', 'expcord_swtrigger');
SetButtonStatus(hTrig, 1);
% Send custom command to Open Ephys
sCmd = get(findobj('tag', TagifyString('expcord_Open Ephys SW Trigger Cmd')), 'string');
if ~isempty(sCmd)
    SendOE(sCmd);
end
pause(0.1)
SetButtonStatus(hTrig, 0)
return

% Update COM counters
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

% Check software trigger(s)
CheckSWTrigger(sCmd);

% Update SW trigger selector
csStrSW = ['None' csStrSW];
hMenu = findobj('tag', 'expcord_SWtrigger');
nVal = get(hMenu, 'value');
if nVal > length(csStrSW)
    nVal = length(csStrSW);
end
set(hMenu, 'value', nVal, 'string', csStrSW);
return

% Default input fields and values
function tDefs = InputDefaults()
tDefs = struct();
tDefs(1).name = 'Recording Directory';
tDefs(1).default = 'C:\Data\';

tDefs(end+1).name = 'Filename prefix';
tDefs(end).default = 'oe_';
tDefs(end).tooltip = 'Prepended to filename';

tDefs(end+1).name = 'Filename postfix';
tDefs(end).default = '';
tDefs(end).tooltip = 'Appended to filename';

tDefs(end+1).name = 'Open Ephys IP:Port';
tDefs(end).default = 'localhost:5556';

tDefs(end+1).name = 'Open Ephys SW Trigger Cmd';
tDefs(end).default = '';
tDefs(end).tooltip = 'StartAcquisition, StopAcquisition\nStartRecord, StopRecord, IsAcquiring, IsRecording, GetRecordingPath, GetRecordingNumber, GetExperimentNumber';
return

% Get figure handle
function hWin = GetWin(varargin)
hWin = findobj('tag', 'spiky_experimentcoordinator');
hWin = hWin(end);
return

% Turn a string into a tag by removing all whitespaces and converting all characters to lowercase.
function sStr = TagifyString(sStr)
sStr = lower(regexprep(sStr, '[^\w'']', ''));
return

