function FV = spiky_Experiment_Coordinator(FV)
% Experiment coordinator
%
% Coordinates an experimental rig using Open Ephys and Arduino components.
% Assumes the Arduino is fully programmed to run in autonomous mode.
%
% TODO:
%   Log arduino commands to file
%

% Set the number of software defined triggers
SetGetSWTriggerNum(2);

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
Spiky.main.ThemeObject([hMonPanel hCtrlPanel]);

% GUI elements
hObj = [];
hText = [];
hFields = [];
hButs = [];

% List of serial ports
tCOMInfo = instrhwinfo('serial');
hText(end+1) = uicontrol(hWin, 'Style', 'edit', 'units', 'normalized', ...
    'Position', [0 0.95 0.3 .05], 'string', 'COM Log', ...
    'foregroundcolor', 'r', 'parent', hMonPanel);
hObj(end+1) = uicontrol(hWin, 'Style', 'popupmenu', 'units', 'normalized', ...
    'Position', [.3 .95 .4 .05], 'String', tCOMInfo.SerialPorts, ...
    'tag', 'expcord_COMselection', 'parent', hMonPanel);
hButs(end+1) = uicontrol(hWin, 'Style', 'pushbutton', 'units', 'normalized', ...
    'string', 'Clear', 'callback', @ClearCOMLog, ...
    'parent', hMonPanel, 'Position', [0.7 .95 .3 .05]);

% Serial monitor
hObj(end+1) = uicontrol(hWin, 'Style', 'edit', 'units', 'normalized', ...
    'enable', 'inactive', 'max', 100, 'parent', hMonPanel, ...
    'Position', [0 .35 1 .6], 'tag', 'expcord_COMmonitor', ...
    'horizontalalignment', 'left', 'fontsize', 8, 'fontname', 'courier');

% Serial monitor counters
hText(end+1) = uicontrol(hWin, 'Style', 'edit', 'parent', hMonPanel, ...
    'units', 'normalized', 'Position', [0 .3 .7 .05], 'string', 'Counters');
hButs(end+1) = uicontrol(hWin, 'Style', 'pushbutton', 'units', 'normalized', ...
    'string', 'Reset', 'callback', {@GetSetCounters, struct()}, 'parent', hMonPanel, ...
    'Position', [0.7 .3 .3 .05]);
hObj(end+2) = uicontrol(hWin, 'Style', 'edit', 'units', 'normalized', ...
    'enable', 'inactive', 'max', 100, 'parent', hMonPanel, ...
    'Position', [0 .05 1 .25], 'tag', 'expcord_COMcounter', ...
    'horizontalalignment', 'left', 'fontsize', 8, 'fontname', 'courier');

% Triggers
vT = fliplr(1:SetGetSWTriggerNum());
for t = vT
    hText(end+1) = uicontrol(hWin, 'Style', 'edit', 'parent', hMonPanel, ...
        'units', 'normalized', 'Position', [0 .05*(t-1) .5 .05], ...
        'string', sprintf('SW Trigger %d', vT(t)));
    hObj(end+1) = uicontrol(hWin, 'Style', 'popupmenu', 'units', 'normalized', ...
        'Position', [.5 .05*(t-1) .5 .05], 'tag', sprintf('expcord_SWtrigger_%d', vT(t)), ...
        'string', 'None', 'callback', @SetGetSWTrigger, 'parent', hMonPanel);
end

% Input fields
uicontrol(hWin, 'Style', 'edit', 'units', 'normalized', ...
    'Position', [.4 0.95 .6 .05], 'tag', '', 'parent', hCtrlPanel);

% Input fields
tDefs = InputDefaults();
nSpacer = 0;
nD = 0;
for d = 1:length(tDefs)
    if isempty(tDefs(d).name)
        nSpacer = nSpacer + .01;
        continue;
    end
    nD = nD + 1;
    % Description
    hText(end+1) = uicontrol(hWin, 'Style', 'edit', 'parent', hCtrlPanel, ...
        'units', 'normalized', 'Position', [0 1-(nD*0.05)-nSpacer .4 .05], ...
        'string', tDefs(d).name);
    % Input field
    if isempty(tDefs(d).tooltip), tDefs(d).tooltip = ''; end
    hFields(end+1) = uicontrol(hWin, 'Style', 'edit', 'parent', hCtrlPanel, ...
        'units', 'normalized', 'Position', [.4 1-(nD*0.05)-nSpacer .6 .05], ...
        'string', tDefs(d).default, 'tooltipString', tDefs(d).tooltip, ...
        'tag', sprintf('expcord_%s', TagifyString(tDefs(d).name)));
end

% Save/load settings buttons
hButs = [];
hButs(end+1) = uicontrol(hWin, 'Style', 'pushbutton', 'units', 'normalized', ...
    'string', 'Load Defaults', 'callback', {@LoadSettings, 'default'}, ...
    'parent', hCtrlPanel, 'Position', [0 1-((nD+1)*0.05)-nSpacer .25 .05]);
hButs(end+1) = uicontrol(hWin, 'Style', 'pushbutton', 'units', 'normalized', ...
    'string', 'Save As Default', 'callback', {@SaveSettings, 'default'}, ...
    'parent', hCtrlPanel, 'Position', [.25 1-((nD+1)*0.05)-nSpacer .25 .05]);
hButs(end+1) = uicontrol(hWin, 'Style', 'pushbutton', 'units', 'normalized', ...
    'string', 'Load Settings', 'callback', @LoadSettings, ...
    'parent', hCtrlPanel, 'Position', [.75 1-((nD+1)*0.05)-nSpacer .25 .05]);
hButs(end+1) = uicontrol(hWin, 'Style', 'pushbutton', 'units', 'normalized', ...
    'string', 'Save Settings', 'callback', @SaveSettings, ...
    'parent', hCtrlPanel, 'Position', [.5 1-((nD+1)*0.05)-nSpacer .25 .05]);

% Status buttons
tDefs = struct();
csButs = {'Arduino COM', 'COM Running', 'OE Acquiring', 'OE Recording', 'Save'};
for b = 1:length(csButs)
    hButs(end+1) = uicontrol(hWin, 'Style', 'togglebutton', 'units', 'normalized', ...
        'backgroundcolor', 'r', 'string', csButs{b}, ...
        'tag', sprintf('expcord_%s', TagifyString(csButs{b})), ...
        'horizontalalignment', 'center', 'callback', @ButtonAction, ...
        'parent', hCtrlPanel, 'Position', [(b*.2)-.2 0 .2 .1]);
end
Spiky.main.ThemeObject([hObj hFields])
Spiky.main.ThemeObject(hText, 'fontsize', 8, 'enable', 'inactive', 'fontweight', 'bold')
LoadSettings('default');
UpdateToggleButtons();
return

% Get the number of defined software triggers
function nNumTrig = SetGetSWTriggerNum(varargin)
global g_expCord_SWTrigNum
if nargin > 0
    g_expCord_SWTrigNum = varargin{1};
end
nNumTrig = g_expCord_SWTrigNum;
return

% Save settings (only edit fields are saved)
function SaveSettings(varargin)
if strcmpi(varargin{end}, 'default')
    [sDir, ~] = fileparts(which(mfilename));
    sFile = fullfile(sDir, 'expcord_default_settings.mat');
end
hEdit = findobj(GetWin(), 'style', 'edit');
if exist('sFile')
    save(sFile, 'hEdit');
else
    uisave({'hEdit'}, 'expcoord_settings');
end
return

% Load settings
function LoadSettings(varargin)
if strcmpi(varargin{end}, 'default')
    [sDir, ~] = fileparts(which(mfilename));
    sFile = fullfile(sDir, 'expcord_default_settings.mat');
end
if ~exist('sFile')
    [sFile, sPath] = uigetfile({'*.mat';'*.*'}, 'Load settings');
    if ~sFile, return; end
    sFile = fullfile(sPath, sFile);
end
if ~exist(sFile, 'file'), return; end
load(sFile, 'hEdit')
for h = 1:length(hEdit)
    hObj = findobj(GetWin(), 'tag', hEdit(h).Tag);
    if isprop(hObj, 'Max')
        if hObj.Max == 1
            hObj.String = hEdit(h).String;
        end
    end
end
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
    case 'COM Running'
        StartStopArduinoProgram();
    case 'OE Acquiring'
        ToggleOEAcquisition();
    case 'OE Recording'
        ToggleOERecording();
    case 'Save'
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
        sRecCmd = [sRecCmd ' AppendText=' sAppendText];
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

% Start/stop the loaded program on the connected Arduino
function StartStopArduinoProgram(varargin)
hCOM = SetGetCOM();
persistent p_bRunning
if ~isobject(hCOM)
    SetCOMStatus(0);
    bStart = 0;
else
    % Toggle COM Running status when no input is provided
    if nargin == 1
        bStart = varargin{1};
    else
        if p_bRunning
            bStart = 0;
        else
            bStart = 1;
        end
    end
    if bStart
        hObj = findobj(GetWin(), 'tag', TagifyString('expcord_Arduino start cmd'));
        fwrite(hCOM, str2num(get(hObj, 'string')))
        p_bRunning = 1;
    else
        hObj = findobj(GetWin(), 'tag', TagifyString('expcord_Arduino stop cmd'));
        fwrite(hCOM, str2num(get(hObj, 'string')))
        p_bRunning = 0;
    end
end
% Set button status
hBut = findobj(GetWin(), 'tag', TagifyString(sprintf('expcord_COM Running')));
SetButtonStatus(hBut, bStart)
return

% Start/stop experiment
function StartStopExp(hBut)
if isobject(hBut)
    bStart = get(hBut, 'Value');
elseif isnumeric(hBut)
    bStart = hBut;
    hBut = findobj(GetWin(), 'string', 'Save', 'style', 'togglebutton');
end
if bStart
    % Clear COM log
    set(findobj('tag', 'expcord_COMmonitor'), 'string', '');
    set(findobj('tag', 'expcord_COMcounter'), 'string', '');
    
    % Start Open Ephys acquisition
    ToggleOEAcquisition(1);
    
    % Start Open Ephys recording
    ToggleOERecording(1);

    % Start program on COM
    StartStopArduinoProgram(1);
    
else
    % Stop program on COM
    StartStopArduinoProgram(0);

    % Stop Open Ephys recording
    ToggleOEAcquisition(0);
end
SetButtonStatus(hBut, bStart)
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

% Close all COM connection
function DisconnectCOM()
global Spiky
% Try to stop running program before disconnecting
try
    StartStopArduinoProgram(0);
end
delete(instrfindall) % if all fails, this cmd closes ALL connections
Spiky.main.sp_disp('Close all COM connections')
SetCOMStatus(0);
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

% Update COM counters
function UpdateCOMCounters(sStr)
tCounters = GetSetCounters();

% Parse the new string
[cTokens, ~] = regexp(sStr, '(\d+).*\[(.*)\]', 'tokens', 'match');
if isempty(cTokens), return; end
nDms = cTokens{1}{1};
sCmd = cTokens{1}{2};

% Replace whitespaces with underscore
sCmd = regexprep(sCmd, '[^\w'']', '_');
if length(sCmd) < 2, return; end
nDms = str2num(nDms);

if all(ismember(sCmd, '0123456789+-.eEdD'))
    sCmd = ['n' sCmd];
end

if isfield(tCounters, sCmd)
    tCounters(1).(sCmd)(1).nCount = tCounters(1).(sCmd).nCount + 1;
    tCounters(1).(sCmd)(1).vDms = [tCounters(1).(sCmd)(1).vDms nDms];
else
    tCounters(1).(sCmd).nCount = 1;
    tCounters(1).(sCmd).vDms = nDms;
end

GetSetCounters(tCounters); % save counters
% Do not change counters after this point!

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

% Check if Stop Saving counter should be triggered
hStop = findobj(GetWin(), 'tag', TagifyString(sprintf('expcord_Stop saving when Counter')));
sStop = get(hStop, 'string');
if strcmpi(sCmd, sStop)
    hCnt = findobj(GetWin(), 'tag', TagifyString(sprintf('expcord_exceeds')));
    nCnt = str2num(get(hCnt, 'string'));
    if tCounters.(sCmd).nCount > nCnt
        StartStopExp(0);
    end
end

% Update SW trigger selector
csStrSW = ['None' csStrSW];
for t = 1:2
    hMenu = findobj('tag', sprintf('expcord_SWtrigger_%d', t));
    nVal = get(hMenu, 'value');
    if nVal > length(csStrSW)
        nVal = length(csStrSW);
    end
    set(hMenu, 'value', nVal, 'string', csStrSW);
end
return

% Set or get current counter structure
function tCounters = GetSetCounters(varargin)
persistent p_tCounters
if nargin > 0
    if isstruct(varargin{end})
        p_tCounters = varargin{end};
    end
end
if isempty(p_tCounters)
    p_tCounters = struct();
end
if isempty(fieldnames(p_tCounters))
    set(findobj('tag', 'expcord_COMcounter'), 'string', '')
end
tCounters = p_tCounters;
return

% Set or get software trigger string
function sSWTrigger = SetGetSWTrigger(varargin)
sSWTrigger = ''; % default
if nargin > 1
    if isobject(varargin{1})
        nVal = get(varargin{1}, 'value');
        csStr = get(varargin{1}, 'string');
        sTrigCh = get(varargin{1}, 'tag');
        if iscell(csStr)
            sStr = csStr{nVal};
        else
            sStr = csStr;
        end
        set(findobj(GetWin(), 'tag', sprintf('expcord_swtrigger%s', sTrigCh(end))), ...
            'string', sStr);
    end
elseif nargin == 1
    if isnumeric(varargin{1})
        sSWTrigger = get(findobj(GetWin(), 'tag', ...
            sprintf('expcord_swtrigger%d', varargin{1})), ...
            'string');
    end
end
return

% Check software trigger
function CheckSWTrigger(sCmd)
for t = 1:SetGetSWTriggerNum()
    sTrig = SetGetSWTrigger(t);
    if strcmp(sTrig, sCmd)
        RunSWTrigger(t);
    end
end
return

% Execute events that should run when a software trigger occurs
function RunSWTrigger(nTrig)
global Spiky
% Send custom command to Open Ephys
sTag = TagifyString(sprintf('expcord_OE SW Trigger %d Cmd', nTrig));
hInp = findobj('tag', sTag);
sCmd = get(hInp, 'string');
if ~isempty(sCmd)
    set(hInp, 'backgroundcolor', 'g')
    SendOE(sCmd);
    pause(.05)
    Spiky.main.ThemeObject(hInp)
end
return

% Default input fields and values
function tDefs = InputDefaults()
tDefs = struct();
tDefs(1).name = 'Recording Directory';
tDefs(end).default = 'C:\Data\';
tDefs(end).tooltip = 'Directory where data files will be saved';

tDefs(end+1).name = 'Filename prefix';
tDefs(end).default = 'oe_';
tDefs(end).tooltip = 'String prepended to filenames';

tDefs(end+1).name = 'Filename postfix';
tDefs(end).default = '';
tDefs(end).tooltip = 'String appended to filenames';

tDefs(end+1).name = 'Open Ephys IP:Port';
tDefs(end).default = 'localhost:5556';
tDefs(end).tooltip = 'Hostname or IP and port of machine running the Open Ephys GUI';

tDefs(end+1).name = ''; % spacer

for t = 1:SetGetSWTriggerNum()
    tDefs(end+1).name = sprintf('SW Trigger %d', t);
    tDefs(end).default = '';
    tDefs(end).tooltip = sprintf('Default Arduino software trigger command %d', t);
end

for t = 1:SetGetSWTriggerNum()
    tDefs(end+1).name = sprintf('OE SW Trigger %d Cmd', t);
    tDefs(end).default = '';
    tDefs(end).tooltip = 'Enter a valid Open Ephys network command';
end

tDefs(end+1).name = ''; % spacer

tDefs(end+1).name = 'Start saving on SW Trigger';
tDefs(end).default = '';
tDefs(end).tooltip = sprintf('Choose a software trigger number between 1 and %d', SetGetSWTriggerNum());

tDefs(end+1).name = 'Stop saving on SW Trigger';
tDefs(end).default = '';
tDefs(end).tooltip = sprintf('Choose a software trigger number between 1 and %d', SetGetSWTriggerNum());

tDefs(end+1).name = ''; % spacer

tDefs(end+1).name = 'Arduino start cmd'; % 4, SESSION_INIT
tDefs(end).default = '4';
tDefs(end+1).name = 'Arduino stop cmd';
tDefs(end).default = '5'; % 5, SESSION_WAIT

tDefs(end+1).name = ''; % spacer

tDefs(end+1).name = 'Stop saving when Counter';
tDefs(end).default = '';
tDefs(end).tooltip = 'Choose a counter which will stop saving';

tDefs(end+1).name = '... exceeds';
tDefs(end).default = '';
tDefs(end).tooltip = 'Enter a number how many counts will be recorded before Save is stopped';

tDefs(end+1).name = ''; % spacer

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

% Clear COM log
function ClearCOMLog(varargin)
set(findobj(GetWin(), 'tag', 'expcord_COMmonitor'), 'string', '')
return
