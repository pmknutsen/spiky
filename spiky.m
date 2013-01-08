function varargout = spiky(varargin)
% SPIKY! AlphaMap .map/.mat browser and spike-sorter
%
% Usage:
%   < spiky >
%   Runs the Spiky! GUI
%
%   < spiky(FUN) >
%   Runs the internal function FUN in Spiky!. Ex:
%     spiky, e.g. spiky('LoadTrial(''A1807004.mat'')')
%
%   < spiky(FUN) >
%   Runs the internal function FUN in Spiky!. Ex:
%     spiky, e.g. spiky('LoadTrial(''A1807004.mat'')')
%

% Set path to spike-sorting software
sPath = which('spiky');
sPath = sPath(1:end-7);
addpath(genpath(sPath), '-end');

% Run local function if called from outside
if nargin > 0
    if strcmp(varargin{1}, 'GetBuildNumber')
        varargout{1} = eval(varargin{1});
    else
        %eval(varargin{1}); % changed 07/11/11
        try
            varargout{1} = eval(varargin{1});
        catch
            eval(varargin{1});
        end
    end
    return
end

% Set FV default settings
FV = SetFVDefaults();

% Check if Spiky! is already running
hPlot = findobj('Tag', 'Spiky');
if ~isempty(hPlot) % activate existing Spiky! window
    figure(hPlot)
    return
end

global g_hSpike_Viewer
g_hSpike_Viewer = figure;
vPos = get(g_hSpike_Viewer, 'position');
set(g_hSpike_Viewer, 'NumberTitle','off', 'Name', ['Spiky! Build# ' GetBuildNumber], 'MenuBar','none', 'UserData', FV, ...
    'Tag', 'Spiky', 'Color', [.2 .2 .2], 'PaperOrientation', 'landscape', 'PaperUnits', 'normalized', ...
    'PaperPosition', [.05 .05 .9 .9], 'InvertHardcopy', 'off', 'position', [vPos(1)-200 vPos(2)-50 vPos(3)+350 vPos(4)+100], ...
    'closeRequestFcn', @ExitSpiky, 'visible', 'off')
movegui(g_hSpike_Viewer, 'center')
set(g_hSpike_Viewer, 'visible', 'on')

% Create toolbar
sPath = which('spiky');
sPath = [sPath(1:end-7) 'icons/'];
hToolbar = uitoolbar('Parent', g_hSpike_Viewer, 'Tag', 'Spiky_Waitbar');

mCData = im2double(imread([sPath 'tool_open.png'])); mCData(mCData == 0) = NaN; % open
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_Open', 'TooltipString', 'Open file', 'ClickedCallback', @OpenFile);
mCData = im2double(imread([sPath 'tool_save.png'])); mCData(mCData == 0) = NaN; % save
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_Save', 'TooltipString', 'Save', 'ClickedCallback', @SaveResults);
mCData = im2double(imread([sPath 'tool_zoom_in.png'])); mCData(mCData == 0) = NaN; % zoom in
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_ZoomOut', 'TooltipString', 'Zoom in', 'ClickedCallback', @ZoomIn, 'separator', 'on');
mCData = im2double(imread([sPath 'tool_zoom_out.png'])); mCData(mCData == 0) = NaN; % zoom out
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_ZoomOut', 'TooltipString', 'Zoom out', 'ClickedCallback', @ZoomOut);
mCData = im2double(imread([sPath 'tool_zoom_reset.png'])); mCData(mCData == 1) = NaN; % zoom reset
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_ZoomReset', 'TooltipString', 'Zoom reset', 'ClickedCallback', @ZoomReset);
mCData = im2double(imread([sPath 'tool_hand.png'])); mCData(mCData == 0) = NaN; % pan tool
uitoggletool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_Pan', 'TooltipString', 'Pan', 'ClickedCallback', @PanWindow);
[mCData, mCM] = imread([sPath 'left.gif']); mCData = ind2rgb(mCData, mCM); mCData(mCData == 1) = NaN; % pan left
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_PanLeft', 'TooltipString', 'Pan left', 'ClickedCallback', @PanLeft);
[mCData, mCM] = imread([sPath 'right.gif']); mCData = ind2rgb(mCData, mCM); mCData(mCData == 1) = NaN; % pan right
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_PanRight', 'TooltipString', 'Pan right', 'ClickedCallback', @PanRight);
tZoomX = load([sPath 'zoomx.mat']); % zoom X
uipushtool('Parent', hToolbar, 'cdata', tZoomX.cdata, 'Tag', 'Spiky_WaitbarAction_ZoomX', 'TooltipString', 'Zoom X-scale', 'ClickedCallback', @ZoomRange);
tZoomY = load([sPath 'zoomy.mat']); % zoom Y
uipushtool('Parent', hToolbar, 'cdata', tZoomY.cdata, 'Tag', 'Spiky_WaitbarAction_ZoomX', 'TooltipString', 'Zoom Y-scale', 'ClickedCallback', @ZoomAmplitude);
mCData = im2double(imread([sPath 'tool_text.png'])); mCData(mCData == 0) = NaN; % set threshold
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_Threshold', 'TooltipString', 'Set threshold', 'ClickedCallback', @SetSpikeThreshold, 'separator', 'on');
[mCData, mCM] = imread([sPath 'tools_table.gif']); mCData = ind2rgb(mCData, mCM); mCData(mCData == 1) = NaN; % experiment variables
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_ExperimentVariables', 'TooltipString', 'Edit experiment variables', 'ClickedCallback', @ExperimentDescriptions);
mCData = im2double(imread([sPath 'tool_waveforms.png'])); mCData(mCData == 0) = NaN; % view waveforms
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_Waveforms', 'TooltipString', 'View sorted waveforms', 'ClickedCallback', @ShowOverlappingSpikeClusters, 'separator', 'on');
mCData = im2double(imread([sPath 'tool_psth.png'])); mCData(mCData == 0) = NaN; % plot psth
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_PSTH', 'TooltipString', 'Plot PSTH', 'ClickedCallback', @PlotPSTH);
mCData = im2double(imread([sPath 'tool_pc.png'])); mCData(mCData == 0) = NaN; % principal components
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_PC', 'TooltipString', 'Show principal components', 'ClickedCallback', @ViewWaveformPCs);

% Set Spiky! to by default NOT be in MERGE MODE
global g_bMergeMode g_bBatchMode
g_bMergeMode = 0;
g_bBatchMode = false;

% Figure menu
%  - File menu
hFile  = uimenu(g_hSpike_Viewer, 'Label', '&File');
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', '&Open...', 'Callback', @OpenFile, 'Accelerator', 'O');
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', 'Open &Directory...', 'Callback', @SetDirectory, 'Accelerator', 'D');
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', 'Open Directory Tree...', 'Callback', @SetDirectoryTree, 'Accelerator', 'T');
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', 'Open Settings...', 'Callback', @OpenSettings);
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', '&Save', 'Callback', @SaveResults, 'separator', 'on', 'Accelerator', 'S');
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', '&Publish to NEX', 'Callback', 'publish_spiky_data');
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', '&Import...', 'Callback', @ImportData, 'Accelerator', 'I');
hFileList = uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', 'Files', 'Tag', 'MenuFileList', 'separator', 'on'); % filelist
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', 'Previous...', 'Callback', 'spiky(''OpenFile([],-1)'');'); % previous
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', '&Next...', 'Callback', 'spiky(''OpenFile([],0)'');', 'Accelerator', 'N'); % next
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', 'Pa&ge Setup...', 'Callback', 'pagesetupdlg(gcbf)', 'separator', 'on');
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', 'Print Pre&view...', 'Callback', 'printpreview(gcbf)');
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', '&Print...', 'Callback', 'printdlg(gcbf)', 'Accelerator', 'P');
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', '&Invert Colors', 'Callback', @InvertColors, 'Checked', 'off');
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', 'Export &Figure...', 'Callback', @FigureExport);
uimenu(g_hSpike_Viewer, 'Parent', hFile, 'Label', 'E&xit Spiky!', 'Callback', @ExitSpiky, 'separator', 'on', 'Accelerator', 'Q');

%  - Edit menu
hEdit  = uimenu(g_hSpike_Viewer, 'Label', '&Edit');
uimenu(g_hSpike_Viewer, 'Parent', hEdit, 'Label', '&Undo', 'Callback', @Undo, 'Accelerator', 'Z');
uimenu(g_hSpike_Viewer, 'Parent', hEdit, 'Label', '&Redo', 'Callback', @Redo, 'Accelerator', 'Y');
uimenu(g_hSpike_Viewer, 'Parent', hEdit, 'Label', '&DAQ Info', 'Callback', @ShowDAQInfo, 'separator', 'on');

%  - View menu
hView  = uimenu(g_hSpike_Viewer, 'Label', '&View');
uimenu(g_hSpike_Viewer, 'Parent', hView, 'Label', '&Show Channels...', 'Callback', @SelectChannels);
uimenu(g_hSpike_Viewer, 'Parent', hView, 'Label', 'Show &Events', 'Callback', @ShowEvents, 'Checked', 'off', 'Accelerator', 'E');
uimenu(g_hSpike_Viewer, 'Parent', hView, 'Label', 'Zoom &In', 'Callback', @ZoomIn, 'separator', 'on');
uimenu(g_hSpike_Viewer, 'Parent', hView, 'Label', 'Zoom &Out', 'Callback', @ZoomOut);
uimenu(g_hSpike_Viewer, 'Parent', hView, 'Label', 'Zoom &Reset', 'Callback', @ZoomReset, 'Accelerator', 'X');
uimenu(g_hSpike_Viewer, 'Parent', hView, 'Label', 'Pan', 'Callback', @PanWindow, 'Accelerator', 'L');
uimenu(g_hSpike_Viewer, 'Parent', hView, 'Label', '&Zoom Range', 'Callback', @ZoomRange, 'Accelerator', 'Z');
uimenu(g_hSpike_Viewer, 'Parent', hView, 'Label', '&Zoom Amplitude', 'Callback', @ZoomAmplitude, 'Accelerator', 'Y');
uimenu(g_hSpike_Viewer, 'Parent', hView, 'Label', '&Normal Time', 'Callback', @NormalTime, 'separator', 'on', 'checked', 'on', 'Tag', 'ShowNormalTime');
uimenu(g_hSpike_Viewer, 'Parent', hView, 'Label', '&Grid', 'Tag', 'Spiky_Menu_ShowGrid', 'Callback', @ShowGrid);
uimenu(g_hSpike_Viewer, 'Parent', hView, 'Label', '&Refresh', 'Callback', @ViewTrialData, 'separator', 'on', 'Accelerator', 'R');

%  - Channels menu
hChannels = uimenu(g_hSpike_Viewer, 'Label', '&Channels');
uimenu(g_hSpike_Viewer, 'Parent', hChannels, 'Label', '&Channel Gains...', 'Callback', @SetGain, 'Accelerator', 'G');
uimenu(g_hSpike_Viewer, 'Parent', hChannels, 'Label', '&Channel Descriptions...', 'Callback', @ChannelDescriptions);
uimenu(g_hSpike_Viewer, 'Parent', hChannels, 'Label', 'Digiti&ze Channel', 'Callback', @DigitizeChannel, 'separator', 'on');
uimenu(g_hSpike_Viewer, 'Parent', hChannels, 'Label', '&Duplicate Channel', 'Callback', @DuplicateChannel);
uimenu(g_hSpike_Viewer, 'Parent', hChannels, 'Label', 'D&elete Channel', 'Callback', @DeleteChannel);
uimenu(g_hSpike_Viewer, 'Parent', hChannels, 'Label', '&PCA Noise Reduction... (B)', 'Callback', @PCACleaning);
uimenu(g_hSpike_Viewer, 'Parent', hChannels, 'Label', '&Invert Channel...', 'Callback', @InvertChannel);
uimenu(g_hSpike_Viewer, 'Parent', hChannels, 'Label', 'Experiment &Variables...', 'Callback', @ExperimentDescriptions, 'separator', 'on', 'Accelerator', 'V');
uimenu(g_hSpike_Viewer, 'Parent', hChannels, 'Label', '&Select Filtered Channels...', 'Callback', @SetFilterChannels, 'separator', 'on', 'Accelerator', 'F');
uimenu(g_hSpike_Viewer, 'Parent', hChannels, 'Label', '&Filter Options...', 'Callback', @FilterOptions );

%  - Spikes menu
hWaveforms  = uimenu(g_hSpike_Viewer, 'Label', '&Spikes');
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', 'Set &Threshold...', 'Callback', @SetSpikeThreshold );
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', '&Auto-Detect Thresholds', 'Callback', @AutoDetectThresholds );
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', '&Run Spike Detection (B)', 'Callback', @DetectSpikes, 'separator', 'on');
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', '&Dejitter Spikes (B)', 'Callback', @DejitterSpikes);
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', 'Remove &Outliers (B)', 'Callback', @RemoveOutlierSpikes);
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', 'Add &Multitrode...', 'Callback', @AddTetrode, 'separator', 'on', 'Accelerator', 'T');
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', 'Remove Multitrode...', 'Callback', @RemoveTetrode );
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', 'Set Max &Jitter...', 'Callback', @SetMaxSpikeJitter );
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', 'Show Waveforms', 'Callback', @ViewWaveforms, 'separator', 'on');
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', '&Channel Statistics', 'Callback', @ChannelStatistics);
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', 'Set Pre-Trigger...', 'Callback', @EditPreTriggerTime, 'separator', 'on');
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', 'Set Post-Trigger...', 'Callback', @EditPostTriggerTime);
uimenu(g_hSpike_Viewer, 'Parent', hWaveforms, 'Label', 'Set Spike Deadtime...', 'Callback', @SetDeadtime);

%  - Sorting menu
hSorting  = uimenu(g_hSpike_Viewer, 'Label', '&Sort');
uimenu(g_hSpike_Viewer, 'Parent', hSorting, 'Label', '&Sort Channel...', 'Callback', @SortSpikes);
uimenu(g_hSpike_Viewer, 'Parent', hSorting, 'Label', 'Sort &All Channels', 'Callback', @SortAllChannels);
uimenu(g_hSpike_Viewer, 'Parent', hSorting, 'Label', '&QuickSort!', 'Callback', @QuickSort);
uimenu(g_hSpike_Viewer, 'Parent', hSorting, 'Label', 'Show S&orted Spikes', 'Callback', @ShowSpikeClusters, 'separator', 'on');
uimenu(g_hSpike_Viewer, 'Parent', hSorting, 'Label', '&Cluster Control', 'Callback', @ShowOverlappingSpikeClusters, 'Accelerator', 'W');
uimenu(g_hSpike_Viewer, 'Parent', hSorting, 'Label', 'Show Aggregration &Tree', 'Callback', @ShowAggregationTree, 'Accelerator', 'A');
uimenu(g_hSpike_Viewer, 'Parent', hSorting, 'Label', 'View &Principal Components', 'Callback', @ViewWaveformPCs);
uimenu(g_hSpike_Viewer, 'Parent', hSorting, 'Label', 'Plot &Rasters', 'Callback', @PlotRasters);
uimenu(g_hSpike_Viewer, 'Parent', hSorting, 'Label', 'Plot &Instantaneous Rate', 'Callback', @PlotInstSpikeRate);
uimenu(g_hSpike_Viewer, 'Parent', hSorting, 'Label', 'Delete sorting data...', 'Callback', @DeleteSortingData, 'separator', 'on');

%  - Analysis menu
hAnalysis = uimenu(g_hSpike_Viewer, 'Label', '&Analysis');
hContinuous = uimenu(g_hSpike_Viewer, 'Parent', hAnalysis, 'Label', '&Continuous');
uimenu(g_hSpike_Viewer, 'Parent', hContinuous, 'Label', 'P&ower Spectral Densities', 'Callback', @ShowSpectralDensity);
uimenu(g_hSpike_Viewer, 'Parent', hContinuous, 'Label', '&Event Triggered Average', 'Callback', @ShowEventTriggeredAverage);
uimenu(g_hSpike_Viewer, 'Parent', hContinuous, 'Label', '&Cross Correlations', 'Callback', @PlotCrossCorrelationsContinuous);
uimenu(g_hSpike_Viewer, 'Parent', hContinuous, 'Label', '&Histogram', 'Callback', @PlotHistogramContinuous);

hDiscrete = uimenu(g_hSpike_Viewer, 'Parent', hAnalysis, 'Label', '&Discrete');
uimenu(hDiscrete, 'Parent', hDiscrete, 'Label', '&Cross Correlations', 'Callback', @PlotCrossCorrelationsDiscrete);
uimenu(hDiscrete, 'Parent', hDiscrete, 'Label', '&Peristimulus Time Histograms (PSTH)', 'Callback', @PlotPSTH);
uimenu(hDiscrete, 'Parent', hDiscrete, 'Label', '&Spiketime Distributions', 'Callback', @PlotSpikeTimeDistributions);

%  - Merge menu
hMerge  = uimenu(g_hSpike_Viewer, 'Label', '&Merge');
uimenu(g_hSpike_Viewer, 'Parent', hMerge, 'Label', '&Distribute Settings', 'Callback', @DistributeSettings);
uimenu(g_hSpike_Viewer, 'Parent', hMerge, 'Label', '&Merge Files', 'Callback', @CreateMergeFile, 'Accelerator', 'M');
uimenu(g_hSpike_Viewer, 'Parent', hMerge, 'Label', '&Load Merge File...', 'Callback', @LoadMergeFile);

%  - Tools menu
hTools  = uimenu(g_hSpike_Viewer, 'Label', '&Tools');
hScripts = uimenu(g_hSpike_Viewer, 'Parent', hTools, 'Label', 'Scripts');
uimenu(g_hSpike_Viewer, 'Parent', hScripts, 'Label', '&Run Script... (B)', 'Callback', @RunScript);
uimenu(g_hSpike_Viewer, 'Parent', hScripts, 'Label', 'Run &Batch Script...', 'Callback', @RunBatchScript);
uimenu(g_hSpike_Viewer, 'Parent', hScripts, 'Label', 'Get Script Help...', 'Callback', @GetScriptHelp);

% Get list of existing scripts in ./scripts directory
sPath = which('spiky');
sPath = CheckFilename([sPath(1:end-7) 'scripts\']);
tFiles = dir(sPath);
bFirst = 1;
for f = 1:length(tFiles)
    if ~isempty(strfind(tFiles(f).name, '.m')) && isempty(strfind(tFiles(f).name, '.m~'))
        sName = strrep(tFiles(f).name(1:end-2), '_', ' ');
        vIndx = strfind(sName, ' ');
        sName([1 vIndx+1]) = upper(sName([1 vIndx+1]));
        if bFirst
            uimenu(g_hSpike_Viewer, 'Label', sName, 'Parent', hScripts, 'Callback', [sprintf('spiky(''RunScript(''''%s'''')'')', tFiles(f).name)], 'Separator', 'on');
            bFirst = 0;
        else
            uimenu(g_hSpike_Viewer, 'Label', sName, 'Parent', hScripts, 'Callback', [sprintf('spiky(''RunScript(''''%s'''')'')', tFiles(f).name)]);
        end
    end
end

uimenu(g_hSpike_Viewer, 'Parent', hTools, 'Label', '&Batch Redo', 'Callback', 'spiky(''BatchRedo([],''''redo'''')'');', 'Accelerator', 'B');
uimenu(g_hSpike_Viewer, 'Parent', hTools, 'Label', '&Autoload New Files...', 'Callback', @AutoloadNewFiles, 'Separator', 'on');
uimenu(g_hSpike_Viewer, 'Parent', hTools, 'Label', '&Keyboard Mode...', 'Callback', @KeyboardMode, 'Accelerator', 'K');
uimenu(g_hSpike_Viewer, 'Parent', hTools, 'Label', '&Clear Persistent Variables', 'Callback', 'clear functions');

%  - Window menu
hWindow  = uimenu(g_hSpike_Viewer, 'Label', '&Window');
uimenu(g_hSpike_Viewer, 'Parent', hWindow, 'Label', 'Close All', 'Callback', @CloseSpikyWindows);
uimenu(g_hSpike_Viewer, 'Parent', hWindow, 'Label', 'Cascade', 'Callback', @CascadeSpikyWindows);
uimenu(g_hSpike_Viewer, 'Parent', hWindow, 'Label', 'Figures', 'Callback', @BringSpikyFiguresToFront);

%  - Help menu
hHelp  = uimenu(g_hSpike_Viewer, 'Label', '&Help');
uimenu(g_hSpike_Viewer, 'Parent', hHelp, 'Label', '&Spiky Website', 'Callback', 'web(''http://code.google.com/p/spiky/'', ''-browser'')');
uimenu(g_hSpike_Viewer, 'Parent', hHelp, 'Label', '&Online Documentation', 'Callback', 'web(''http://code.google.com/p/spiky/wiki/'', ''-browser'')', 'Accelerator', 'H');
uimenu(g_hSpike_Viewer, 'Parent', hHelp, 'Label', '&About Spiky!', 'Callback', @AboutSpiky, 'separator', 'on');
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PanRight(varargin)
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
vChild = findobj(get(findobj('Tag', 'Spiky'), 'children'), 'Type', 'axes'); % axes handles
vX = get(vChild(end), 'xlim');
FV.vXlim = vX+(diff(vX)/4); % pan right by 25%
SetStruct(FV); ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PanLeft(varargin)
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
vChild = findobj(get(findobj('Tag', 'Spiky') ,'children'), 'Type', 'axes'); % axes handles
vX = get(vChild(end), 'xlim');
FV.vXlim = vX-(diff(vX)/4); % pan left by 25%
SetStruct(FV); ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomReset(varargin)
[FV, hWin] = GetStruct;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomIn(varargin)
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
vChild = findobj(get(findobj('Tag', 'Spiky'), 'children'), 'Type', 'axes'); % axes handles
vX = get(vChild(end), 'xlim');
FV.vXlim = [mean(vX)-(diff(vX)/4) mean(vX)+(diff(vX)/4)];
SetStruct(FV); ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomOut(varargin)
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
vChild = findobj(get(findobj('Tag', 'Spiky') ,'children'), 'Type', 'axes'); % axes handles
vX = get(vChild(end), 'xlim');
FV.vXlim = [mean(vX)-(diff(vX)*2) mean(vX)+(diff(vX)*2)];
SetStruct(FV); ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomRange(varargin)
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
if FV.bPanOn, pan off, end
hold on; hLin = plot(NaN,NaN);

% If panning is currently enabled we need to turn that off in order
% to detect mouse presses in window
set(findobj(get(gcf,'children'), 'type', 'axes'), 'ButtonDownFcn', '')

% Point 1
waitforbuttonpress
global pnt1
pnt1 = get(gca,'CurrentPoint'); % button down detected

% Callback that draws the line between points
set(gcf, 'windowButtonMotionFcn', @ZoomRangeUpdateLine)

try
    waitforbuttonpress;
catch
    return;
end
pnt2 = get(gca,'CurrentPoint'); % button down detected

set(gcf, 'windowButtonMotionFcn', ''); % reset motion callback

hLine = findobj('Tag', 'SpikyXZoomIndicator');
if isempty(hLine), return, end

FV.vXlim = sort(get(hLine(3), 'xdata'));
delete(hLine)

SetStruct(FV); ViewTrialData


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomRangeUpdateLine(varargin)
global pnt1
pnt2 = get(gca, 'CurrentPoint'); % button down detected

hLine = findobj('Tag', 'SpikyXZoomIndicator');
hLine = flipud(hLine);
if isempty(hLine)
    hLine(1) = line(0,0);
    hLine(2) = line(0,0);
    hLine(3) = line(0,0);
    nYPos = get(gca, 'ylim');
    set(gca,'color', [.2 .1 .1])
    set(hLine(1), 'Tag', 'SpikyXZoomIndicator', 'color', [.6 .6 .6], 'ydata', [nYPos(1)/2 nYPos(1)/2])
    set(hLine(2), 'Tag', 'SpikyXZoomIndicator', 'color', [.6 .6 .6], 'ydata', [(nYPos(1)/2)-nYPos(1)/4 (nYPos(1)/2)+nYPos(1)/4], 'xdata', [pnt1(1,1) pnt1(1,1)])
    set(hLine(3), 'Tag', 'SpikyXZoomIndicator', 'color', [.6 .6 .6], 'ydata', [(nYPos(1)/2)-nYPos(1)/4 (nYPos(1)/2)+nYPos(1)/4], 'xdata', [pnt1(1,1) pnt1(1,1)])
end
set(hLine(1), 'xdata', [pnt1(1,1) pnt2(1,1)])
set(hLine(3), 'xdata', [pnt2(1,1) pnt2(1,1)])
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PCACleaning(varargin)
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
global g_bBatchMode
% Build a matrix containing the continuous data vectors of selected
% electrode. Note that the selected electrode must contain at least 3
% channels (i.e. be a tri-trode, tetrode or more. Stereotrodes and
% individual electrodes cannot be cleaned!)
if ~g_bBatchMode
    if isfield(FV, 'tExperimentVariables')
        nIndx = find(strcmp({FV.tExperimentVariables.sVariable}, 'bPCACleanedChannels'));
        if ~isempty(nIndx)
            if FV.tExperimentVariables(nIndx).sValue
                switch questdlg('This file has already been PCA cleaned. If it is necessary to redo the cleaning it is recommended you first re-load the original files from disk.', ...
                        'Spiky!', 'Continue', 'Abort', 'Abort')
                    case 'Abort', return
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
    sAns = questdlg('Do you want to plot and compare PCA cleaned traces?', 'Spiky!', 'No', 'Yes', 'Never plot', 'Yes');
else
    sAns = 'No';
end
if strcmpi(sAns, 'Never plot'), p_bNeverPlotPCA = 1; end

if strcmpi(sAns, 'Yes') && ~p_bNeverPlotPCA && ~g_bBatchMode
    hFig = figure;
    set(hFig, 'color', [.2 .2 .2], 'name', 'Spiky! PCA Cleaning', 'NumberTitle', 'off')
    
    for c = 1:size(mRaw,2)
        hAx = axes('position', [.05 (1/size(mRaw,2))*(c-1) .95 1/size(mRaw,2)]);
        set(hAx, 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'xtick', [])
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
    sAns = questdlg('Keep cleaned data? Keeping will overwrite the data vectors in this .spb file. Original, raw data is not over-written.', 'Spiky!', 'Yes', 'No', 'Cancel', 'Yes');
end
switch sAns
    case 'Yes'
        % Replace raw data with cleaned data in FV
        for i = 1:length(csChannels)
            FV.tData.(csChannels{i}) = mClean(:,i)';
        end

        % Log this change in FV.tExperimentVariables
        if isfield(FV, 'tExperimentVariables')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowEvents(varargin)
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
% update menu item
global g_hSpike_Viewer % handle to Spiky! window
hMenuItem = findobj(g_hSpike_Viewer, 'Label', 'Show &Events'); % handle of 'Show Event'
switch get(hMenuItem, 'Checked')
    case 'on' % turn off if on
        set(hMenuItem, 'Checked', 'off');
    case 'off' % turn on if off
        set(hMenuItem, 'Checked', 'on');
end
ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveResults(varargin)
% Save trial settings
global g_bMergeMode g_hSpike_Viewer
if ~CheckDataLoaded, return, end
sPointer = get(g_hSpike_Viewer,'Pointer');
set(g_hSpike_Viewer,'Pointer','watch')
[FV, hWin] = GetStruct;
sPath = [FV.sLoadedTrial(1:end-4) '.spb'];
if ~g_bMergeMode % remove raw data so its not duplicated on disk
    %FV = rmfield(FV, 'tData');
end
%save(sPath, 'FV', '-MAT', '-V6') % no compression
save(sPath, 'FV', '-MAT')
% Update the title of the Spiky! window to reflect that settings were saved
global g_hSpike_Viewer
if length(FV.sLoadedTrial) > 70
    sLoadedTrial = ['...' FV.sLoadedTrial(end-69:end)];
else sLoadedTrial = FV.sLoadedTrial; end
set(g_hSpike_Viewer, 'Name', CheckFilename(sprintf('Spiky! - %s', sLoadedTrial)))
set(g_hSpike_Viewer,'Pointer',sPointer)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OpenSettings(varargin)
% OpenSettings loads the .spb file of a .mat file. This file contains various
% settings and results
persistent p_bAlwaysUseDuplicate
if isempty(p_bAlwaysUseDuplicate) p_bAlwaysUseDuplicate = 0; end
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
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
load(sPath, 'FV', '-MAT')

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
                        sAns = questdlg('A duplicate of the raw data exists in an existing Spiky! generated file. Do you wish to use this data, or reload all raw data from the .daq file? Note that reloading raw data will erase imported fields.', ...
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
                        cFieldnamesDAQ{f}
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AutoloadNewFiles(varargin)
[FV, hWin] = GetStruct;
if ~isfield(FV, 'sDirectory')
    warndlg('You must first select the directory that should be monitored in File->Open Directory', 'Spiky!')
else
    map_autoload(CheckFilename([FV.sDirectory '\']))
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetFilterChannels(varargin)
% Select which electrodes carry EMG/LFP data
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
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

for i = 1:length(csDisplayChannels)
    sBoxStr = sprintf('%s (%s)', csChannelDescriptions{i}, csDisplayChannels{i});
    if any(strcmp(FV.csEMGChannels, csDisplayChannels{i}))
        uicontrol(hFig, 'Style', 'checkbox', 'Position', [10 nCL 150 20], 'String', sBoxStr, 'HorizontalAlignment', 'left', 'backgroundcolor', [.8 .8 .8], 'value', 1);
    else
        uicontrol(hFig, 'Style', 'checkbox', 'Position', [10 nCL 150 20], 'String', sBoxStr, 'HorizontalAlignment', 'left', 'backgroundcolor', [.8 .8 .8], 'value', 0);
    end
    nCL = nCL - 25;
end
uicontrol(hFig, 'Style', 'pushbutton', 'Position', [50 nCL 50 20], ...
    'Callback', 'global g_vSelCh; g_vSelCh=flipud(get(get(gcf,''children''),''value'')); g_vSelCh=[g_vSelCh{:}];close(gcf)', ...
    'String', 'OK' ); % OK button
uiwait % wait for user to close window or click OK button
g_vSelCh = g_vSelCh(1:end-1);
FV.csEMGChannels = csDisplayChannels(find(g_vSelCh)); % selected channels
FV.csDisplayChannels = unique([FV.csDisplayChannels FV.csEMGChannels]); % selected channels
clear global g_vSelCh
SetStruct(FV)
ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FilterOptions(varargin)
[FV, hWin] = GetStruct;
cPrompt = {'High-pass (Hz)','Low-pass (Hz)', 'Rectify (1=yes, 0=no)'};
cAnswer = inputdlg(cPrompt, 'Filter Options',1, {num2str(FV.nEMGHiPass), num2str(FV.nEMGLoPass), num2str(FV.nEMGRectify)});
if isempty(cAnswer), return, end
FV.nEMGHiPass = str2num(cAnswer{1});
FV.nEMGLoPass = str2num(cAnswer{2});
FV.nEMGRectify = str2num(cAnswer{3});
SetStruct(FV)
ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetDirectory(varargin)
% Set the current directory
FV = SetFVDefaults(); % initialize FV with default settings
global g_bBatchMode
sDirectory = uigetdir;           % path to session directory (without last slash at end)
if sDirectory == 0, return, end
cd(sDirectory)
% Get list of DAQ files
sFiles = dir(CheckFilename(sprintf('%s\\*.daq', sDirectory)));
if isempty(sFiles)
    warndlg('No DAQ files were found in the chosen directory.', 'No files found')
    return
end
global g_hSpike_Viewer
hFileList = findobj(g_hSpike_Viewer, 'Tag', 'MenuFileList');
% Clear current file list
hFileListChildren = findobj(g_hSpike_Viewer, 'Parent', hFileList);
delete(hFileListChildren)
for f = 1:size(sFiles, 1)
    uimenu(g_hSpike_Viewer, 'Parent', hFileList, 'Label', sFiles(f).name, ...
        'UserData', CheckFilename([sDirectory '\\' sFiles(f).name]), ...
        'Callback', sprintf('spiky(''OpenFile([], %d)'');', f) );
end
FV.sDirectory = sDirectory;
SetStruct(FV)
OpenFile([],1); % Load first trial
g_bBatchMode = false;
ViewTrialData   % Update GUI
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetDirectoryTree(varargin)
% Load list of all DAQ files from all descending directories below a top directory
% Set the top directory of the tree
FV = SetFVDefaults();
sDirectory = uigetdir;
global g_hSpike_Viewer
if sDirectory == 0, return, end

% Iterate recursively over all sub-directories and extract all .DAQ files
[cPaths, cFiles] = GetFilePaths(sDirectory, '.daq');

% Display error if no files were found
if isempty(cFiles)
    warndlg('No DAQ files were found in the chosen directory.', 'Spiky!')
    return
end

% Reset current file list
hFileList = findobj(g_hSpike_Viewer, 'Tag', 'MenuFileList');
hFileListChildren = findobj(g_hSpike_Viewer, 'Parent', hFileList);
delete(hFileListChildren)

% Update GUI
for f = 1:length(cFiles)
    uimenu(g_hSpike_Viewer, 'Parent', hFileList, 'Label', cFiles{f}, ...
        'UserData', CheckFilename([cPaths{f} '\\' cFiles{f}]), ...
        'Callback', sprintf('spiky(''OpenFile([], %d)'');', f) );
end
FV.sDirectory = sDirectory;
SetStruct(FV)
OpenFile([],1); % Load first trial
ViewTrialData   % Update GUI

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cPaths, cFiles] = GetFilePaths(sBaseDir, sSuffix)
% --- Get paths of videos in a directory and its sub-directories
% Recursive search for .bin files from a selected directory
% inputs: sSuffix    Find files this file extension (e.g. '.avi')
%         sBaseDir   Base directory
cPaths = {};
cFiles = {};
tDirList = dir(sBaseDir);
for t1 = 3:length(tDirList)
    if tDirList(t1).isdir % depth = 1
        % recursively call this function
        sBaseDirRecurs = CheckFilename([sBaseDir '/' tDirList(t1).name]);
        [cPaths2, cFiles2] = GetFilePaths(sBaseDirRecurs, sSuffix);
        cPaths = {cPaths{:}, cPaths2{:}};
        cFiles = {cFiles{:}, cFiles2{:}};
    else
        % Add file names to paths cell
        if strcmp(tDirList(t1).name(end-3:end), sSuffix)
            %sFilename = CheckFilename([sBaseDir '/' tDirList(t1).name]);
            cPaths{end+1} = sBaseDir;
            cFiles{length(cPaths)} = tDirList(t1).name;
        end
    end
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ViewTrialData(varargin)
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
% Dont display any trial data if Spiky! is in MERGE MODE
global g_bMergeMode g_bBatchMode
if ~g_bMergeMode && ~g_bBatchMode % destroy waitbar in case we are not in Merge or Batch mode
    SpikyWaitbar(1,1);
end 

hFig = findobj('Tag', 'Spiky');

% Delete all current axes objects
delete(findobj(hFig, 'type', 'axes'))

% Delete all current uicontextmenu objects
delete(findobj(hFig, 'type', 'uicontextmenu'))

% Reset click callback (unless we're in pan mode)
hPan = pan(hFig);
if strcmp(get(hPan, 'Enable'), 'off')
    set(hFig, 'WindowButtonDownFcn', '');
end
zoom off

% Select channels
if isempty(g_bMergeMode), g_bMergeMode = 0; end
if isempty(FV.csDisplayChannels) && ~g_bMergeMode, SelectChannels; end
[FV, hWin] = GetStruct;

% Check that all channel names exist in FV.csChannels
cFields = fieldnames(FV.tData);
for i = 1:length(cFields)
    if ~isempty(strfind(cFields{i}, '_TimeBegin'))
        FV.csChannels{end+1} = cFields{i}(1:end-10);
        if isfield(FV.tData, cFields{i}(1:end-10))
            sField = cFields{i}(1:end-10);
            vMinMax = [min(FV.tData.(sField)) max(FV.tData.(sField))];
            [vMinMax, FV] = AdjustChannelGain(FV, vMinMax, sField); % mV
            FV.tYlim(1).(FV.csChannels{end}) = vMinMax;
        else
            FV.tYlim(1).(FV.csChannels{end}) = []; % assume its digital
        end
    end
end
FV.csChannels = unique(FV.csChannels);

% Check that all DisplayChannels exist (and remove those that dont)
vRemIndx = [];
for i = 1:length(FV.csDisplayChannels)
    if ~isfield(FV.tData, FV.csDisplayChannels{i}), vRemIndx = [vRemIndx i]; end
end
FV.csDisplayChannels(vRemIndx) = [];
SetStruct(FV)

hSubplots = [];
bShowDigitalEvents = strcmp(get(findobj(hWin, 'Label', 'Show &Events'),'checked'), 'on');
if isempty(FV.csDigitalChannels), bShowDigitalEvents = 0; end
for i = 1:length(FV.csDisplayChannels)
    sCh = FV.csDisplayChannels{i};
    vCont = FV.tData.(sCh); % continuous trace (V)
    if length(vCont) <= 1 && ~FV.bPlotRasters
        uiwait(warndlg(sprintf('Cannot display channel %s as it is not a vector.', sCh)))
        continue
    end

    % Get channel description
    sYLabel = GetChannelDescription(sCh);

    if ( ~FV.bPlotRasters || ~isfield(FV.tSpikes, sCh) ) || ~g_bMergeMode
        [vCont, FV] = AdjustChannelGain(FV, vCont, sCh); % mV
        [FV,hWin] = GetStruct; % reload FV since channel gains may have changed in previous line
        
        nFs = FV.tData.([sCh '_KHz']) * 1000; % sampling frequency (Hz)
        %nFs = 1.0178*2500;
        nBeginTime = FV.tData.([sCh '_TimeBegin']); % start of sampling (sec)
        nEndTime = FV.tData.([sCh '_TimeEnd']); % start of sampling (sec)
        vTime = linspace(nBeginTime, nEndTime, length(vCont));

        % Alternative way to compute vTime  - Corrected 022312
        %vTime = nBeginTime:(1/nFs):( nBeginTime + (1/nFs)*(length(vCont)-1) );
        
        % Limit vector to FV.vXlim
        if ~isempty(FV.vXlim)
            vXlim(1) = max([1 round((FV.vXlim(1)-nBeginTime)/diff(vTime(1:2)))]);
            vXlim(2) = min([length(vTime) round((FV.vXlim(2)-nBeginTime)/diff(vTime(1:2)))]);
            if ~(vXlim(2) <= vXlim(1))
                vTime = vTime(vXlim(1):vXlim(2));
                vCont = vCont(vXlim(1):vXlim(2));
            end
        end
        % Decimate continuous trace
        nLen = 5000;
        vShowIndx = unique(round(linspace(1,length(vTime),nLen)));
    end

    if ~FV.bPlotRasters || ~isfield(FV.tSpikes, sCh)
        %if (nLen) < length(vCont)
            vTimeAllTrace = vTime;
            vTime = vTime(vShowIndx);
            vContAllTrace = vCont;
            vCont = vCont(vShowIndx);
        %end

        % filter and rectify EMG signals
        if any(strcmp(FV.csEMGChannels, sCh))
            if nFs < 2500
                [vCont, vTime, nNewFs] = FilterChannel(vCont, vTime, nFs, FV.nEMGLoPass, FV.nEMGHiPass, FV.nEMGRectify, 'none');
            else
                [vCont, vTime, nNewFs] = FilterChannel(vCont, vTime, nFs, FV.nEMGLoPass, FV.nEMGHiPass, FV.nEMGRectify, 'decimate');
            end
            FV.tYlim.(sCh) = [min(vCont) max(vCont)];
        end
    end

    % continuous trace
    nSubs = length(FV.csDisplayChannels);

    if bShowDigitalEvents, nSubs = nSubs + 1; end
    nSubHeight = (.91-(nSubs*.02)) / nSubs;
    nSubY = .075 + [.02*(i-1)] + [nSubHeight*(i-1)];
    hSubplots(end+1) = axes('position', [.06 nSubY .91 nSubHeight+.01]);
    
    % Plot spike rasters
    if FV.bPlotRasters && isfield(FV.tSpikes, sCh)
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
                %nRAdj = abs(diff(FV.vXlim)*.5); % adjust extended range by 50% for panning purposes
                %vIndx = vSpiketimes < [min(FV.vXlim)-nRAdj] | vSpiketimes > [max(FV.vXlim)+nRAdj];
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
                if diff(FV.vXlim) < 1 % 1 ms spikes when timescale is < 10 sec
                    vIndices(2:2:end) = vIndices(1:2:end) + 0.001;
                elseif diff(FV.vXlim) < 5 % 2 ms spikes when timescale is < 10 sec
                    vIndices(2:2:end) = vIndices(1:2:end) + 0.003;
                elseif diff(FV.vXlim) < 10 % 5 ms spikes when timescale is < 20 sec
                    vIndices(2:2:end) = vIndices(1:2:end) + 0.01;
                elseif diff(FV.vXlim) < 20 % 5 ms spikes when timescale is < 20 sec
                    vIndices(2:2:end) = vIndices(1:2:end) + 0.02;
                elseif diff(FV.vXlim) < 50 % 5 ms spikes when timescale is < 20 sec
                    vIndices(2:2:end) = vIndices(1:2:end) + 0.04;
                else
                    vIndices(2:2:end) = vIndices(1:2:end) + 0.2;
                end
                
                % Rearrange indices
                vYs = repmat([nRow nRow NaN NaN], 1, length(vIndices)/4);
                % Limit each spike to 1 msec width
                nSpikeHeight = 50; % adjust this number to change height of spikes
                hLin = plot(vIndices, vYs, 'linewidth', nSpikeHeight/length(vUnits), 'color', mCol);
                mUserData = [vIndices vYs'];
            end
            
            hold on
            hMenu = uicontextmenu;
            set(hMenu, 'userdata', mUserData, 'Tag', [sCh '-' num2str(vUnits(nU))]);
            set(hMenu, 'Tag', [sCh '-' num2str(vUnits(nU))]);
            hItems(1) = uimenu(hMenu, 'Label', '&Copy to Figure', 'Callback', 'figure;mD=get(get(gcbo, ''parent''),''userdata'');plot(mD(:,1),mD(:,2),''linewidth'',1)');
            hItems(end+1) = uimenu(hMenu, 'Label', '&Hide', 'Callback', @SelectChannels, 'Tag', sCh);
            hItems(end+1) = uimenu(hMenu, 'Label', '&Delete Section', 'Callback', @DeleteSection, 'Tag', sCh, 'Separator', 'on');
            hItems(end+1) = uimenu(hMenu, 'Label', '&Shift Time', 'Callback', @ShiftBeginTime, 'Tag', sCh);
            hItems(end+1) = uimenu(hMenu, 'Label', 'Show &Markers', 'Callback', 'set( get(gcbo, ''userdata''), ''marker'', ''o'')', 'userdata', hLin, 'Separator', 'on');
            hItems(end+1) = uimenu(hMenu, 'Label', 'Meas&ure', 'Callback', @MeasureLine, 'Tag', sCh);
            hItems(end+1) = uimenu(hMenu, 'Label', 'Digiti&ze Channel (Manual)', 'Callback', @DigitizeChannel, 'Tag', sCh);
            hItems(end+1) = uimenu(hMenu, 'Label', '&Properties', 'Separator', 'on', 'Callback', @ChannelProperties, 'Tag', sCh);
            hItems(end+1) = uimenu(hMenu, 'Label', '&Set Electrode Position', 'Callback', @SetElectrodePosition, 'Tag', sCh);
            hItems(end+1) = uimenu(hMenu, 'Label', '&Set Receptive Field', 'Callback', @SetReceptiveField, 'Tag', sCh);
            set(hLin, 'uicontextmenu', hMenu)
            
            csUnits{end+1} = num2str(vUnits(nU));
        end
    else % plot continuous trace
        hLin = plot(vTime, vCont, 'color', FV.mColors(i,:));
        % attach context menu to line to enable additional options
        hMenu = uicontextmenu;
        if exist('vTimeAllTrace')
            set(hMenu, 'userdata', [vTimeAllTrace(:) vContAllTrace(:)], 'Tag', sCh);
        end
        set(hMenu, 'Tag', sCh);
        hItems(1) = uimenu(hMenu, 'Label', '&Copy to Figure', 'Callback', 'figure;mD=get(get(gcbo, ''parent''),''userdata'');plot(mD(:,1),mD(:,2));xlabel(''Time (s)'')');
        hItems(end+1) = uimenu(hMenu, 'Label', '&Hide', 'Callback', @SelectChannels, 'Tag', sCh);
        hItems(end+1) = uimenu(hMenu, 'Label', '&Delete Section', 'Callback', @DeleteSection, 'Tag', sCh, 'Separator', 'on');
        hItems(end+1) = uimenu(hMenu, 'Label', '&Shift Time', 'Callback', @ShiftBeginTime, 'Tag', sCh);
        hItems(end+1) = uimenu(hMenu, 'Label', 'Show &Markers', 'Callback', 'set( get(gcbo, ''userdata''), ''marker'', ''o'')', 'userdata', hLin, 'Separator', 'on');
        hItems(end+1) = uimenu(hMenu, 'Label', 'Meas&ure', 'Callback', @MeasureLine, 'Tag', sCh);
        hItems(end+1) = uimenu(hMenu, 'Label', 'Digiti&ze Channel', 'Callback', @DigitizeChannelAuto, 'Tag', sCh, 'Separator', 'on');
        hItems(end+1) = uimenu(hMenu, 'Label', 'Digiti&ze Manually', 'Callback', @DigitizeChannel, 'Tag', sCh);
        hItems(end+1) = uimenu(hMenu, 'Label', '&Properties', 'Separator', 'on', 'Callback', @ChannelProperties, 'Tag', sCh);
        set(hLin, 'uicontextmenu', hMenu)
    end
    
    % Set axis labels
    if i == 1, xlabel('Time (sec)'); end
    if FV.bPlotRasters && isfield(FV.tSpikes, sCh)
        hLabel = ylabel(sYLabel, 'FontSize', 8);
    else
        hLabel = ylabel([sYLabel ' (mV)'], 'FontSize', 8);
    end
    set(hLabel, 'Interpreter', 'none')
        
    set(hSubplots(end), 'Tag', sCh, 'color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'FontSize', 7)
    if g_bMergeMode
        if isempty(FV.vXlim), axis tight
        else set(hSubplots(end), 'xlim', FV.vXlim); end
    else
        nX_s = vTime(1); nX_e = vTime(end);
        set(hSubplots(end), 'xlim', [nX_s nX_e]);
    end

    % Set Y axis limits
    if ~isfield(FV, 'bPlotInstSpikeRate'), FV.bPlotInstSpikeRate = 0; end
    if FV.bPlotInstSpikeRate && isfield(FV.tSpikes, sCh)
        % do nothing for now
    elseif FV.bPlotRasters && isfield(FV.tSpikes, sCh)
        set(hSubplots(end), 'ylim', [.5 nRow+.5], 'ytick', 1:1:nRow, 'yticklabel', csUnits);
    elseif isfield(FV.tYlim, sCh)
        if ~isempty(FV.tYlim.(sCh)) && ~any(isnan(FV.tYlim.(sCh)))
            nYd = (diff(FV.tYlim.(sCh)) * .1) .* [-1 1];
            vYLim = FV.tYlim.(sCh) + nYd;
            if vYLim(1) == vYLim(2)
                nD = vYLim(1) * .05 - (10^-6);
                set(hSubplots(end), 'ylim',  vYLim - [-nD nD]); % use ylim +/- 5%
            else
                set(hSubplots(end), 'ylim',  vYLim); % use saved ylim * 5%
            end
        end
    end

    % Plot threshold line(s)
    if isfield(FV.tSpikeThresholdsPos, sCh) && ~FV.bPlotRasters % positive
        nThresh = FV.tSpikeThresholdsPos.(sCh);
        hold on;
        hThresh = plot([vTime(1) vTime(end)], [nThresh nThresh], ':', 'color', [.7 .7 .7]); hold off
        hCntxtMenu  = uicontextmenu;
        set(hThresh, 'UIContextMenu', hCntxtMenu, 'Tag', [sCh '_pos']);
        uimenu(hCntxtMenu, 'Label', 'Remove Threshold Line', 'Callback', @RemoveSpikeThreshold)
        uimenu(hCntxtMenu, 'Label', 'Set Threshold Manually', 'Callback', @SetThresholdManually)
    end
    if isfield(FV.tSpikeThresholdsNeg, sCh) && ~FV.bPlotRasters % negative
        nThresh = FV.tSpikeThresholdsNeg.(sCh);
        hold on;
        hThresh = plot([vTime(1) vTime(end)], [nThresh nThresh], ':', 'color', [.7 .7 .7]); hold off
        hCntxtMenu  = uicontextmenu;
        set(hThresh, 'UIContextMenu', hCntxtMenu, 'Tag', [sCh '_neg']);
        uimenu(hCntxtMenu, 'Label', 'Remove Threshold Line', 'Callback', @RemoveSpikeThreshold)
        uimenu(hCntxtMenu, 'Label', 'Set Threshold Manually', 'Callback', @SetThresholdManually)
    end

    % Plot spikes on top of continuous trace (only if rasters are not displayed)
    if isfield(FV.tSpikes, sCh) && ~FV.bPlotRasters
        %axes(hSubplots(end))
        vSpiketimes = FV.tSpikes.(sCh).spiketimes; % nFs resolution
        vSpiketimes = vSpiketimes ./ nFs;
        vSpiketimes(vSpiketimes < vTime(1) | vSpiketimes > vTime(end)) = [];
        vVolts = repmat(nThresh, length(vSpiketimes), 1);
        hold on; hDot = plot(hSubplots(end), vSpiketimes', vVolts', 'w.');        
        set(hDot, 'ButtonDownFcn', @IdentifySpike)
    end
    % Channel legend
    %if ishandle(hSubplots(end)), axes(hSubplots(end)), end
end

% Plot digital events
if bShowDigitalEvents
    % Create subplot
    if ~exist('nSubHeight', 'var')
        i = 0;
        nSubs = 1;
        nSubHeight = (.91-(nSubs*.02)) / nSubs;
    end
    nSubY = .1 + [.015*i] + [nSubHeight*i];
    hSubplots(end+1) = axes('position', [.06 nSubY .91 nSubHeight]);
    hold on
    
    % Iterate and plot digital events
    csEvents = {};
    nY = 0;
    for nF = 1:length(FV.csDigitalChannels) % iterate over fieldnames
        % Ignore DAQ_Start, DAQ_Stop and DAQ_Trigger (default fields generated by DA Toolbox)
        if any(strcmp(FV.csDigitalChannels(nF), {'DAQ_Start', 'DAQ_Stop', 'DAQ_Trigger'}))
            continue
        end
        
        csEvents{end+1} = [FV.csDigitalChannels{nF} '_Up']; % event name
        nFs = FV.tData.([FV.csDigitalChannels{nF} '_KHz']);
        vUpTimes = FV.tData.(csEvents{end}); % UP times, sec
        vUpTimes = sort(vUpTimes); % Sort
        nY = nY + 1;
        
        % Continue if there are no UP times
        if isempty(vUpTimes), continue, end
        
        vDownTimes = FV.tData.([FV.csDigitalChannels{nF} '_Down']); % DOWN times, sec
        vDownTimes = sort(vDownTimes); % Sort
        if isempty(vDownTimes), vDownTimes = FV.tData.([FV.csDigitalChannels{nF} '_TimeEnd']); end

        % Remove DOWN times that occur before UP events
        vDownTimes(vDownTimes <= vUpTimes(1)) = [];

        % If up and down times differ in length, crop the longest
        if length(vUpTimes) > length(vDownTimes)
            vUpTimes = vUpTimes(1:length(vDownTimes));
        elseif length(vUpTimes) < length(vDownTimes)
            vDownTimes = vDownTimes(1:length(vUpTimes));
        end

        % Continue to next channel if there are no events
        if isempty(vDownTimes) || isempty(vUpTimes), continue, end

        % Remove events outside of axis limits
        if ~isempty(FV.vXlim)
            %nRAdj = abs(diff(FV.vXlim)*.5); % adjust extended range by 50% for panning purposes
            %vIndx = find((vUpTimes < [min(FV.vXlim)-nRAdj] | vUpTimes > [max(FV.vXlim)+nRAdj]) ...
            %    & (vDownTimes < [min(FV.vXlim)-nRAdj] | vDownTimes > [max(FV.vXlim)+nRAdj]));
            vIndx = find((vUpTimes < [min(FV.vXlim)] | vUpTimes > [max(FV.vXlim)]) ...
                & (vDownTimes < [min(FV.vXlim)] | vDownTimes > [max(FV.vXlim)]));
            vUpTimes(vIndx) = [];
            vDownTimes(vIndx) = [];
        end

        % Continue to next channel if there are no events
        if isempty(vDownTimes) || isempty(vUpTimes), continue, end

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
        if length(vIndices) > 1000
            vPlotIndx = 1:round(length(vIndices)/1000):length(vIndices);
        else
            vPlotIndx = 1:length(vIndices);
        end
        plot(vIndices(vPlotIndx), vYs(vPlotIndx), 'linewidth', 5, 'color', [1 1 0])
        
    end
    
    % Plot file start and end markers
    if isfield(FV.tData, 'FileStart') % only in Merge Mode
        for i = 1:length(FV.tData.FileStart)
            nX = FV.tData.FileStart(i).Timestamp;
            if ~isempty(FV.vXlim)
                nRAdj = abs(diff(FV.vXlim)*.5); % adjust extended range by 50% for panning purposes
                if nX < (min(FV.vXlim)-nRAdj) || nX > (max(FV.vXlim)+nRAdj), continue, end
            end
            vY = [0 nF+.5];
            hLin = plot([nX nX],  vY, 'color', 'g');
            set(hLin, 'DisplayName', FV.tData.FileStart(i).File, 'ButtonDownFcn', 'hLeg=legend(gcbo);set(hLeg,''interpreter'',''none'',''TextColor'',''w'',''FontSize'',10); legend boxoff;')
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
            set(hLin, 'DisplayName', FV.tData.FileEnd(i).File, 'ButtonDownFcn', 'hLeg=legend(gcbo);set(hLeg,''interpreter'',''none'',''TextColor'',''w'',''FontSize'',10); legend boxoff;')
        end
    end

    % Set axis properties
    set(hSubplots(end), 'Tag', 'Events', 'color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'FontSize', 7)

    % Y-tick labels
    cYTickLabels = {};
    for c = 1:length(FV.csDigitalChannels)
        % Ignore DAQ_Start, DAQ_Stop and DAQ_Trigger (default fields generated by DA Toolbox)
        if any(strcmp(FV.csDigitalChannels(c), {'DAQ_Start', 'DAQ_Stop', 'DAQ_Trigger'}))
            continue
        end
        cYTickLabels{end+1} = FV.csDigitalChannels{c};
        if isfield(FV, 'tChannelDescriptions')
            if isfield(FV.tChannelDescriptions, 'sChannel')
                nIndx = find(strcmp({FV.tChannelDescriptions.sChannel}, cYTickLabels{end}));
                if isempty(nIndx), continue, end
                if ~isempty(FV.tChannelDescriptions(nIndx(1)).sDescription)
                    cYTickLabels{end} = FV.tChannelDescriptions(nIndx(1)).sDescription;
                end
            end
        end
    end
    
    axis(hSubplots(end), 'tight')    
    if isempty(cYTickLabels)
        set(hSubplots(end), 'ylim', [0 1], 'ytick', [], 'yticklabel', [], 'color', [.2 .1 .1])
        text(mean(get(hSubplots(1), 'xlim')), .5, 'No digital events found', 'color', 'w', 'horizontalalignment','center')
    else
        if all(ishandle(hSubplots))
            set(hSubplots(end), 'ylim', [0.65 length(cYTickLabels)+.35], ...
                'ytick', 1:length(cYTickLabels), 'yticklabel', cYTickLabels, 'color', [.2 .1 .1])
        end
    end
    box on
    
    % Context menu in Events axis
    hMenu = uicontextmenu;
    uimenu(hMenu, 'Label', 'Event &Statistics', 'Callback', @ShowEventStatistics);
    uimenu(hMenu, 'Label', 'Undigitize Channel', 'Callback', @UndigitizeChannel);
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
hHandles = get(hFig,'children');
hAxes = findobj(hHandles,'type','axes');
hOthers = setdiff(hHandles, hAxes);
try % try here because sometimes returned handles arent axes after all...
    mPos = cell2mat(get(hAxes, 'position')); % axes positions
    [vY, vOrder] = sort(mPos(:,2)); % order by Y position
    set(hFig, 'children', [flipud(hOthers); hAxes(vOrder)])
end

SetStruct(FV)
PanMode()
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vSpiketimes = DropDeadtimeSpikes(vSpiketimes)
[FV, hWin] = GetStruct;
% round off spiketimes to deadtime resolution and keep first of all unique values
vSpiketimesRound = round(vSpiketimes ./ (FV.nDeadTimeMs/1000)) .* (FV.nDeadTimeMs/1000); % sec
[vSpiketimesRound, vRoundIndx] = unique(vSpiketimesRound, 'first');
vSpiketimes = vSpiketimes(vRoundIndx);
vISI = diff(vSpiketimes); % sec
vDropIndx = find(vISI <= (FV.nDeadTimeMs / 1000)); % ms
vSpiketimes(vDropIndx+1) = [];
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PanWindow(varargin)
[FV, hWin] = GetStruct;
global g_hSpike_Viewer % handle to Spiky! window
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PanMode(varargin)
[FV, hWin] = GetStruct;
if FV.bPanOn
    hPan = pan;
    set(hPan, 'enable', 'on', 'Motion', 'horizontal', 'ActionPostCallback', @UpdatePannedWindow)
    linkaxes(findobj(get(gcf,'children'), 'type', 'axes'), 'x')
else
    hPan = pan;
    set(hPan, 'enable', 'off')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdatePannedWindow(varargin)
% Wait until x-axis limits have changed
vOrigLims = get(gca, 'xlim');

% Get new x limits
[FV, hWin] = GetStruct;
vCurrLim = get(gca, 'xlim');
FV.vXlim = vCurrLim; % pan left by 25%
SetStruct(FV); ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vCont, vTime, nNewFs] = FilterChannel(vCont, vTime, nFs, nLoPass, nHiPass, bRectify, sOption)
vTimeOrig = vTime;

% Interpolate NaN indices
vNaNIndx = isnan(vCont);
vCont(vNaNIndx) = interp1(find(~vNaNIndx), vCont(~vNaNIndx), find(vNaNIndx), 'linear');

switch lower(sOption)
    case 'decimate'
        % decimate (i.e. resample signal)
        nNewFs = nLoPass * 2.5;
        nStep = ceil(nFs / nNewFs);
        nNewFs = nFs/nStep; % new Fs must be an integer
        if ~isempty(vTime), vTime = vTime(1:nStep:end); end
        vCont = vCont(1:nStep:end);
    otherwise
        nNewFs = nFs;
        vTime = [];
end

% High pass filter
if nHiPass > 0
    [b,a] = butter(3, nHiPass/(double(nNewFs)/2) , 'high'); % high-pass
    vCont  = single(filtfilt(b, a, double(vCont)));
end

% Low-pass filter
nLPPass = nLoPass/(double(nNewFs)/2);
if nLPPass < 1 && nLPPass > 0
    [b,a] = butter(3, nLoPass/(double(nNewFs)/2) , 'low'); % low-pass
    if length(vCont) > (3 * length(a))
        vCont  = single(filtfilt(b, a, double(vCont)));
    end
end

% Rectify
if bRectify, vCont = abs(vCont); end

% Return NaN's where they were removed above
vCont(vNaNIndx) = NaN;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function csSelected = SelectChannels(varargin)
% Select the channels that should be displayed in the GUI
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;

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
            ViewTrialData
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
global g_vSelCh;
if length(csDisplayChannels) == 1
    g_vSelCh = 1;
else
    if isempty(csDisplayChannels), return, end
    hFig = figure;
    nCL = length(csDisplayChannels)*25+35;
    vHW = get(g_hSpike_Viewer, 'position');
    set(hFig, 'position', [vHW(1)+40 (vHW(2)+vHW(4)-nCL-50) 220 nCL+25], 'menu', 'none', 'Name', 'Select channels', 'NumberTitle', 'off')
    for i = 1:length(csDisplayChannels)
        sBoxStr = sprintf('%s (%s)', csChannelDescriptions{i}, csDisplayChannels{i});
        if any(strcmp(FV.csDisplayChannels, csDisplayChannels{i}))
            uicontrol(hFig, 'Style', 'checkbox', 'Position', [10 nCL 170 20], 'String', sBoxStr, 'HorizontalAlignment', 'left', 'backgroundcolor', [.8 .8 .8], 'value', 1);
        else
            uicontrol(hFig, 'Style', 'checkbox', 'Position', [10 nCL 170 20], 'String', sBoxStr, 'HorizontalAlignment', 'left', 'backgroundcolor', [.8 .8 .8], 'value', 0);
        end
        nCL = nCL - 25;
    end
    % 'Select all' checkbox
    uicontrol(hFig, 'Style', 'checkbox', 'Position', [10 nCL 170 20], 'String', 'Select all', ...
        'HorizontalAlignment', 'left', 'backgroundcolor', [.8 .8 .8], 'value', 0, 'callback', ...
        'set(get(gcf,''children''),''value'', 1)');
    nCL = nCL - 25;
    uicontrol(hFig, 'Style', 'pushbutton', 'Position', [50 nCL 100 20], ...
        'Callback', 'global g_vSelCh; g_vSelCh=flipud(get(get(gcf,''children''),''value'')); g_vSelCh=[g_vSelCh{:}];close(gcf)', ...
        'String', '      OK      ' ); % OK button
    uiwait % wait for user to close window or click OK button
    g_vSelCh = g_vSelCh(1:end-2);
end
% Assign channels to FV or just return them
bReturn = 0;
csSelected = csDisplayChannels(find(g_vSelCh));
if ~isempty(varargin)
    if strcmp(varargin{1}, 'return') || isempty(g_vSelCh)
        bReturn = 1;
    end
end
% Return selected channels back to calling function
if bReturn, return
else % Assign channels to FV and update display
    FV.csDisplayChannels = csSelected; % selected channels
    clear global g_vSelCh
    SetStruct(FV)
    ViewTrialData
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OpenFile(varargin)
% OpenFile loads both trial data and the associated settings/results file. This
% function is called from all menu items in the GUI and is a higher-level function
% than LoadTrial. Generally, OpenFile is the only function that calls LoadTrial.
[FV, hWin] = GetStruct;
global g_hSpike_Viewer g_bBatchMode;
persistent p_bNeverApply;

% Decide which file to load next
if isempty(varargin{2})
    % Request which file to load
    if isfield(FV, 'sDirectory')
        [sFile, sPath] = uigetfile( {'*.daq', 'Data Acquisition Toolbox Files (*.daq)'; ...
            '*MergeFile*.spb','Merge files (*MergeFile*.spb)'; ...
            '*.*','All Files (*.*)'}, ...
            'Select data file', CheckFilename([FV.sDirectory '\']));
    else
        [sFile, sPath] = uigetfile( {'*.daq', 'Data Acquisition Toolbox files (*.daq)'; ...
            'MergeFile.spb','Merge files (MergeFile.spb)'; ...
            '*.*','All Files (*.*)'}, ...
            'Select data file');
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
    % Check if filetype is supported
    if all(~strcmp(sFile(end-2:end), {'daq', 'spb'}))
        uiwait(warndlg(sprintf('The filetype .%s is not supported', sFile(end-2:end))))
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
    hFileList = findobj(g_hSpike_Viewer, 'Tag', 'MenuFileList');
    csPaths = flipud(get(get(hFileList, 'Children'), 'UserData'));
    switch varargin{2}
        case -1
            if isempty(csPaths)
                warndlg('No file list has been loaded.', 'Spiky!');
                return
            end
            % Get index of current file
            nNewIndx = find(strcmp(csPaths, FV.sLoadedTrial)) - 1;
            if nNewIndx < 1
                warndlg('You are already at the first file in the current file list.', 'Spiky!')
                return
            end
            sNextTrial = csPaths{nNewIndx};
        case 0
            if isempty(csPaths)
                warndlg('No file list has been loaded.', 'Spiky!');
                return
            end
            nNewIndx = find(strcmp(csPaths, FV.sLoadedTrial)) + 1;
            if nNewIndx > length(csPaths)
                warndlg('You are already at the last file in the current file list.', 'Spiky!')
                return
            end
            sNextTrial = csPaths{nNewIndx};
        otherwise
            if length(csPaths) < varargin{2}
                warndlg('You have reached the last file in the current file list.', 'Spiky!')
                return
            end
            if iscell(csPaths)
                sNextTrial = csPaths{varargin{2}};
            else
                sNextTrial = csPaths(varargin{2},:);
            end
    end
end

% Ask first, then save
if ~isempty(strfind(get(g_hSpike_Viewer, 'Name'), '*')) && ~isempty(FV.sLoadedTrial)
    switch questdlg('Save changes to current file?', 'Spiky!', 'Yes', 'No', 'Cancel', 'Yes')
        case 'Yes', SaveResults
        case 'Cancel', return
    end
end

SetStruct(FV); % update current path

% Exit merge mode if applicable
global g_bMergeMode
if g_bMergeMode, StopMergeMode(); end

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

% Ask to re-use settings
if isempty(FV.sLoadedTrial)
    sAns = 'No';
else
    if isempty(p_bNeverApply), p_bNeverApply = 0; end
    if ~p_bNeverApply && ~g_bBatchMode
        sAns = questdlg('Should current settings be applied to the new file?', 'Spiky!', 'Yes', 'No', 'Never apply', 'No');
    else
        sAns = 'No';
    end
    if strcmpi(sAns, 'Never apply'), p_bNeverApply = 1; sAns = 'No'; end
end

switch sAns
    case 'Yes'
        LoadTrial(sNextTrial); % load new trial data
        [FV, hWin] = GetStruct; % get changed structure (new data/settings)
        % replace default settings with settings of previous file
        csFieldnames = fieldnames(FV_old);
        for fn = 1:length(csFieldnames)
            FV.(csFieldnames{fn}) = FV_old.(csFieldnames{fn});
        end
        SetStruct(FV);
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
        SetStruct(FV);
        % Load trial data. Note that the LoadTrial function will replace
        % the default settings with saved settings/results if they exist
        LoadTrial(sNextTrial);
    otherwise
        return;
end

ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV,hWin] = GetStruct
global g_hSpike_Viewer
hWin = g_hSpike_Viewer;
if ishandle(g_hSpike_Viewer)
    FV = get(g_hSpike_Viewer, 'UserData');
else
    FV = SetFVDefaults();
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetStruct(FV, varargin)
% varargin{1}    'nosaveflag'   Dont reset save flag
global g_hSpike_Viewer
if isempty(g_hSpike_Viewer), g_hSpike_Viewer = findobj('Tag', 'Spiky'); end
set(g_hSpike_Viewer, 'UserData', FV);
bSetSaveFlag = 1;
if ~isempty(varargin)
    if strcmpi(varargin, 'nosaveflag'), bSetSaveFlag = 0; end
end
if bSetSaveFlag
    if length(FV.sLoadedTrial) > 70
        sLoadedTrial = ['...' FV.sLoadedTrial(end-69:end)];
    else sLoadedTrial = FV.sLoadedTrial; end
    set(g_hSpike_Viewer, 'Name', CheckFilename(sprintf('Spiky! - %s*', sLoadedTrial)))
end
% set current working directory to that in active database
if isfield(FV, 'sDirectory')
    cd(FV.sDirectory);
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetGain(varargin)
[FV, hWin] = GetStruct;
if ~CheckDataLoaded, return, end

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
centerfig(hFig, findobj('Tag', 'Spiky'))
set(hFig, 'visible', 'on')

for nCh = 1:length(FV.csChannels)
    % channel name
    sName = FV.csChannels{nCh};
    % append channel description to its name
    if isfield(FV, 'tChannelDescriptions')
        if isfield(FV.tChannelDescriptions, 'sDescription')
            nIndx = find(strcmp({FV.tChannelDescriptions.sChannel}, sName));
            if ~isempty(nIndx)
                if ~isempty(FV.tChannelDescriptions(nIndx).sDescription)
                    sName = [sName ' (' FV.tChannelDescriptions(nIndx).sDescription ')'];
                end
            end
        end
    else
        FV.tChannelDescriptions = struct([]);
    end

    uicontrol(hFig, 'Style', 'text', 'Position', [10 nCL 180 20], 'String', sName, 'HorizontalAlignment', 'left', 'backgroundcolor', [.8 .8 .8]);
    if isfield(FV.tGain, FV.csChannels{nCh})
        nGainThis = FV.tGain.(FV.csChannels{nCh});
    else nGainThis = FV.nGain; end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetMaxSpikeJitter(varargin)
% Set the maximal duration (ms) that two spikes on different channels
% should be considered to be the same same
[FV, hWin] = GetStruct;
FV.nMaxCrossChannelSpikeJitter = GetInputNum('Set maximal spiketime jitter across electrodes (ms)', 'Spiky!', FV.nMaxCrossChannelSpikeJitter);
SetStruct(FV)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EditPreTriggerTime(varargin)
% Change the duration of the signal to keep before threshold-crossings triggers
[FV, hWin] = GetStruct;
FV.nSpikeTrigPreMs = GetInputNum('Set pre-trigger time (ms)', 'Spiky!', FV.nSpikeTrigPreMs);
SetStruct(FV)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetDeadtime(varargin)
% Set the length of time following threshold-crossings events that new
% spikes cannot be detected.
[FV, hWin] = GetStruct;
FV.nDeadTimeMs = GetInputNum(sprintf('Set post-spike deadtime (ms):\n(note that this parameter only affects visually displayed spikes)'), 'Spiky!', FV.nDeadTimeMs);
SetStruct(FV)
ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EditPostTriggerTime(varargin)
% Change the duration of the signal to keep after threshold-crossings triggers
[FV, hWin] = GetStruct;
FV.nSpikeTrigPostMs = GetInputNum('Set post-trigger time (ms)', 'Spiky!', FV.nSpikeTrigPostMs);
SetStruct(FV)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bResult = LoadTrial(snTrial)
% LoadTrial loads data from .daq or .mat file generated by the MAP converter. It
% does NOT load settings from associated .spk files.
% Syntax:   LoadTrial(snTrial) where snTrial is the string name of the file to load
[FV, hWin] = GetStruct;
bResult = 0;

sFile = CheckFilename(snTrial);
if ~exist(sFile, 'file')
    uiwait(warndlg('Data file does not exist'))
    bLoaded = 0;
    return
end

switch lower(sFile(end-3:end))
    case '.mat' % load .mat file
        tData = load(sFile, '-MAT');
    case '.daq' % load .daq file
        try
        [mData, vTime, vAbsTime, tEvents, tDAQInfo] = daqread(CheckFilename(sFile));
        catch
            sStr = sprintf('Error when reading file:\n%s\n\n%s', sFile, lasterr);
            uiwait(warndlg(lasterr, 'modal'));
            sp_disp(sStr)
            return
        end
        mData = single(mData); % convert to single precision to conserve memory
        
        if isempty(mData)
            sStr = 'An error occurred during loading of .DAQ file: File appears to be empty.';
            uiwait(warndlg(sStr, 'modal'))
            sp_disp(sStr)
            return
        end
        tData = struct([]);

        % Import all DAQ channels
        cDaqChannels = {tDAQInfo.ObjInfo.Channel.HwChannel};
        cDaqChannelNames = {tDAQInfo.ObjInfo.Channel.ChannelName};
        for c = 1:length(cDaqChannels)
            % Channel name
            sChName = ['DAQ_' num2str(cDaqChannels{c})]; % use channel's hardware ID
            % Channel description (if alternative channel name exists)
            if isempty(str2num(cDaqChannelNames{c})) % use channel name
                sChDescr = num2str(cDaqChannelNames{c}); % use channel's hardware ID
                if ~isfield(FV, 'tChannelDescriptions')
                    FV.tChannelDescriptions = struct([]);
                end
            else sChDescr = ''; end
            FV.tChannelDescriptions(end+1).sChannel = sChName;
            FV.tChannelDescriptions(end).sDescription = sChDescr;

            % Begin time is counted as the number of seconds elapsed since last midnight
            vTime = tDAQInfo.ObjInfo.InitialTriggerTime;
            vHourToSec = vTime(4)*60*60;
            vMinToSec = vTime(5)*60;
            vSecSinceMidnight = vHourToSec + vMinToSec + vTime(6);

            % Auto-digitize bi-modal signals (putative digital inputs)
            %  1) round
            nMax = max(mData(:,c));
            %  2) normalize range to 0 -> 10
            vNorm = mData(:,c) - min(mData(:,c));
            vNorm = (vNorm / nMax) * 10;
            %  3) round
            vRnd = round( vNorm );
            %  4) get number of unique values. As a heuristic, we assume a
            %     digital signal has less than 5
            nLenU = length(unique(vRnd));
            %  5) get threshold; (max-min) / 2. If signal is a TTL,
            %  assume threshold to be > 1 V
            nThresh = (nMax - min(mData(:,c)) ) / 2;            
            if (nLenU < 5) && nThresh > 1
                %  6) digitize through @DigitizeChannel
                nFs = tDAQInfo.ObjInfo.SampleRate; % Hz
                nTimeBegin = vSecSinceMidnight;
                nTimeEnd = (size(mData, 1) / tDAQInfo.ObjInfo.SampleRate) + vSecSinceMidnight;
                [vUpTimes vDownTimes] = DigitizeChannel(mData(:,c), nThresh, nFs, nTimeBegin, nTimeEnd);

                tData(1).([sChName '_Up']) = vUpTimes; % sec
                tData.([sChName '_Down']) = vDownTimes; % sec

                % Save threshold value used for digitization
                if ~isfield(FV, 'tEventThresholds')
                    FV.tEventThresholds = struct([]);
                    nIndx = 1;
                else nIndx = length(FV.tEventThresholds) + 1; end
                FV.tEventThresholds(nIndx).sChannel = sChName;
                FV.tEventThresholds(nIndx).nThreshold = nThresh;
                
                tData.([sChName '_EventThreshold']) = nThresh; % keep this for backwards compatibility...
                
                FV.csDigitalChannels = unique([FV.csDigitalChannels sChName]);
            end
            
            tData(1).(sChName) = mData(:, c)'; % save raw data
            tData.([sChName '_KHz']) = tDAQInfo.ObjInfo.SampleRate / 1000; % kHz
            tData.([sChName '_KHz_Orig']) = tDAQInfo.ObjInfo.SampleRate / 1000;
            tData.([sChName '_TimeBegin']) = vSecSinceMidnight;
            tData.([sChName '_TimeEnd']) = (size(mData, 1) / tDAQInfo.ObjInfo.SampleRate) + vSecSinceMidnight;
        end

        % Add DAQ events as digital signals
        tEventLog = tDAQInfo.ObjInfo.EventLog;
        for e = 1:length(tEventLog)
            sChName = ['DAQ_' tEventLog(e).Type];
            sChDescr = '';
            if isfield(tData, sChName)
                tData.([sChName '_Up']) = [tData.([sChName '_Up']) tEventLog(e).Data];
                tData.([sChName '_Down']) = [tData.([sChName '_Down']) tEventLog(e).Data];
            else
                tData.([sChName '_Up']) = (tEventLog(e).Data.RelSample / tDAQInfo.ObjInfo.SampleRate) + vSecSinceMidnight; % sec
                tData.([sChName '_Down']) = tData.([sChName '_Up']);
                tData.([sChName '_KHz']) = tDAQInfo.ObjInfo.SampleRate / 1000;
                tData.([sChName '_KHz_Orig']) = tDAQInfo.ObjInfo.SampleRate / 1000;
                tData.([sChName '_TimeBegin']) = vSecSinceMidnight;
                tData.([sChName '_TimeEnd']) = (size(mData, 1) / tDAQInfo.ObjInfo.SampleRate) + vSecSinceMidnight;
                FV.csDigitalChannels = unique([FV.csDigitalChannels sChName]);
            end
            FV.tChannelDescriptions(end+1).sChannel = sChName;
            FV.tChannelDescriptions(end).sDescription = sChDescr;
        end
end

% Store data
FV.tData = tData;
FV.sLoadedTrial = sFile;

% Names of all channels
cFields = fieldnames(FV.tData);
FV.csChannels = {};
for i = 1:length(cFields)
    if ~isempty(strfind(cFields{i}, '_TimeBegin'))
        FV.csChannels{end+1} = cFields{i}(1:end-10);
        if isfield(tData, cFields{i}(1:end-10))
            FV.tYlim(1).(FV.csChannels{end}) = [min(tData.(cFields{i}(1:end-10))) max(tData.(cFields{i}(1:end-10)))];
        else FV.tYlim(1).(FV.csChannels{end}) = []; end % assume its digital
    end
end
FV.vXlim = [];

SetStruct(FV) % save settings obtained thus far

% Load .spb file if it exists
sPath = [FV.sLoadedTrial(1:end-4) '.spb'];
if exist(sPath, 'file')
    OpenSettings([FV.sLoadedTrial(1:end-4) '.spb'])
    [FV, hWin] = GetStruct;
end

% Put back tData as this was removed with the OpenSettings() call above
if isempty(FV.tData)
    FV.tData = tData;
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

if length(sFile) > 70
    sLoadedTrial = ['...' sFile(end-69:end)];
else sLoadedTrial = sFile; end
set(hWin, 'Name', sprintf('Spiky! - %s', sLoadedTrial)) % update window name
bResult = 1;
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AboutSpiky(varargin)
cAboutText = { 'Maintainer:', ...
    ' Per Magne Knutsen <pknutsen@ucsd.edu>', ...n
    ' Department of Physics', ...
    ' University of California, San Diego', ...
    ' 92122 La Jolla, USA', ...
    ' ', ...
    'Contributors:', ...
    ' Per Magne Knutsen, UCSD', ...
    ' Idan Steinberg, Weizmann Institute', ...
    ' Amir Monovitch, Weizmann Institute', ...
    ' ', ...
    'Acknowledgements:', ...
    'Spiky! includes modified code from the open-source Chronux project (http://www.chronux.org) for all of its spike-sorting routines. Chronux is released under the GNU General Public License.', ...
    ' ', ...
    'Musial PG, Baker SN, Gerstein GL, King EA, Keating JG (2002) Signal-to-noise ratio improvement in multiple electrode recording. J Neurosci Methods 115:29-43', ...
    ' ', ...
    'Fee MS, Mitra PP, Kleinfeld D (1996) Automatic sorting of multiple unit neuronal signals in the presence of anisotropic and non-Gaussian variability. J Neurosci Methods 69:175-188', ...
    };
msgbox(cAboutText, 'About Spiky!')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bLoaded = CheckDataLoaded(varargin)
[FV, hWin] = GetStruct;
if ~isfield(FV, 'tData') % check data has been loaded
    uiwait(warndlg('No data has been loaded'))
    bLoaded = 0;
    return
else bLoaded = 1; end % data has been loaded
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetSpikeThreshold(varargin)
[FV, hWin] = GetStruct;
if ~CheckDataLoaded, return, end
zoom off
try
    [nX, nY] = ginput(1);
catch
    return
end
zoom xon
sTag = get(gca, 'Tag');
vContData = FV.tData.(sTag);
[vContData, FV] = AdjustChannelGain(FV, vContData, sTag); % Adjust gain (mV)
if any(strcmp(FV.csEMGChannels, sTag))
    nFs = FV.tData.([sTag '_KHz']) * 1000; % channel sampling rate
    [vContData, vTime, nFs] = FilterChannel(vContData, [], nFs, FV.nEMGLoPass, FV.nEMGHiPass, FV.nEMGRectify, 'nodecimate');
end
nMed = median(vContData);
if nY == 0 || nY == nMed
    % Threshold line can't be zero!
    warndlg('Threshold cannot be zero.')
    return;
elseif nY < nMed
    % Negative threshold line
    if isfield(FV.tSpikeThresholdsNeg, sTag) && isfield(FV.tSpikes, sTag)
        % Ask whether to discard currently detected spikes
        if strcmp(questdlg('Replacing the current threshold line will discard all detected spikes. Continue?', 'Spiky!', 'Yes', 'No', 'No'), 'No')
            return;
        else FV.tSpikes = rmfield(FV.tSpikes, sTag); end
    end
    FV.tSpikeThresholdsNeg(1).(sTag) = nY; % mV
else
    % Positive threshold line
    if isfield(FV.tSpikeThresholdsPos, sTag) && isfield(FV.tSpikes, sTag)
        % Ask whether to discard currently detected spikes
        if strcmp(questdlg('Replacing the current threshold line will discard all detected spikes. Continue?', 'Spiky!', 'Yes', 'No', 'No'), 'No')
            return;
        else FV.tSpikes = rmfield(FV.tSpikes, sTag); end
    end
    FV.tSpikeThresholdsPos(1).(sTag) = nY; % mV
end
SetStruct(FV)
ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RemoveSpikeThreshold(sCh, sDir)
[FV, hWin] = GetStruct;
sTag = get(gco, 'Tag');
sCh = sTag(1:end-4);

% Ask whether to discard currently detected spikes
if isfield(FV.tSpikes, sCh)
    if strcmp(questdlg('Removing this threshold line will discard all detected spikes on this channel. Continue?', 'Spiky!', 'Yes', 'No', 'No'), 'No')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetThresholdManually(sCh, sDir)
[FV, hWin] = GetStruct;
sTag = get(gco, 'Tag');
sCh = sTag(1:end-4);

% Ask whether to discard currently detected spikes
if isfield(FV.tSpikes, sCh)
    if strcmp(questdlg('Adjusting the threshold manually will erase all detected spikes and sorting data on this channel. Continue?', 'Spiky!', 'Yes', 'No', 'No'), 'No')
        return;
    else FV.tSpikes = rmfield(FV.tSpikes, sCh); end
end

% Ask for new threshold value
switch lower(sTag(end-2:end))
    case 'neg'
        cAns = inputdlg('New threshold value:', 'Spiky!', 1, {num2str(FV.tSpikeThresholdsNeg.(sCh))});
        if isempty(cAns), return, end
        FV.tSpikeThresholdsNeg.(sCh) = str2num(cAns{1});
    case 'pos'
        cAns = inputdlg('New threshold value:', 'Spiky!', 1, {num2str(FV.tSpikeThresholdsPos.(sCh))});
        if isempty(cAns), return, end
        FV.tSpikeThresholdsPos.(sCh) = str2num(cAns{1});
    otherwise
        return;
end
SetStruct(FV)
ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AutoDetectThresholds(varargin)
% Automatically select 4 x STD as negative spike threshold
if ~CheckDataLoaded, return, end
[FV,hWin] = GetStruct;

% Select threshold to find thresholds on
csSelected = FV.csDisplayChannels;
hWait = waitbar(0, 'Automatically determining spike thresholds...');

for nCh = 1:length(csSelected)
    hWait = waitbar(nCh/length(csSelected), hWait);
    sCh = csSelected{nCh};
    vCont = FV.tData.(sCh);
    [vCont, FV] = AdjustChannelGain(FV, vCont, sCh); % Adjust gain (mV)

    % Filter and rectify EMG signals
    nFs = FV.tData.([sCh '_KHz']) * 1000; % channel sampling rate
    if any(strcmp(FV.csEMGChannels, sCh))
        [vCont, vTime, nFs] = FilterChannel(vCont, [], nFs, FV.nEMGLoPass, FV.nEMGHiPass, FV.nEMGRectify, 'nodecimate');
    else
        [vCont, vTime, nFs] = FilterChannel(vCont, [], nFs, 6000, 600, 0, 'nodecimate');
    end
    
    %vThresholds = prctile(vCont, [1 99]);
    % Negative threshold (values < 1th percentile)
    %FV.tSpikeThresholdsNeg(1).(sCh) = vThresholds(1);
    FV.tSpikeThresholdsNeg(1).(sCh) = nanmean(vCont) - (4*nanstd(vCont));
    % Positive threshold (values > 99th percentile)
    %FV.tSpikeThresholdsPos(1).(sCh) = vThresholds(2);
end
close(hWait)
SetStruct(FV)
ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DetectSpikes(varargin)
% Note: This function supports both positive and negative thresholds, while
% GetWaveforms currently does not!!   24/02/2007
if ~CheckDataLoaded, return, end
global g_bMergeMode
[FV, hWin] = GetStruct;

% Iterate across all channels with threshold lines(s)
csChNeg = fieldnames(FV.tSpikeThresholdsNeg);
csChPos = fieldnames(FV.tSpikeThresholdsPos);
csChUnique = unique([csChNeg;csChPos]);
if length(csChUnique) > 1
    if ~g_bMergeMode
        bResult = SpikyWaitbar(0, length(csChUnique));
    end
end
for nCh = 1:length(csChUnique)
    sCh = csChUnique{nCh};
    % Negative threshold line
    if isfield(FV.tSpikeThresholdsNeg, sCh), vThresh(1) = FV.tSpikeThresholdsNeg.(sCh);
    else vThresh(1) = NaN; end
    % Positive threshold line
    if isfield(FV.tSpikeThresholdsPos, sCh), vThresh(2) = FV.tSpikeThresholdsPos.(sCh);
    else vThresh(2) = NaN; end
    % Filter continuous signal
    vCont = FV.tData.(sCh);
    [vCont, FV] = AdjustChannelGain(FV, vCont, sCh); % Adjust gain (mV)
    nFs = FV.tData.([sCh '_KHz']) * 1000; % channel sampling rate
    % Use custom filter, if applicable
    vContRaw = vCont;
    if any(strcmp(FV.csEMGChannels, sCh))
        [vCont, vTime, nFs] = FilterChannel(vCont, [], nFs, FV.nEMGLoPass, FV.nEMGHiPass, FV.nEMGRectify, 'nodecimate');
    else
        % Use default filter (0.6  - 6 kHz)
        [vCont, vTime, nFs] = FilterChannel(vCont, [], nFs, 6000, 600, 0, 'nodecimate');
    end
    % Run spike detection
    nBeginTime = FV.tData.([sCh '_TimeBegin']) * 1000; % begin time
    tSpikes = GetWaveforms(vCont, vContRaw, vThresh, nFs, FV.nSpikeTrigPreMs, FV.nSpikeTrigPostMs, 0.01, nBeginTime);
    FV.tSpikes(1).(sCh) = tSpikes;
    FV.tSpikes.(sCh).dejittered = 0;
    if length(csChUnique) > 1
        if ~g_bMergeMode
            bResult = SpikyWaitbar(nCh, length(csChUnique));
            if ~bResult
                bResult = SpikyWaitbar(length(csChUnique), length(csChUnique));
                return;
            end
        end
    end
end
if length(csChUnique) > 1
    if ~g_bMergeMode
        bResult = SpikyWaitbar(length(csChUnique), length(csChUnique));
    end
end
SetStruct(FV)
ViewTrialData();
global g_bBatchMode
if ~g_bBatchMode, BatchRedo([], 'DetectSpikes'); end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tSpikes = GetWaveforms(vTrace, vRawTrace, vThresh, nFS, nPreTrigMs, nPostTrigMs, nDeadTimeMs, nBeginTime)
% vTrace        Continuous trace used for spike detection
% vRawTrace     Continuous trace used for spike extraction (waveforms)
% vThresh       Threshold [NEG POS]
%               PS: Currently only negative thresholds are supported!
% nFs           Sampling rate
% nPreTrigMs    Time to keep before spike onset (ms)
% nPostTrigMs   Time to keep after spike onset (ms)
% nDeadTimeMs   Time after spike when no other spike can be detected (ms)
% hWait         Absolute time of first sample (ms)

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

vSpiketimes = sort([vSpiketimes_Neg; vSpiketimes_Pos]);

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
mWaveforms = vRawTrace(mIndx); % waveforms
tSpikes.waveforms = double(mWaveforms); % spikes must of double() type (required by sorting functions)
% Keep spiketimes after adjusting by nBeginTime
nBeginTime = nBeginTime * (nFS/1000); % samples
tSpikes.spiketimes = vSpiketimes + nBeginTime;
tSpikes.Fs = nFS;
tSpikes.threshV = [vThresh(1) inf]; % thresholds [hi lo]
tSpikes.threshT = nPreTrig + 2;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DejitterSpikes(varargin)
if ~CheckDataLoaded, return, end
[FV,hWin] = GetStruct;
csFieldnames = fieldnames(FV.tSpikes);
if length(csFieldnames) > 1
    hWait = waitbar(0, 'Dejittering spikes ...');
end
for i = 1:length(csFieldnames)
    if length(csFieldnames) > 1, waitbar(i/length(csFieldnames), hWait); end
    % Dont try to dejitter stereo/tetrodes
    if ~isempty(strfind(csFieldnames{i},'__')), continue, end
    %if isfield(FV.tSpikes.(csFieldnames{i}), 'dejittered')
    %    if FV.tSpikes.(csFieldnames{i}).dejittered, continue, end
    %end
    if isempty(FV.tSpikes.(csFieldnames{i}).waveforms), continue, end
    tSpikes = FV.tSpikes.(csFieldnames{i});
    % in case several thresholds are given (eg. for merged files) check
    % they are all the same; if not, dejittering cannot be made
    if length(unique(tSpikes.threshV(:,1))) > 1
        sStr = ['Channel ' csFieldnames{i} ' contains more than one threshold value. Therefore, the channel cannot be dejittered.'];
        uiwait(warndlg(sStr))
        continue
    end
    tSpikes.threshV = tSpikes.threshV(1,:);
    tSpikes.threshT = tSpikes.threshT(1);
    try
        tSpikes = ss_dejitter(tSpikes, round(tSpikes.Fs(1)/10000)*2); % dejitter spikes
    catch
        sStr = sprintf('Failed to dejitter channel %s. The error returned was:\n\n%s', csFieldnames{i}, lasterr);
        hFig = warndlg(sStr);
        uiwait(hFig)
        continue
    end
    % Sometimes ss_dejitter fails; if so, show an error msg
    if isempty(tSpikes.waveforms)
        sStr = ['Failed to dejitter channel ' csFieldnames{i}];
        hFig = warndlg(sStr);
        uiwait(hFig)
    else % Dejittering succeeded
        FV.tSpikes.(csFieldnames{i}) = tSpikes;
        FV.tSpikes.(csFieldnames{i}).dejittered = 1;
        sStr = ['Dejittered channel ' csFieldnames{i}];
    end
    sp_disp(sStr)
end
if length(csFieldnames) > 1, close(hWait); end
SetStruct(FV)
global g_bBatchMode
if ~g_bBatchMode BatchRedo([], 'DejitterSpikes'); end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RemoveOutlierSpikes(varargin)
% Remove waveforms that are considered outliers using K-means based outlier
% detection algorithm
if ~CheckDataLoaded, return, end
[FV,hWin] = GetStruct;
csFieldnames = fieldnames(FV.tSpikes);
if length(csFieldnames) > 1,
    hWait = waitbar(0, 'Removing outliers from all channels...');
    set(hWait, 'position', get(hWait, 'position')+[0 150 0 0])
end
for nCh = 1:length(csFieldnames)
    if length(csFieldnames) > 1, waitbar(nCh/length(csFieldnames), hWait), end
    tSpikes = FV.tSpikes.(csFieldnames{nCh});
    tSpikes = ss_outliers(tSpikes, (1-1/(size(tSpikes.waveforms,2)*100))); % remove outliers
    FV.tSpikes.(csFieldnames{nCh}) = tSpikes;
end
if length(csFieldnames) > 1, close(hWait), end
SetStruct(FV)
global g_bBatchMode
if ~g_bBatchMode BatchRedo([], 'RemoveOutlierSpikes'); end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowAggregationTree(varargin)
% Show overclustered assingments and results after joining clusters
% TODO: Needs to be modified if users has joined/split clusters manually
if ~CheckDataLoaded, return, end
[FV,hWin] = GetStruct;
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
set(hFig, 'color', [.2 .2 .2], 'name', 'Spiky! Aggregration Tree', 'NumberTitle', 'off', 'menubar', 'none')
aggtree(FV.tSpikes.(sCh));

hHeader = header(['Channel ' sCh ' - Spike Cluster Aggregation Tree'], 12);
set(hHeader, 'color', 'w', 'interpreter', 'none')

SetStruct(FV)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bSorted = CheckIfSorted(varargin)
[FV,hWin] = GetStruct;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowSpikeClusters(varargin)
[FV,hWin] = GetStruct;
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
for u = 1:nUnits % one row per unit
        
    % open new window if 4 units/rows have been reached
    if p == 10
        hMainWin = figure;
        set(hMainWin, 'Color', [.2 .2 .2], 'name', 'Spiky! Spike Clusters', 'NumberTitle', 'off')
        p = 1;
        % scale window width to number of columns
        vPos = get(hMainWin, 'position');
        if nCols == 10, set(hMainWin, 'position', vPos .* [1 1 2 1]), end
        centerfig(hMainWin)
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
    vIndxKeep = randperm(length(vIndx));
    vIndxKeep = vIndxKeep(1:min([1000 length(vIndx)]));
    
    % Always show spikes with extreme amplitudes
    vSums = sum(abs(mWaveforms(vIndx,:)),2);
    vIndxKeep = [vIndxKeep find(vSums > nanmean(vSums)+(3*nanstd(vSums)))'];
    
    %   1) Plot waveforms of unit
    nX = .05; nY = .08;
    nW = .90; nH = .8;
    nSpace = 0;
    axes('position', [nX+(p-1)*(nW/nCols)+nSpace nY+(nH/3)*2 nW/nCols-nSpace nH/3]);
    set(gca, 'uicontextmenu', hMenu(end))
    
    vYY = mWaveforms(vIndx(vIndxKeep),:)';
    vXX = repmat([([1:size(vYY,1)]/nFs)*1000]', 1, size(vYY, 2));
    if vUnits(u) == 0, mCol = [.5 .5 .5]; % outlier
    else mCol = FV.mColors(u,:); end
    plot(vXX, vYY, '-', 'color', mCol, 'linewidth', .5)
    axis tight
    set(gca, 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'fontsize', 6)
    title(['Unit# ' num2str(vUnits(u)) ' - ' num2str(length(vIndx)) ' spikes'], 'fontsize', 6, 'color', [.6 .6 .6])
    if p>1, set(gca, 'xtick', [])
    else xlabel('ms'); ylabel('mV'); end

    %   2) Plot 2D histogram of unit
    axes('position', [nX+(p-1)*(nW/nCols)+nSpace nY+nH/3 nW/nCols-nSpace nH/3]);
    set(gca, 'uicontextmenu', hMenu(end))    
    hist2d(mWaveforms(vIndx,:),200);
    %colormap(pink(256))
    colormap(hot(256))
    set(gca, 'xtick', [], 'ytick', [], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6])

    %   3) ISI histogram
    hAx = axes('position', [nX+(p-1)*(nW/nCols)+nSpace nY nW/nCols-nSpace nH/3]);
    set(hAx, 'uicontextmenu', hMenu(end))
    vTheseSpiketimes = DropDeadtimeSpikes((vSpiketimes(vIndx) ./ nFs)) .*1000 ;
    nMaxISI = 50; % ms
    nSpkLen = length(find(diff(vTheseSpiketimes) <=50 ));
    if nSpkLen > 2500, nStep = 0.25;
    elseif nSpkLen > 500, nStep = 0.5;
    else nStep = 1; end
    vRange = 0:nStep:nMaxISI;
    [vN, vX] = hist(diff(vTheseSpiketimes), [vRange vRange(end)+1]); % msec, linear scale
    vN = vN(1:end-1); % spikes/bin
    vX = vX(1:end-1);
    bar(vX, vN, 'facecolor', [.85 .85 .85], 'edgecolor', [.85 .85 .85])
    hold on
    % Plot deadtime
    plot([FV.nDeadTimeMs FV.nDeadTimeMs], [0 max(vN)+10], 'r:')
    set(hAx, 'Color', [.1 .1 .1], 'xtick', [0:10:vRange(end)], 'xlim', [0 nMaxISI/2], ...
        'ylim', [0 max(vN)+1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'fontsize', 6)
    if p==1 ylabel('N'); end
    xlabel('ISI (ms)')
    hHeader = header(['Clustered spikes: ' GetChannelDescription(sCh)], 12);
    set(hHeader, 'color', 'w', 'interpreter', 'none')
    % Load channel pop-down meny
    csDescriptions = cell(size(csFieldnames));
    for c = 1:length(csFieldnames)
        csDescriptions{c} = GetChannelDescription(csFieldnames{c});
    end
    
    uicontrol(hMainWin, 'Style', 'popupmenu', 'units', 'normalized', ...
        'Position', [0 .95 .2 .05], 'String', ...
        csDescriptions, 'userdata', csFieldnames, 'callback', ...
        'cUD=get(gcbo,''userdata'');nV=get(gcbo,''value'');close(gcf);spiky([''ShowSpikeClusters([],'' '''''''' cUD{nV} '''''''' '')''])')
    drawnow
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowOverlappingSpikeClusters(varargin)
if ~CheckDataLoaded, return, end
[FV,hWin] = GetStruct;
if ~isempty(varargin)
    if isnumeric(varargin{1}), sCh = '';
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
sFigTitle = 'Spiky! Cluster Control';
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
set(hMainWin, 'Color', [.2 .2 .2], 'Tag', sFigTag, 'name', sFigTitle, 'NumberTitle', 'off', 'MenuBar', 'figure')
hold on
set(hAx, 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'Tag', sAxTag)
drawnow
% Remove all checkboxes from window
delete(findobj(hMainWin, 'Style', 'checkbox'));
for u = 1:length(vUnits) % one row per unit
    % Line color
    if vUnits(u) == 0, mCol = [.2 .2 .2]; mColTxt = [.4 .4 .4]; % outliers
    else mCol = FV.mColors(u,:); mColTxt = FV.mColors(u,:); end

    if ~ishandle(hMainWin), return, end % exit if window was closed
    % Create new line objects for this unit if none already exist
    sLineTag = num2str(vUnits(u));
    hLines = findobj(hAx, 'Tag', sLineTag);
    if isempty(hLines)
        vIndx = find(tSpikes.hierarchy.assigns == vUnits(u));
        vIndxKeep = randperm(length(vIndx)); % keep just 1000 waveforms and drop the rest
        vIndxKeep = vIndxKeep(1:min([2000 length(vIndx)]));
        
        % Always show spikes with extreme amplitudes
        vSums = sum(abs(mWaveforms(vIndx,:)),2);
        vIndxKeep = [vIndxKeep find(vSums > nanmean(vSums)+(3*nanstd(vSums)))'];
        
        vYY = mWaveforms(vIndx(vIndxKeep),:)';
        vXX = repmat([([1:size(vYY,1)]/nFs)*1000]', 1, size(vYY, 2));
        if ~ishandle(hMainWin), return, end % exit if window was closed
        figure(hMainWin);
        hLines = plot(vXX, vYY, '-', 'color', mCol, 'linewidth', .5);
        set(hLines, 'tag', sLineTag);
    else
        set(hLines, 'color', mCol)
    end

    if vUnits(u) == 0, sUnitName = 'Outliers';
    else sUnitName = ['Unit ' num2str(vUnits(u))]; end
    hCheckbox = uicontrol(hMainWin, 'Style', 'checkbox', 'units', 'normalized', 'Position', [.85 .85-([u-1]*.04) .15 .05], 'String', sUnitName, 'HorizontalAlignment', 'left', 'backgroundcolor', [.2 .2 .2], 'value', 1, 'foregroundColor', mColTxt);
    set(hCheckbox, 'Tag', num2str(vUnits(u)), 'callback', @ToggleWaveforms);
    switch get(hLines(1), 'visible')
        case 'on'
            set(hCheckbox, 'value', 1)
        case 'off'
            set(hCheckbox, 'value', 0)
    end

    xlabel('ms'); ylabel('mV')
    hold on
    drawnow
end

% Plot threshold line(s)
cFields = fieldnames(FV.tSpikeThresholdsNeg);
nIndx = find(strcmpi(cFields, sCh));
if ~ishandle(hAx) return; end
if ~isempty(FV.tSpikeThresholdsNeg) && ~isempty(fieldnames(FV.tSpikeThresholdsNeg)) % negative threshold
    nT_neg = FV.tSpikeThresholdsNeg.(cFields{nIndx});
    plot(get(hAx, 'xlim'), [nT_neg nT_neg], 'w--')
end
if ~isempty(FV.tSpikeThresholdsPos) && ~isempty(fieldnames(FV.tSpikeThresholdsPos)) % positive threshold
    nT_pos = FV.tSpikeThresholdsPos.(cFields{nIndx});
    plot(get(hAx, 'xlim'), [nT_pos nT_pos], 'w--')
end

if ~ishandle(hMainWin) return, end

uicontrol(hMainWin, 'Style', 'popupmenu', 'units', 'normalized', ...
    'Position', [.1 .85 .2 .05], 'String', ...
    {'Reassign spike' 'Reassign group' 'Assign as outlier' 'Merge clusters' 'Merge visible' 'Split cluster' 'Undo last action'}, 'callback', @ClusterControl)

uicontrol(hMainWin, 'Style', 'text', 'units', 'normalized', ...
    'Tag', 'StatusField', 'Position', [.31 .855 .53 .04], 'backgroundcolor', [.1 .1 .1], ...
    'foregroundcolor', 'w', 'horizontalalignment', 'left', 'visible', 'off')

hHeader = header(['Clustered spikes: ' GetChannelDescription(sCh)], 12);
set(hHeader, 'color', 'w', 'interpreter', 'none')

% Pop-up menu
csDescriptions = cell(size(csFieldnames));
for c = 1:length(csFieldnames)
    csDescriptions{c} = GetChannelDescription(csFieldnames{c});
end

uicontrol(hMainWin, 'Style', 'popupmenu', 'units', 'normalized', ...
    'Position', [0 .95 .2 .05], 'String', ...
    csDescriptions, 'userdata', csFieldnames, 'callback', ...
    'cUD=get(gcbo,''userdata'');nV=get(gcbo,''value'');close(gcf);spiky([''ShowOverlappingSpikeClusters('' '''''''' cUD{nV} '''''''' '')''])')

% Set the click Callback function
hCursor = datacursormode(hMainWin);
datacursormode on
set(hCursor, 'UpdateFcn', @ShowCurrentLine, 'SnapToDataVertex', 'on');
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClusterControl(varargin)
[FV,hWin] = GetStruct;
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
            nSpikeCluster = str2num(get(hLine, 'tag')); % selected spike
            vYdata = get(hLine, 'ydata');
            nNewCluster = 0;
            if strcmp(questdlg(sprintf('Move spike from unit %d to unit %d?', nSpikeCluster, nNewCluster), 'Spiky!', 'OK', 'Cancel', 'OK'), 'OK')
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
            if strcmp(questdlg(sprintf('Move spike from unit %d to unit %d?', nSpikeCluster, nNewCluster), 'Spiky!', 'OK', 'Cancel', 'OK'), 'OK')
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
        g_lin = plot([g_pnt1(3) g_pnt1(3)], [g_pnt1(3) g_pnt1(3)], 'w');
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
            % Decimate waveform (reduces execution time)
            %vY_spk = interp1(vX_spk, vY_spk, vX_spk(1:10:end),'nearest');
            %vX_spk = vX_spk(1:10:end);
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
                    PInt = inv(A)*b; % intersection point
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
        % Highlight intersecting waveforms
        if ~isempty(hSelectedHandles)
            set(hSelectedHandles, 'marker', 'square', 'markersize', 4, 'MarkerEdgeColor', [1 1 1], 'LineWidth', 2)
            vSpikeClusters = str2num(cell2mat(get(hSelectedHandles, 'tag'))); % unit IDs of selected spikes
            set(hStatus, 'string', 'Select a spike from the new cluster', 'visible', 'on')
            waitfor(figure('visible', 'off', 'Tag', 'UIHoldFigure')) % Create an invisible figure
            hLine = findobj(gcf, 'type', 'line', 'marker', 's');
            nNewCluster = str2num(get(hLine, 'tag'));
            if strcmp(questdlg(sprintf('Reassign selected spikes from units [ %s] to unit %d?', sprintf('%d ', unique(vSpikeClusters)), nNewCluster), 'Spiky!', 'OK', 'Cancel', 'OK'), 'OK')
                for h = 1:length(hSelectedHandles)
                    vYdata = get(hSelectedHandles(h), 'ydata');
                    vIndx = find(FV.tSpikes.(sChannel).hierarchy.assigns == vSpikeClusters(h));
                    mWaveforms = FV.tSpikes.(sChannel).waveforms(vIndx, :);
                    mYdata = repmat(vYdata, size(mWaveforms,1), 1);
                    nIndx = vIndx(find(all(mYdata == mWaveforms, 2)));
                    FV.tSpikes.(sChannel).hierarchy.assigns(nIndx) = nNewCluster;
                end
                SetStruct(FV)
                close(hFig)
                ShowOverlappingSpikeClusters(sChannel);
                return
            end
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
        if nSpikeClusterFrom == nSpikeClusterTo, warndlg('You cant merge a cluster into itself.', 'Spiky!'); return, end

        if nSpikeClusterTo == 0 || nSpikeClusterFrom == 0
            warndlg('You cannot merge outliers with sorted clusters. You can only reassign individual outliers to existing clusters.')
            return
        end

        if strcmp(questdlg(sprintf('Merge cluster %d into %d?', nSpikeClusterFrom, nSpikeClusterTo), 'Spiky!', 'OK', 'Cancel', 'OK'), 'OK')
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

        if strcmp(questdlg(sprintf('Merge all visible clusters into cluster %d?', nSpikeClusterTo), 'Spiky!', 'OK', 'Cancel', 'OK'), 'OK')
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
        end
        
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
        elseif strcmp(questdlg(sprintf('Split cluster %d into its parent clusters %d and %d?', nSpikeClusterSplit, nSpikeClusterSplit, mTree(vRowIndx(end), 2)), 'Spiky!', 'OK', 'Cancel', 'OK'), 'OK')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ToggleWaveforms(varargin)
% Toggle visibility of a selected group of waveforms
hCallbackHandle = gco; % callback handle
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
    if ~iscell(cYData), cYData = mat2cell(cYData); end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SortSpikes(varargin)
% Sort spikes on selected channel
if ~CheckDataLoaded, return, end
[FV,hWin] = GetStruct;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SortAllChannels(varargin)
% Run spike sorting on all thresholded channels
if ~CheckDataLoaded, return, end
[FV,hWin] = GetStruct;
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
    tSpikes.Fs = tSpikes.Fs(1);
    % Over-clustering
    tSpikes = ss_kmeans(tSpikes);
    % Aggregation
    tSpikes = ss_energy(tSpikes);
    tSpikes = ss_aggregate(tSpikes);
    % Save results
    FV.tSpikes.(sCh) = tSpikes;
    SetStruct(FV)
end
if length(csChannels) > 1, close(hWait); end
ShowOverlappingSpikeClusters;
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sCh, bResult] = SelectChannelNumber(varargin)
if ~CheckDataLoaded, return, end
sCh = []; bResult = 0;
[FV,hWin] = GetStruct;

% Strings representing all individual and stereo/tetrode electrodes
sTit = 'Select channel';
if isempty(varargin)
    csChannels = FV.csChannels;
else
    csChannels = varargin{1};
    if length(varargin) > 1
        sTit = varargin{2};
    end
end
csChannelNames = csChannels;

% Get channel descriptions, if available
csDescriptions = {};
for ch = 1:length(csChannels)
    csDescriptions{ch} = GetChannelDescription(csChannels{ch});
end

% Append channel descriptions to channel DAQ names
for c = 1:length(csChannels)
    sDAQ = csChannels{c};
    sDescr = csDescriptions{c};
    if ~isempty(sDescr)
        csChannels{c} = [sDAQ '  (' sDescr ')'];
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
        'Name', sTit, 'NumberTitle', 'off', 'CloseRequestFcn', 'delete(gcf); global g_bClosePressed; g_bClosePressed = 1;')
    centerfig(hFig, hWin)
    uicontrol(hFig, 'Style', 'popupmenu', 'Position', [10 5 300 20], 'String', csChannels, 'Callback', 'global g_nCh; g_nCh=get(gco, ''value'');', 'value', 1);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ViewWaveforms(varargin)
if ~CheckDataLoaded, return, end
[FV,hWin] = GetStruct;
nSubplots = 0; % get number of channels for which spikes exist
csFieldnames = fieldnames(FV.tSpikes);
for i = 1:length(csFieldnames)
    if ~isempty(FV.tSpikes.(csFieldnames{i})), nSubplots = nSubplots + 1; end;
end
if length(csFieldnames)==0
    hWarn = warndlg('No waveforms have been extracted on any channel. You must run spike detection before viewing waveforms.');
    uiwait(hWarn); return
end
hFig = figure; colormap(pink(256));
set(hFig, 'color', [.2 .2 .2], 'name', 'Spiky! Waveforms', 'NumberTitle', 'off', 'menubar', 'none')
p = 1;
% Find how many channels should be shown
nChTot = 0;
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
    subplot(2,nSubplots,p); %p = p+1;
    vXX = repmat([([1:size(FV.tSpikes.(csFieldnames{i}).waveforms,2)]/nFs)*1000]', 1, size(FV.tSpikes.(csFieldnames{i}).waveforms, 1));
    vYY = FV.tSpikes.(csFieldnames{i}).waveforms';
    plot(vXX, vYY, ':', 'linewidth', .5, 'color', [.7 .7 .7]);
    title([sCh ' - ' num2str(size(vXX,2)) ' spikes'], 'interpreter', 'none', 'color', 'w');
    ylabel('Amplitude')
    axis tight;
    set(gca, 'xlim', [0 max(vXX(:,1))], 'xtick', 0:.5:max(vXX(:,1)), 'xcolor', [.7 .7 .7], 'ycolor', [.7 .7 .7], 'color', [.1 .1 .1]);
    subplot(2,nSubplots,p+nSubplots); p = p + 1;
    hist2d(FV.tSpikes.(csFieldnames{i}).waveforms);
    set(gca, 'xtick', [], 'xcolor', [.7 .7 .7], 'ycolor', [.7 .7 .7]);
    ylabel('Amplitude (mV)')
    xlabel('ms')
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ViewWaveformPCs(varargin)
if ~CheckDataLoaded, return, end
[FV,hWin] = GetStruct;
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
if any(unique(FV.tSpikes.(sCh).hierarchy.assigns) == 0) % outliers exist
    vUnique = unique(FV.tSpikes.(sCh).hierarchy.assigns);
    vIndx = find(vUnique > 0);
    nIndx = find(vUnique == 0);
    mCols = [.3 .3 .3];
    mCols(vIndx, 1:3) = FV.mColors(2:length(vIndx)+1,:);
else % no outliers exist
    mCols = FV.mColors(1:length(unique(FV.tSpikes.(sCh).hierarchy.assigns)),:);
end

%ssg_databrowse3d(FV.tSpikes.(sCh), mCols, FV.tSpikes.(sCh).hierarchy.assigns); % old call
mColsTemp = zeros(max(FV.tSpikes.(sCh).hierarchy.assigns), 3);
mColsTemp(unique(FV.tSpikes.(sCh).hierarchy.assigns), :) = mCols;
FV.tSpikes.(sCh).overcluster.colors = mColsTemp(1:end,:);
ssg_databrowse3d(FV.tSpikes.(sCh), FV.tSpikes.(sCh).hierarchy.assigns, unique(FV.tSpikes.(sCh).hierarchy.assigns))
set(gcf, 'color', [.2 .2 .2], 'Name', 'Spiky! Principal Components', 'NumberTitle', 'off', 'ToolBar', 'figure')
hAx = get(gcf, 'Children');
set(hAx(strcmp(get(hAx, 'type'), 'axes')), 'color', [.1 .1 .1])
hHeader = header(['Channel ' sCh ' - Click in axis labels to change scale'], 10);
set(hHeader, 'color', 'w', 'interpreter', 'none')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function KeyboardMode(varargin)
[FV, hWin] = GetStruct;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AddTetrode(varargin)
% Add a multi-trode (2, 3 or 4 wires)
% AddTetrode(sCh1__2) will use the channel names supplied
% AddTetrode() will ask user which channels to combine
if ~CheckDataLoaded, return, end
global g_vElectrodes
g_vElectrodes = [1 1 1 1];
[FV, hWin] = GetStruct;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RemoveTetrode(varargin)
[FV,hWin] = GetStruct;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FV = SetFVDefaults()
FV.sLoadedTrial = '';
FV.nGain = 1;                   % default hardware gain
FV.csDisplayChannels = {};      % channels to be displayed (analogue, continuous)
FV.csChannels = {};             % names of all channels in datafile
FV.csDigitalChannels = {};      % list of digital channels
FV.mColors = [.1 .1 .9; .1 .75 .1; .9 .1 .1; .1 .75 .75; .75 .1 .75; .75 .75 .1; 1 .5 .25; .6 .5 .4; .9 .9 .1];
FV.mColors = [FV.mColors; (FV.mColors.^2).^2; sqrt(sqrt(FV.mColors)); sqrt(FV.mColors); FV.mColors.^2; sqrt(sqrt(sqrt(FV.mColors))); ((FV.mColors.^2).^2).^2];
FV.tSpikeThresholdsNeg = struct([]);
FV.tSpikeThresholdsPos = struct([]);
FV.tGain = struct([]);
FV.tYlim = struct([]);
FV.tSpikes = struct([]);
FV.nSpikeTrigPreMs = 0.75; % ms
FV.nSpikeTrigPostMs = 1.75; % ms
FV.vXlim = [];
FV.nDeadTimeMs = 0.01; % ms
FV.nEMGHiPass = 300;    % high-pass frequency (Hz)
FV.nEMGLoPass = 10000;  % low-pass frequency (Hz)
FV.nEMGRectify = 0;     % rectify data (1=yes, 0=no)
FV.csEMGChannels = {};
FV.bPlotRasters = 0;
FV.bPanOn = 0;
FV.nMaxCrossChannelSpikeJitter = 0.5; % maximal allowed spiketimmer jitter across electrodes (ms)
% Experiment variables
FV.tExperimentVariables = struct([]);
FV.tExperimentVariables(1).sVariable = 'sAnimal';
FV.tExperimentVariables(1).sValue = '000000';
FV.tExperimentVariables(2).sVariable = 'sDate';
FV.tExperimentVariables(2).sValue = 'dd-mm-yyyy';
FV.tExperimentVariables(3).sVariable = 'sExperiment';
FV.tExperimentVariables(3).sValue = 'undefined';
g_bBatchMode = false;
FV.sBatchLock = 'off';
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DistributeSettings(varargin)
% Save the current settings as settings for all other files in current directory.
% Overwrite settings of other files if they exist by default
[FV,hWin] = GetStruct;
if ~CheckDataLoaded, return, end
global g_hSpike_Viewer

% Get current file list
hFileList = findobj(g_hSpike_Viewer, 'Tag', 'MenuFileList');
csPaths = flipud(get(get(hFileList, 'Children'), 'UserData'));
if isempty(csPaths)
    warndlg('No additional files have been loaded into your workspace (see File->Files). You first need to load a directory, a directory tree or a custom list of files into your workspace before distributing settings and merging results.', 'Spiky!');
    return
end

if ischar(csPaths)
    uiwait(warndlg('Cannot distribute settings since there is only a single file loaded.', 'Spiky!'))
    return
end

% Warn that this action will overwrite data of other files
%switch questdlg('The current settings will be distributed to all files loaded into your workspace (File->Files). This action will overwrite existing settings and sorting results of other files. Imported fields will be kept. Are you sure you want to continue?', 'Spiky!', 'Yes', 'No', 'Cancel', 'Cancel')
%    case 'No', return
%    case 'Cancel', return
%end

sThisTrial = FV.sLoadedTrial;
FV = rmfield(FV, {'tSpikes' 'tData'}); % remove data fields
FV.tSpikes = struct([]);
FV.tData = struct([]);
FV.vXlim = []; % remove the x-limits since they might not apply to other files

% Iterate over all files in directory
SpikyWaitbar(0, length(csPaths));
for f = 1:length(csPaths)
    SpikyWaitbar(f, length(csPaths));

    % set sLoadedTrial
    sNewFile = CheckFilename(csPaths{f});
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
    save(CheckFilename([FV.sLoadedTrial(1:end-4) '.spb']), 'FV', '-MAT', '-V6')
end
sp_disp(sprintf('Done distributing current settings to %d files in current directory', length(csPaths)))
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CreateMergeFile(varargin)
% Create a merge file containing spikes and EMGs of all files in this
% directory. Function asks which data to merge, and whether to merge with
% existing merge-data already on disks.
global g_hSpike_Viewer g_vSelCh
[FV,hWin] = GetStruct;

if ~CheckDataLoaded(), return, end

% Ask user what data to merge
hFig = figure;
cFields = {'Spikes' 'Events' 'Imported channels' 'Filtered channels'};
nCL = length(cFields)*25+5;
vHW = get(g_hSpike_Viewer, 'position');
set(hFig, 'position', [vHW(1:2) 200 nCL+25], 'menu', 'none', 'Name', 'Merge Data', 'NumberTitle', 'off')
centerfig(hFig, g_hSpike_Viewer)
for i = 1:length(cFields)
    uicontrol(hFig, 'Style', 'checkbox', 'Position', [10 nCL 150 20], 'String', [' ' cFields{i}], ...
        'HorizontalAlignment', 'left', 'backgroundcolor', [.8 .8 .8], 'value', 1);
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

% Ask whether to save over existing file if it exists
sPath = CheckFilename([FV.sDirectory '\MergeFile.spb']);
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
tEMG = struct([]);      % temporary variable that holds filtered EMG data
tImported = struct([]); % temporary variable that holds imported data
FV_merge.sDirectory = FV.sDirectory;
sDiffSampleRateAction = 'none';

if strcmpi(sAppendReplace, 'Append')
    % Get existing merge data
    % Note that we don't append to existing spike data, since spikes may
    % have been dejittered (in which case matrices wont merge readily).
    % Thus, whenever spikes are selected to be merged, they will always
    % replace ALL spikes in the merge structure, regardless.
    FV_old = load(CheckFilename([FV_merge.sDirectory '\MergeFile.spb']), '-MAT');
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
    warndlg('No additional files have been loaded into your workspace (see File->Files). You first need to load a directory, a directory tree or a custom list of files into your workspace before distributing settings and merging results.', 'Spiky!');
    return
end
bResult = SpikyWaitbar(0, length(csPaths));
for f = 1:length(csPaths)
    bWaitResult = SpikyWaitbar(f, length(csPaths));
    if ~bWaitResult % waitbar was closed
        StopMergeMode();
        SpikyWaitbar(length(csPaths), length(csPaths));
        uiwait(warndlg('Merge mode has been cancelled.'))
        return
    end
    bResult = LoadTrial(csPaths{f}); % load trial data and settings
    drawnow
    [FV,hWin] = GetStruct;
    
    % Get spike mergedata
    if bSpikes
        DetectSpikes(); % detect spikes on thresholded channels
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
        % get spikes from all thresholded channels
        csChannels = fieldnames(FV.tSpikes);
        % Collect spiketime and associated information
        for nCh = 1:length(csChannels)
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
                            sAns = questdlg(sMsg, 'Spiky!', 'Ignore Spikes', 'Resample Spikes', 'Cancel Merging', 'Resample Spikes');
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
            if nFs1 ~= nFs2 % resample or ignore if samplerates are different in this file
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
            
            % Add fields that indicate start and end of file
            if ~isfield(FV_merge.tData, 'FileStart'), FV_merge.tData(1).FileStart = {}; end
            if ~isfield(FV_merge.tData, 'FileEnd'), FV_merge.tData(1).FileEnd = {}; end
            
            FV_merge.tData(1).FileStart(end+1).Timestamp = FV.tData.([FV.csDigitalChannels{c} '_TimeBegin']); % file start marker
            FV_merge.tData(1).FileStart(end).File = FV.sLoadedTrial;
            FV_merge.tData(1).FileEnd(end+1).Timestamp = FV.tData.([FV.csDigitalChannels{c} '_TimeEnd']); % file end marker
            FV_merge.tData(1).FileEnd(end).File = FV.sLoadedTrial;
        end
    end % end of events merging
    
    % Merge imported fields
    if bImported
        csFields = fieldnames(FV.tData);
        csChannels = {};
        for nCh = 1:length(csFields) % find imported fields
            if ~isempty(strfind(csFields{nCh}, '_Imported')), csChannels{end+1} = csFields{nCh}(1:end-9); end
        end
        csChannels = unique(csChannels);
        for nCh = 1:length(csChannels) % iterate over imported channels
            vCont = FV.tData.(csChannels{nCh});
            [vCont, FV] = AdjustChannelGain(FV, vCont, csChannels{nCh}); % Adjust gain (mV)
            % time vector
            nBegin = FV.tData.([csChannels{nCh} '_TimeBegin']); % sampling start, sec
            nEnd = FV.tData.([csChannels{nCh} '_TimeEnd']); % sampling end, sec
            nFs = FV.tData.([csChannels{nCh} '_KHz']) * 1000; % sampling frequency Hz
            vTime = (nBegin+1/nFs):(1/nFs):(nBegin+length(vCont)/nFs); % absolute time, sec
            % move data into temporary variable
            if ~isfield(tImported, csChannels{nCh}), nIndx = 1;
            else nIndx = length(tImported.(csChannels{nCh})) + 1; end
            tImported(1).(csChannels{nCh})(nIndx).cont = vCont;
            tImported.(csChannels{nCh})(nIndx).timebegin = nBegin; % store begin time
            tImported.(csChannels{nCh})(nIndx).timeend = nEnd; % store end time
            tImported.(csChannels{nCh})(nIndx).time = vTime; % store continuous time vector, sec
            tImported.(csChannels{nCh})(nIndx).fs = nFs; % sampling rate
        end
    end
    
    % Merge filtered signals (EMG, LFP, ECoG etc)
    % Note: Any signal that has been Filtered through the GUI will be exported to
    % the Merge File (store in temporary var first and join with FV_merge after outer loop)
    if bFiltered
        csChannels = FV.csEMGChannels;
        for nCh = 1:length(csChannels)
            % check first that channel exists
            if ~isfield(FV.tData, csChannels{nCh}) continue, end
            % raw signal
            sCh = csChannels{nCh};
            vCont = FV.tData.(csChannels{nCh});
            [vCont, FV] = AdjustChannelGain(FV, vCont, sCh); % Adjust gain (mV)
            % time vector
            nBegin = FV.tData.([csChannels{nCh} '_TimeBegin']); % sampling start, sec
            nEnd = FV.tData.([csChannels{nCh} '_TimeEnd']); % sampling end, sec
            nFs = FV.tData.([csChannels{nCh} '_KHz']) * 1000; % sampling frequency Hz
            vTime = (nBegin+1/nFs):(1/nFs):(nBegin+length(vCont)/nFs); % absolute time, sec
            % filter and rectify signal
            [vCont, vTime, nNewFs] = FilterChannel(vCont, vTime, nFs, FV.nEMGLoPass, FV.nEMGHiPass, FV.nEMGRectify, 'decimate');
            % store signal in temporary variable
            if ~isfield(tEMG, csChannels{nCh}), nIndx = 1;
            else nIndx = length(tEMG.(csChannels{nCh})) + 1; end
            tEMG(1).(csChannels{nCh})(nIndx).cont = vCont;
            tEMG.(csChannels{nCh})(nIndx).timebegin = nBegin; % store begin time
            tEMG.(csChannels{nCh})(nIndx).timeend = nEnd; % store end time
            tEMG.(csChannels{nCh})(nIndx).time = vTime; % store continuous time vector, sec
            tEMG.(csChannels{nCh})(nIndx).fs = nNewFs; % new sampling rate
        end
    end
    
    % Keep ALL channel descriptions (regardless of whether they were included in the merge file or not)
    % As first step, check for and ignore duplicates
    for i = 1:length(FV.tChannelDescriptions)
        sCh = FV.tChannelDescriptions(i).sChannel;
        if isempty(FV_merge.tChannelDescriptions)
            FV_merge.tChannelDescriptions = FV.tChannelDescriptions;
        else
            j = find(strcmp({FV_merge.tChannelDescriptions.sChannel}, sCh));
            if isempty(j)
                FV_merge.tChannelDescriptions(end+1).sChannel = FV.tChannelDescriptions(i).sChannel;
                FV_merge.tChannelDescriptions(end).sDescription = FV.tChannelDescriptions(i).sDescription;
            end
        end
    end
end

% Move imported data into FV_merge after creating one continuous vector
% containing the concatenated signal from all files. Missing data segments
% are filled with NaNs
if bImported
    cFields = fieldnames(tImported);
    for nF = 1:length(cFields)
        % initialize continuous vector (min begintime to max endtime)
        vTimeBegin = [tImported.(cFields{nF}).timebegin]; % file endtimes, sec
        vTimeEnd = [tImported.(cFields{nF}).timeend]; % file endtimes, sec
        vFs = [tImported.(cFields{nF}).fs];
        if length(unique(vFs)) > 1
            uiwait(warndlg('Sorry, cannot export all imported channels. Sampling rate must be same across all files. Use of different sampling rates for the same channel is currently not supported. Will continue to try to export remaining imported channels.'))
            continue
        end
        nD = 1/vFs(1); % duration of one sampling points
        vContAll = [];
        nAbsBeginTime = round(min(vTimeBegin) / nD - 1); % begin time of first trial, samples
        %nAbsBeginTime = round(min(vTimeBegin) / nD); % begin time of first trial, samples
        vIndxAll = [];
        for i = 1:length(vTimeEnd)
            nBegin = tImported.(cFields{nF})(i).timebegin; % sec
            %nBeginTime = nBegin; % NEW
            nBegin = round(nBegin / nD); % samples
            %nRelBeginToFirst = round((nBeginTime-min(vTimeBegin)) / nD) + 1; % samples NEW
            vCont = tImported.(cFields{nF})(i).cont;
            vIndx = (nBegin:[nBegin+length(vCont)-1])-nAbsBeginTime;
            
            %vIndx = nRelBeginToFirst:( nRelBeginToFirst+length(vCont) - 1 );
            vIndxAll = [vIndxAll vIndx];
            vContAll(vIndx) = vCont; % samples
            vIndxPrev = vIndx;
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

% Move FILTERED data into FV_merge after creating one continuous vector
% containing the concatenated signal from all files. Missing data segments
% are filled with NaNs
if bFiltered
    cFields = fieldnames(tEMG);
    for nF = 1:length(cFields)
        % initialize continuous vector (min begintime to max endtime)
        vTimeBegin = [tEMG.(cFields{nF}).timebegin]; % file endtimes, sec
        vTimeEnd = [tEMG.(cFields{nF}).timeend]; % file endtimes, sec
        vFs = [tEMG.(cFields{nF}).fs];
        if length(unique(vFs)) > 1
            uiwait(warndlg('Sorry, cannot export all imported channels. Sampling rate must be same across all files. Use of different sampling rates for the same channel is currently not supported. Will continue to try to export remaining imported channels.'))
            continue
        end
        nD = 1/vFs(1); % duration of one sampling points
        vContAll = [];
        nAbsBeginTime = round(min(vTimeBegin) / nD - 1);
        vIndxAll = [];
        for i = 1:length(vTimeEnd)
            nBegin = tEMG.(cFields{nF})(i).timebegin; % sec
            nBegin = round(nBegin / nD); % samples
            vCont = tEMG.(cFields{nF})(i).cont;
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

% Set FV as FV_merge
if strcmpi(sAppendReplace, 'Replace')
    % In Replace mode, reset FV to defaults
    FV = SetFVDefaults();
elseif strcmpi(sAppendReplace, 'Append')
    % In Append mode, load and use as basis existing mergefile
    FV_old = load(CheckFilename([FV_merge.sDirectory '\MergeFile.spb']), '-MAT');
    FV = FV_old.FV;
end
if bSpikes, FV.tSpikes = FV_merge.tSpikes; end
FV.tData = FV_merge.tData;
FV.sDirectory = FV_merge.sDirectory;
FV.csDisplayChannels = FV_merge.csDisplayChannels;
FV.csDigitalChannels = FV_merge.csDigitalChannels;
FV.tChannelDescriptions = FV_merge.tChannelDescriptions;
FV.csEMGChannels = {}; % empty since we filtered above
FV.bPlotRasters = 1;
FV.sLoadedTrial = CheckFilename([FV_merge.sDirectory '\MergeFile.spb']);
SetStruct(FV)
save(FV.sLoadedTrial, 'FV', '-MAT', '-V6') % Save merged results to disk
if length(FV.sLoadedTrial) > 70
    sLoadedTrial = ['...' FV.sLoadedTrial(end-69:end)];
else sLoadedTrial = FV.sLoadedTrial; end
set(g_hSpike_Viewer, 'Name', sprintf('Spiky! - %s', sLoadedTrial))
ViewTrialData();
ChannelStatistics(); % show channel statistics

% Check merge file for potential problems
vStart = FV.tData.DAQ_Start_Up; % start time duplicates
if length(unique(vStart)) < length(vStart)
    uiwait(warndlg(sprintf('There is a potential problem with the merged files:\nThere are duplicates of file Start times')))
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StartMergeMode(varargin)
% Initialize workspace for being in MERGE MODE
% print message that we are in a MERGE MODE
[FV,hWin] = GetStruct;
global g_hSpike_Viewer g_bMergeMode
figure(g_hSpike_Viewer)
g_bMergeMode = 1;
% Update menu items
set(findobj(g_hSpike_Viewer, 'Label', 'E&xit Merge Mode'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', 'Set &Threshold...'), 'Enable', 'off')
set(findobj(g_hSpike_Viewer, 'Label', '&Auto-Detect Thresholds...'), 'Enable', 'off')
set(findobj(g_hSpike_Viewer, 'Label', '&Run Spike Detection'), 'Enable', 'off')
set(findobj(g_hSpike_Viewer, 'Label', 'Add &Multitrode...'), 'Enable', 'off')
set(findobj(g_hSpike_Viewer, 'Label', 'Remove Multitrode...'), 'Enable', 'off')
set(findobj(g_hSpike_Viewer, 'Label', 'Set Max &Jitter...'), 'Enable', 'off')
set(findobj(g_hSpike_Viewer, 'Label', 'Set Pre-Trigger...'), 'Enable', 'off')
set(findobj(g_hSpike_Viewer, 'Label', 'Set Post-Trigger...'), 'Enable', 'off')
%set(findobj(g_hSpike_Viewer, 'Label', 'Set Spike Deadtime...'), 'Enable', 'off')
set(findobj(g_hSpike_Viewer, 'Label', '&Create Merge File'), 'Enable', 'off')
set(findobj(g_hSpike_Viewer, 'Label', '&Distribute Settings'), 'Enable', 'off')
set(findobj(g_hSpike_Viewer, 'Label', 'Show P&ower Spectral Densities'), 'Enable', 'off')
set(findobj(g_hSpike_Viewer, 'Label', 'Digitize Channel'), 'Enable', 'off')
set(findobj(g_hSpike_Viewer, 'Label', 'Digitize Manually'), 'Enable', 'off')
hRaster = findobj(g_hSpike_Viewer, 'Label', 'Plot &Rasters');
set(hRaster, 'Enable', 'off', 'Checked', 'on')
FV.bPlotRasters = 1;
SetStruct(FV)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StopMergeMode(varargin)
% Clean up and exit from merge mode
global g_hSpike_Viewer g_bMergeMode
g_bMergeMode = 0;
% Update menu items
set(findobj(g_hSpike_Viewer, 'Label', 'Set &Threshold...'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', '&Auto-Detect Thresholds...'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', '&Run Spike Detection'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', 'Define &Multitrode...'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', 'Set Max &Jitter...'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', 'Set Pre-Trigger...'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', 'Set Post-Trigger...'), 'Enable', 'on')
%set(findobj(g_hSpike_Viewer, 'Label', 'Set Spike Deadtime...'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', '&Create Merge File'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', '&Distribute Settings'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', 'Show P&ower Spectral Densities'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', 'Digitize Channel'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', 'Digitize Manually'), 'Enable', 'on')
set(findobj(g_hSpike_Viewer, 'Label', 'Plot &Rasters'), 'Enable', 'on')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadMergeFile(varargin)
% Load a file containing data of already merged files
[FV,hWin] = GetStruct;
% Check if file path is given in function input
if ~isempty(varargin{1}) && ischar(varargin{1}) && length(varargin) == 2
    sPath = varargin{1};
    sFile = varargin{2};
else
    % Select merge file
    [sFile, sPath] = uigetfile( {'*.spb','Spiky! files (*.spb)'}, 'Select merged file');
    if sFile == 0, return, end
end
FV.sDirectory = sPath;
% Load merge file
StartMergeMode();
sFilePath = CheckFilename([sPath sFile]);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FigureExport(varargin)
% Export the Spiky! window to a graphics file
global g_hSpike_Viewer % handle to Spiky! window
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
    'Export Spiky! window as');
if all(sFile == 0), return, end % user cancelled save dialog
saveas(g_hSpike_Viewer, [sPath sFile])
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InvertColors(varargin)
% Toggle the InvertHardcopy parameter of the Spiky! window. This feature
% decides whether the background color of the window should be included in
% exports.
global g_hSpike_Viewer % handle to Spiky! window
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChannelStatistics(varargin)
% Show basic channel statistics, such as number of detected spikes
% Format:
% Channel_name  Fs  #spikes
[FV,hWin] = GetStruct;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomAmplitude(varargin)
% Zoom on y-axis in selected axis
if ~CheckDataLoaded, return, end
[FV,hWin] = GetStruct;
global g_hSpike_Viewer % handle to Spiky! window
sPointer = get(g_hSpike_Viewer,'Pointer');
set(g_hSpike_Viewer,'Pointer','crosshair')
zoom off
waitfor(g_hSpike_Viewer, 'CurrentPoint') % wait for user to click in figure
hAx = gca;
set(hAx, 'color', [.2 .1 .1])
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
    sCh = get(hAx,'Tag');
    FV.tYlim.(sCh) = get(hAx, 'ylim');
    set(hAx, 'color', [.1 .1 .1])
end
SetStruct(FV)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotRasters(varargin)
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
% update menu item
global g_hSpike_Viewer % handle to Spiky! window
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotInstSpikeRate(varargin)
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
% update menu item
global g_hSpike_Viewer % handle to Spiky! window
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExitSpiky(varargin)
% Quit and exit Spiky!
[FV,hWin] = GetStruct;
global g_hSpike_Viewer
% Check whether there are changes to be saved
if ~isempty(strfind(get(g_hSpike_Viewer, 'Name'), '*'))
    switch questdlg('Save changes to current file?', 'Spiky!', 'Yes', 'No', 'Cancel', 'Yes')
        case 'Yes', SaveResults;
        case 'Cancel', return
    end
end
% Close Spiky! windows
delete(g_hSpike_Viewer)
CloseSpikyWindows
% Clear global variables
clear g_hSpike_Viewer g_bMergeMode
clear OpenSettings % clear persistent variables in OpenSettings()
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function QuickSort(varargin)
% QuickSort! is a function for automatically detectingspikes, dejittering, removing
% outliers and sorting all currently displayed channels. Its a tool for quickly
% giving an impression of the data. Some steps are skipped in Merge Mode. It assumes
% thresholds has already been define.
[FV, hWin] = GetStruct;
global g_bMergeMode
if ~CheckDataLoaded, return, end
if ~isfield(FV, 'tData'), return, end % check data has been loaded
%try
if ~g_bMergeMode % skip spike detection if we're in Merge Mode
    DetectSpikes()
    if isempty(FV.tSpikeThresholdsNeg) % check thresholds have been set
        uiwait(warndlg('Threshold have not been defined. QuickSort! was aborted.', 'Spiky!'))
        return
    end
end
DejitterSpikes()
RemoveOutlierSpikes()
SortAllChannels()
%catch uiwait(warndlg(sprintf('QuickSort! failed with error messsage:\n%s', lasterr), 'Spiky!')); end % if the try failed, it usually means user aborted
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowSpectralDensity(varargin)
% Compute and display the power spectral densities of displayed analoge
% channels. This tool is useful for detecting interesting frequencies in
% the signals, such as line-noise or LFP oscillations.
[FV,hWin] = GetStruct;
if ~CheckDataLoaded, return, end

hFig = figure;
set(hFig, 'color', [.2 .2 .2], 'name', 'Spiky! Spectral Densities', 'NumberTitle', 'off')
hAx = axes('position', [.1 .1 .74 .8]);
hold on
sCh = FV.csDisplayChannels;
for nCh = 1:length(sCh)
    vData = eval(['FV.tData.' FV.csDisplayChannels{nCh}]);
    if isempty(vData) continue; end
    % Fill in NaN's with nearest non-NaN neighbour
    vData(isnan(vData)) = interp1(find(~isnan(vData)), vData(~isnan(vData)), find(isnan(vData)), 'linear');
    nFs = eval(['FV.tData.' FV.csDisplayChannels{nCh} '_KHz']); % Fs in KHz
    [vPxx, vF] = pwelch(vData, kaiser(length(vData),4), 0, 2^16, nFs*1000);
    vF_i = logspace(log10(1), log10(max(vF)), 250);
    vPxx = interp1(vF, vPxx, vF_i, 'pchip');
    mCol = FV.mColors(nCh,:);
    hLines = plot(vF_i, vPxx, 'color', mCol);
    % Toggle view check box
    set(hLines, 'tag', char(nCh));
    hCheckbox = uicontrol(hFig, 'Style', 'checkbox', 'units', 'normalized', ...
        'Position', [.85 .85-([nCh-1]*.06) .15 .05], 'String', sCh(nCh), ...
        'HorizontalAlignment', 'left', 'backgroundcolor', [.2 .2 .2], ...
        'value', 1, 'foregroundColor', mCol);
    set(hCheckbox, 'Tag', char(nCh), 'callback', @ToggleWaveforms);
end
xlabel('Frequency (Hz)');
ylabel('Power (arb.)')
hHeader = header('Power Spectral Densities', 12);
set(hHeader, 'color', 'w', 'interpreter', 'none')
set(hAx, 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'xscale', 'log', 'yscale', 'li', 'xtick', [.1 1 10 100 1000 10000], 'xticklabel', [.1 1 10 100 1000 10000])
axes(hAx); axis tight
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sTxt] = ShowCurrentLine(obj, event_obj)
sTxt = {['']};

sTypeOfCurrObj = get(gco, 'Type');
if ~strcmp(lower(sTypeOfCurrObj), 'line'), return, end

% Reset color on unclicked lines
hToggleHandles = findobj(gcf, 'Tag', get(gco,'tag'), 'Type', 'line');
cCol = get(hToggleHandles(2:end),'color');
mCol = reshape([cCol{:}],3,length(hToggleHandles)-1)';
set(hToggleHandles(2:end), 'color', median(mCol,1), 'linewidth', 1, 'marker', 'none');

% Reset marker on all line objects
hHandles = findobj(gcf, 'Type', 'line');
set(hHandles, 'marker', 'none')

set(gco, 'color', median(mCol,1), 'linewidth', 2, 'marker', 's', 'markeredgecolor', 'w', 'markersize', 4)

% If an UIHoldFigure figure is open, close it/them
close(findobj('Tag', 'UIHoldFigure'));

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IdentifySpike(varargin)
[FV,hWin] = GetStruct;
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
text(vXY(1,1), vXY(1,2), sTxt, 'color', vCol, 'fontsize', 7)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImportData(varargin)
[FV,hWin] = GetStruct;
persistent p_cDataFields p_sGainField p_sFsField p_sONField p_sDefFile
sPath = [FV.sLoadedTrial(1:end-4) '.spb'];

% Filename without extension
vI = unique([strfind(sPath, '/') strfind(sPath, '\')]);
if isempty(vI)
    sFilename = sPath(1:end-4);
else
    sFilename = sPath(vI(end)+1:end-4);
end
    
[sFile, sPath] = uigetfile( {[sFilename '*.mat'], ['Matching MAT files (' sFilename '*.mat)'];
    [sFilename '.*'],  ['Matching files (' sFilename '.*)']; ...
    '*.daq',  'Data Acquisition Toolbox files (*.daq)'; ...
    '*.txt',  'Text files (*.txt)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Select data file', p_sDefFile);

if sFile == 0, return, end
p_sDefFile = sFile;

% Set current directory to sPath
cd(sPath)

switch sFile(end-3:end)
    case '.txt'
        uiwait(warndlg('Import of text files not yet functional in this version of Spiky!')); return
    case '.daq'
        % Import data from Data Acquisition Toolbox files
        [mData, vTime, vAbsTime, tTvents, tDAQInfo] = daqread([sPath sFile]);

        cFields = {tDAQInfo.ObjInfo.Channel(:).HwChannel};
        vFields = [tDAQInfo.ObjInfo.Channel(:).HwChannel];

        % Select data vector(s)
        [sDataField, nManualValue] = SelectImportVariable(cFields, 'Set data vector(s)', {'Set data vector(s)'}, 0, NaN, [], 1);
        vData = mData(:, vFields == sDataField);

        % Get sampling rate (kHz)
        nFs = tDAQInfo.ObjInfo.SampleRate / 1000;

        % Select gain
        [sGainField, nManualValue] = SelectImportVariable(cFields, 'Select gain field', {'Select gain field'}, 1, 10000, [], 0);
        if ~isnan(nManualValue), nGain = nManualValue;
        else nGain = tImportedData.(sGainField); end

        % Reset FV structure
        FV = SetFVDefaults();

        % Assign new values to FV
        FV.sDirectory = sPath;
        FV.sLoadedTrial = [sPath sFile];
        FV.tGain(1).sDataField = nGain;
        FV.tData = struct([]);
        FV.tData(1).DAQImport = vData';
        FV.tData(1).('DAQImport_KHz') = nFs;
        FV.tData(1).('DAQImport_KHz_Orig') = nFs;
        FV.tData(1).('DAQImport_TimeBegin') = 0;
        FV.tData(1).('DAQImport_TimeEnd') = 1;

    case '.mat'
        % Import vectorized data from .mat files
        tImportedData = load([sPath sFile], '-MAT');
        cFields = fieldnames(tImportedData);
        
        % Select data vector(s)
        [cDataFields, nManualValue] = SelectImportVariable(cFields, 'Set data vector(s)', {'Set data vector(s)'}, 0, NaN, p_cDataFields, 1);
        p_cDataFields = cDataFields;
        if isempty(cDataFields), return, end % user closed window
        
        % Check data vector is a vector (minimum length 2)
        for c = 1:length(cDataFields)
            if any(size(tImportedData.(cDataFields{c})) == 1) % vector data
                if length(tImportedData.(cDataFields{c})) < 2
                    uiwait(warndlg(sprintf('Import failed. Input field %s has a length of less than 2 number. Must have a length greater than this.', cDataFields{c})));
                    return
                end
            end
        end
        
        % Select sampling rate (kHz)
        [cFsField, nManualValue] = SelectImportVariable(cFields, 'Set sampling rate (Hz)', {'Set sampling rate (Hz)'}, 1, 40000, p_sFsField, 0);
        if isempty(cFsField), return, end % user closed window
        p_sFsField = cFsField{1};
        if ~isnan(nManualValue), nFs = double(nManualValue) / 1000;
        else nFs = double(tImportedData.(p_sFsField)) / 1000; end
        
        % Select gain
        [cGainField, nManualValue] = SelectImportVariable(cFields, 'Set gain', {'Set gain'}, 1, 10000, p_sGainField, 0);
        if isempty(cGainField), return, end % user closed window
        p_sGainField = cGainField{1};
        if ~isnan(nManualValue), nGain = double(nManualValue);
        else nGain = double(tImportedData.(p_sGainField)); end
        
        % Ask whether data should be appended to existing channels.
        % Note: Append is not supported with matrix import (i.e. importing multiple channels at once)
        bAppend = 1;
        for c = 1:length(cDataFields)
            if numel(tImportedData.(cDataFields{c})) > sum(size(tImportedData.(cDataFields{c})))
                bAppend = 0;
            end
        end
        
        if ~isempty(FV.sLoadedTrial) && bAppend
            sAns = questdlg('Open in new workspace or append to existing?', 'Import data', 'New Workspace', 'Append', 'Cancel', 'Append');
            switch sAns
                case 'Append'
                    % If data is to be appended, get onset time
                    bAppend = 1;
                    cFields = fieldnames(FV.tData);
                    cDescriptions = {};
                    vIndx = [];
                    for i = 1:length(cFields)
                        if ~isempty([strfind(cFields{i}, '_Up') strfind(cFields{i}, '_Down')])
                            if ~isempty(FV.tData.(cFields{i}))
                                vIndx(end+1) = i;
                                % Get channel description
                                nStartIndx = [strfind(cFields{i}, '_Up') strfind(cFields{i}, '_Down')];
                                sStrMatch = cFields{i}(1:nStartIndx-1);
                                nMatchIndx = find(strcmpi({FV.tChannelDescriptions.sChannel}, sStrMatch));
                                cDescriptions{i} = FV.tChannelDescriptions(nMatchIndx).sDescription;
                            end
                        end
                    end
                    
                    % Get alignment event
                    [cONField, nManualValue, nButton] = SelectImportVariable([cFields(vIndx); 'Append to existing channels'], ...
                        'Align to event', {'Align Start', 'Align End'}, 1, 0, p_sONField, 0, ...
                        [cDescriptions(vIndx)'; ' ']);
                    if isempty(cONField), return, end % user closed window
                    p_sONField = cONField{1};
                    switch nButton
                        case 1 % align start of signal to trigger
                            sAlign = 'start';
                        case 2 % align end of signal to trigger
                            sAlign = 'end';
                    end
                    
                case 'New Workspace', bAppend = 0;
                case 'Cancel', return;
            end
        else bAppend = 0; end

         % Dont append and reset FV structure
        if ~bAppend
            FV = SetFVDefaults();
            FV.sDirectory = sPath;
            FV.tData = struct([]);
            FV.sLoadedTrial = [sPath sFile];
            %FV.tData(1).([sDataField '_TimeBegin']) = 0;
        end
        
        % Iterate over all fields to be imported
        for c = 1:length(p_cDataFields)
            sDataField = p_cDataFields{c};
            
            % Get data and create inner loop for data matrices
            mData = tImportedData.(sDataField);
            
            % Transpose data vector, if needed, so that column(s) equal channels.
            vSize = size(mData);
            if vSize(2) > vSize(1)
                mData = mData';
            end
            
            for di = 1:size(mData, 2)
                
                % Temporal alignment offset
                if ~isnan(nManualValue)
                    nOnsetTime = nManualValue;
                else
                    if ~strcmpi(cONField, 'Append to existing channels')
                        switch lower(sAlign)
                            case 'start' % align start of signal to trigger
                                nOnsetTime = FV.tData.(p_sONField)(1);
                            case 'end' % align end of signal to trigger
                                nOffsetTime = FV.tData.(p_sONField)(1); % end of signal, sec
                                %nDur = length(tImportedData.(sDataField)(:)) /  (nFs*1000);
                                nDur = length(mData(:,di)) /  (nFs*1000);
                                nOnsetTime = nOffsetTime - nDur; % start of signal, sec
                        end
                    end
                end
                
                if strcmpi(cONField, 'Append to existing channels')
                    % Append new data to an existing field
                    if isfield(FV.tData, sDataField)
                        %FV.tData.(sDataField) = [FV.tData(1).(sDataField) double(tImportedData.(sDataField)(:)')];
                        FV.tData.(sDataField) = [FV.tData(1).(sDataField) double(mData(:, di))];
                        %FV.tData(1).([sDataField '_TimeEnd']) = FV.tData(1).([sDataField '_TimeEnd']) + ...
                        %    (length(tImportedData.(sDataField)(:)') / (nFs*1000)); % sec
                        FV.tData(1).([sDataField '_TimeEnd']) = FV.tData(1).([sDataField '_TimeEnd']) + ...
                            (length(mData(:, di)) / (nFs*1000)); % sec
                    else
                        uiwait(warndlg(sprintf('The channel %s does not exist and therefore cannot be appended with new data.', sDataField)));
                        continue
                    end
                else
                    % Create field for data
                    sChName = sprintf('%s%d', sDataField, di);
                    FV.tGain(1).(sChName) = nGain;
                    FV.tData(1).([sChName '_KHz']) = nFs;
                    FV.tData(1).([sChName '_KHz_Orig']) = nFs;
                    %FV.tData(1).([sChName '_TimeEnd']) = nOnsetTime + (length(tImportedData.(sDataField)(:)') / (nFs*1000)); % sec
                    FV.tData(1).([sChName '_TimeEnd']) = nOnsetTime + (length(mData(:, di)) / (nFs*1000)); % sec
                    %FV.tData(1).(sChName) = double(tImportedData.(sDataField)(:)');
                    FV.tData(1).(sChName) = double(mData(:, di));
                    FV.tData(1).([sChName '_Imported']) = 1; % Mark field as imported
                    FV.csDisplayChannels{end+1} = sChName;
                    FV.tData(1).([sChName '_TimeBegin']) = nOnsetTime; % sec
                end
            end % end of channels loop
            
        end
    otherwise
        uiwait(warndlg('No valid file selected')); return
end
SetStruct(FV)
ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
if bManualOpt, cFields{end+1} = 'Set Manually'; end

% Modify fieldnames to include description strings
if ~isempty(varargin)
    cDescr = varargin{1};
    cFieldsMod = {};
    for i = 1:length(cFields)-1
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

if bMultiSel nMax = 10;
else nMax = 1; end

hFig = figure;
set(hFig, 'Menu', 'none', 'Name', sTitle, 'NumberTitle', 'off', 'Position', [200 600 300 length(cFields)*20+25]);
hList = uicontrol('Position', [0 26 300 length(cFields)*20], 'Style', 'listbox', ...
    'String', cFieldsMod, 'FontSize', 10, 'Tag', 'ImportListBox', 'Max', nMax);
if ~isempty(vDefField), set(hList, 'Value', vDefField); end

for b = 1:length(csButton)
    nW = 300 / length(csButton);
    nX = nW * (b - 1);
    uicontrol('Position', [nX 0 nW 25], 'Style', 'pushbutton', 'String', csButton{b}, 'FontSize', 10, ...
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
    else nManualValue = NaN; end
end

clear global nSelectedIndx
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotSpikeTimeDistributions(varargin)
[FV,hWin] = GetStruct;
% Plot rasters of sorted units
if ~CheckDataLoaded, return, end
[FV,hWin] = GetStruct;
hFig = figure; % create figure
set(hFig, 'color', [.2 .2 .2], 'Name', 'Spiky! Spiketime Distributions', 'NumberTitle', 'off', 'menubar', 'none')
hAx = axes('position', [0.1 .075 .88 .9]);
set(hAx, 'color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'fontsize', 8)
hold on
% Iterate over channels
csChannels = fieldnames(FV.tSpikes);
nRow = 0;
sUnits = {}; nU_max = 0;
for nCh = 1:length(csChannels)
    % Iterate over units
    tSpikes = FV.tSpikes.(csChannels{nCh});
    if ~isfield(tSpikes, 'hierarchy'), continue, end
    vUnits = unique(tSpikes.hierarchy.assigns);
    nFs = tSpikes.Fs(1);
    mBoxes = []; mCols = [];
    for nU = 1:length(vUnits)
        nRow = nRow + 1;
        vIndx = find(tSpikes.hierarchy.assigns == vUnits(nU));
        vSpiketimes = tSpikes.spiketimes(vIndx); % samples
        vSpiketimes = vSpiketimes ./ nFs; % sec
        mBoxes = [mBoxes; vSpiketimes repmat(nU, length(vSpiketimes), 1)];
        if vUnits(nU) == 0, mCols(end+1, :) = [.4 .4 .4];
        else mCols(end+1, :) = FV.mColors(nU,:); end
        %plot(vSpiketimes, repmat(nRow,length(vSpiketimes),1), '.', 'color', mCols(end,:), 'linewidth', .5);
        sUnits{end+1} = [csChannels{nCh} '_' num2str(vUnits(nU))];
    end
    boxplot(mBoxes(:,1), mBoxes(:,2), 'orientation', 'horizontal', 'colors', mCols, 'positions', nU_max+(1:nU))
    nU_max = nU_max + mBoxes(end,2);
end
if nRow == 0 % if no units have been sorted
    close(hFig)
    uiwait(warndlg('No units have been sorted.'))
    return
end
xlabel('Time (sec)')
set(hAx, 'ylim', [0 nRow+1], 'ytick', 1:nRow, 'yticklabel', sUnits, 'ygrid', 'off', 'xgrid', 'on')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotCrossCorrelationsDiscrete(varargin)
[FV,hWin] = GetStruct;
[sCh, bResult] = SelectChannelNumber(fieldnames(FV.tSpikes)');

hFig = figure;
set(hFig, 'color', [.2 .2 .2], 'name', 'Spiky! Discrete Cross Correlations', 'NumberTitle', 'off');
vXLim = [-.05 .05];

% Iterate over units (plot auto-corrs and cross-corrs)
vUnits = unique(FV.tSpikes.(sCh).hierarchy.assigns);
nFs = FV.tData.([sCh '_KHz']) * 1000;
for u1 = 1:length(vUnits)
    vIndx1 = FV.tSpikes.(sCh).hierarchy.assigns == vUnits(u1);
    vSpiketimes1 = FV.tSpikes.(sCh).spiketimes(vIndx1) / nFs; % sec
    % limit number of spikes to 10000
    nMaxLen = 10000;
    if length(vSpiketimes1) > nMaxLen
        vRand = randperm(length(vSpiketimes1));
        vSpiketimes1 = vSpiketimes1(vRand(1:nMaxLen));
    end
    for u2 = 1:length(vUnits)
        if u1 > u2, continue, end
        vIndx2 = FV.tSpikes.(sCh).hierarchy.assigns == vUnits(u2);
        vSpiketimes2 = FV.tSpikes.(sCh).spiketimes(vIndx2) / nFs; % sec
        % limit number of spikes to 10000
        if length(vSpiketimes2) > nMaxLen
            vRand = randperm(length(vSpiketimes2));
            vSpiketimes2 = vSpiketimes2(vRand(1:nMaxLen));
        end
        nX = .08; nY = .1; nW = .88; nH = .85; nNN = length(vUnits);
        hAx = axes('position', [nX+(u2-1)*(nW/nNN) (1-nY)-u1*(nH/nNN) nW/nNN nH/nNN], 'Color', [.1 .1 .1]);

        % Compute cross-correlation
        try
            [vC, vT] = ccorr(vSpiketimes1', vSpiketimes2', vXLim, 'n', [], 1, 1, 200);
        catch
            continue
        end
        vT = vT*1000; % ms

        if all(vUnits([u1 u2]) == 0), vCol = [.5 .5 .5]; % outliers
        elseif u1 == u2, vCol = FV.mColors(u1,:);
        else vCol = [.7 .7 .7]; end
        plot(vT, vC, 'color', vCol)
        hold on
        plot([0 0], [0 max(vC)*2], ':', 'color', [.6 .6 .6])
        set(hAx, 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'xlim', vXLim*1000, ...
            'fontsize', 7, 'ylim', [0 max(vC)+.0001] )
        if u1 > u2, set(hAx, 'xtick', []); end
        if u1 == 1
            hTit = title(sprintf(' Unit %d ', vUnits(u2)));
            set(hTit, 'color', FV.mColors(u2,:), 'fontsize', 8, 'fontweight', 'bold', 'backgroundcolor', [.1 .1 .1])
        end
        xlabel('ms'); ylabel('')
    end
end
hHeader = header(['Channel ' sCh ' - Spiketrain Cross Correlations'], 12);
set(hHeader, 'color', 'w', 'interpreter', 'none')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotCrossCorrelationsContinuous(varargin)
[FV,hWin] = GetStruct;

% TODO
% Plot cross correlations between all continuous data traces; rows vs axes
csChannels = FV.csDisplayChannels;

hFig = figure;
set(hFig, 'color', [.2 .2 .2], 'name', 'Spiky! Continuous Cross Correlations', 'NumberTitle', 'off')

% Iterate over channels and generate subplot
for ch1 = 1:length(csChannels)
    vTrace1 = FV.tData.(csChannels{ch1});
    nFs1 = FV.tData.([csChannels{ch1} '_KHz']) * 1000; % sampling frequency (Hz)
    nBeginTime1 = FV.tData.([csChannels{ch1} '_TimeBegin']); % start of sampling (sec)
    vTime1 = (nBeginTime1+1/nFs1):(1/nFs1):(nBeginTime1+length(vTrace1)/nFs1); % time vector (sec)

    for ch2 = 1:length(csChannels)
        vTrace2 = FV.tData.(csChannels{ch2});
        nFs2 = FV.tData.([csChannels{ch2} '_KHz']) * 1000; % sampling frequency (Hz)
        nBeginTime2 = FV.tData.([csChannels{ch2} '_TimeBegin']); % start of sampling (sec)
        vTime2 = (nBeginTime2+1/nFs2):(1/nFs2):(nBeginTime2+length(vTrace2)/nFs2); % time vector (sec)

        % Abort if the two traces dont have same length, Fs or BeginTime
        if (length(vTrace1) ~= length(vTrace2)) ... % length
                || (nFs1 ~= nFs2) ... % Fs
                || (nBeginTime1 ~= nBeginTime2) % BeginTime
            warndlg('Two traces dont start at same time or have different sampling rates. Cant cross correlate these')
            continue
        end

        % axes
        nX = .08; nY = .1; nW = .88; nH = .85;
        nNN = length(csChannels);
        hAx = axes('position', [nX+(ch2-1)*(nW/nNN) (1-nY)-ch1*(nH/nNN) nW/nNN nH/nNN], 'Color', [.1 .1 .1]);

        % Cross correlate
        nMaxLagSec = .1; % sec
        nMaxLagSamp = nFs1/(1/nMaxLagSec); % samples
        [vC, vLags] = xcorr(vTrace1, vTrace2, nMaxLagSamp, 'coeff');

        % X axis values
        vX = -(nMaxLagSamp/nFs1):(1/nFs1):(nMaxLagSamp/nFs1);

        % Plot
        if ch1 == ch2
            plot(vX, vC, 'w', 'linewidth', 1.5)
            set(hAx, 'Color', [.25 .25 .25])
        else
            plot(vX, vC, 'y')
            set(hAx, 'Color', [.1 .1 .1])
        end
        set(hAx, 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'fontsize', 7, 'ylim', [-.5 1] )

        % Labels above subplots (first row only)
        if ch1 == 1
            hTit = title(sprintf('%s', csChannels{ch2}));
            set(hTit, 'color', FV.mColors(ch2,:), 'fontsize', 8, 'fontweight', 'bold', 'backgroundcolor', [.1 .1 .1],'interpreter','none')
        else
            % if we're not in the first OR last row , remove x ticks
            if ch1 ~= length(csChannels)
                set(hAx, 'xticklabel', [])
            end
        end

        % Labels next to subplots (first column only)
        if ch2 == 1
            hTit = ylabel(sprintf('%s', csChannels{ch1}));
            set(hTit, 'color', FV.mColors(ch1,:), 'fontsize', 8, 'fontweight', 'bold', 'backgroundcolor', [.1 .1 .1],'interpreter','none')
        else
            % if we're not in the first column, remove y ticks
            set(hAx, 'yticklabel', [])
        end

        % If we're in the last row, add x label (time)
        if ch1 == length(csChannels)
            xlabel('s')
        end

        drawnow
    end
end

linkaxes(get(hFig,'children'),'xy')
zoom xon

hHeader = header(['Channel Cross Correlations'], 12);
set(hHeader, 'color', 'w', 'interpreter', 'none')

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowEventTriggeredAverage(varargin)
[FV,hWin] = GetStruct;

% Select trigger event
[sEventCh, bResult] = SelectChannelNumber(FV.csDigitalChannels', 'Select trigger event');
if ~bResult, return, end

% Select continuous channel
vIndx = [];
for ch = 1:length(FV.csChannels)
    if isempty(find(strcmp(FV.csChannels(ch), FV.csDigitalChannels), 1))
        vIndx(end+1) = ch;
    end
end
[sContCh, bResult] = SelectChannelNumber(FV.csChannels(vIndx)', 'Select continuous signal');
if ~bResult, return, end

% Get event up times
vEventTimes = FV.tData.([sEventCh '_Up']); % sec, abs time

% Error handling
% i) Check that we have any events
if isempty(vEventTimes)
    waitfor(warndlg(sprintf('No event triggers were detected.\nCannot display event triggered average.')));
    return
end
% ii)Check that we have data
if ~isfield(FV.tData, [sContCh '_KHz'])
    waitfor(warndlg(sprintf('Missing data field: %s. Analysis aborted.', [sContCh '_KHz'])));
    return
end

% Get continuous signal
nContFs = FV.tData.([sContCh '_KHz']) * 1000; % Hz
vCont = FV.tData.(sContCh);
vContBegin = FV.tData.([sContCh '_TimeBegin']);

% Pre/post times and stimulus delay
persistent p_nStimDel p_nPre p_nPost p_nDetrend p_nFirstPulse p_nLastPulse p_nDerivative
persistent p_bAbsolute p_bLowPassHz p_bHiPassHz p_bHilbert
if isempty(p_nStimDel), p_nStimDel = 0; end
if isempty(p_nPre), p_nPre = 1; end
if isempty(p_nPost), p_nPost = 2; end
if isempty(p_nDetrend), p_nDetrend = 0; end
if isempty(p_nFirstPulse), p_nFirstPulse = 1; end
if isempty(p_nLastPulse), p_nLastPulse = length(vEventTimes); end
if isempty(p_nDerivative), p_nDerivative = 0; end
if isempty(p_bAbsolute), p_bAbsolute = 0; end
if isempty(p_bHiPassHz), p_bHiPassHz = 0; end
if isempty(p_bLowPassHz), p_bLowPassHz = 1000; end
if isempty(p_bHilbert), p_bHilbert = 0; end

cPrompt = {'Pre-event duration (s)','Post-event duration (s)', 'Detrend (1=yes, 0=no)', ...
    'Stimulus delay (ms)', 'First pulse', sprintf('Last pulse (max %d)', length(vEventTimes)), ...
    'Derivative', 'Use absolute values (1=yes, 0=no)', 'High pass (Hz)', 'Low pass (Hz)', 'Average Hilbert amplitude (1=yes, 0=no)'};

cAnswer = inputdlg(cPrompt,'Averaging options', 1, ...
    {num2str(p_nPre), num2str(p_nPost), num2str(p_nDetrend), num2str(p_nStimDel), ...
    num2str(p_nFirstPulse), num2str(min([p_nLastPulse length(vEventTimes)])), num2str(p_nDerivative), ...
    num2str(p_bAbsolute), num2str(p_bHiPassHz), num2str(p_bLowPassHz), num2str(p_bHilbert)});
if isempty(cAnswer), return, end
p_nPre  = str2num(cAnswer{1}); % sec
p_nPost  = str2num(cAnswer{2}); % sec
p_nDetrend  = str2num(cAnswer{3});
p_nStimDel = str2num(cAnswer{4}); % ms
p_nFirstPulse = str2num(cAnswer{5});
p_nLastPulse = str2num(cAnswer{6});
p_nDerivative = str2num(cAnswer{7});
p_bAbsolute = str2num(cAnswer{8});
p_bHiPassHz = str2num(cAnswer{9});
p_bLowPassHz = str2num(cAnswer{10});
p_bHilbert = str2num(cAnswer{11});

% Adjust event times for stimulus delay
vEventTimes = vEventTimes + (p_nStimDel/1000); % sec

% Event times relative to start of continuous signal
vRelTimes = vEventTimes - vContBegin;

% Get sample numbers where events occured in continuous signal
vOnsets = round(vRelTimes * nContFs); % samples

% Derivate
if p_nDerivative > 0
    vCont = diff(vCont, p_nDerivative);
end

% Low/Hi pass and rectify
vTime = linspace(0, (1/nContFs)*length(vCont), length(vCont));
[vCont, vTime2, nContFs] = FilterChannel(vCont, vTime, nContFs, p_bLowPassHz, p_bHiPassHz, p_bAbsolute, 'none');

% Convert unit to Intensity/s^-d
vCont = vCont * (nContFs^p_nDerivative);

% Compute event triggered average
SpikyWaitbar(0, 20);
nLen = length(vCont);
mAll = [];
for o = p_nFirstPulse:p_nLastPulse
    SpikyWaitbar((o/length(vOnsets))*20, 20);

    %nStimStart = (nLen - vOnsets(o));
    nStart = vOnsets(o) - round(p_nPre * nContFs);
    nEnd   = vOnsets(o) + round(p_nPost * nContFs);
    if nStart < 1 || nEnd > nLen
        continue
    end
    vThis = vCont(nStart:nEnd);
    
    % Detrend
    if p_nDetrend, vThis = detrend(vThis); end

    % Compute the hilbert transform and
    if p_bHilbert
        % Subtract set-point (< 2 Hz)
        vSetPoint = filter_series(double(vThis(:)), nContFs, 2);
        vThis = vThis - vSetPoint';
        
        % Replace NaNs with 0
        vNaNIndx = find(isnan(vThis));
        vThis(vNaNIndx) = 0;
        
        % Hilbert transform
        vHilb = hilbert(vThis);
        
        % Hilbert magnitude and phase
        vThis = abs(vHilb);
    end
    
    mAll(end+1, :) = vThis;
end
SpikyWaitbar(20, 20);

% Compute statistics
vMean = nanmean(mAll);
vMedian = nanmedian(mAll);
vErrMean = nanstd(mAll) ./ sqrt(size(mAll, 1));
vErrMedian = vErrMean * 1.25;

vTime = linspace(-p_nPre, p_nPost, size(mAll, 2)); % sec

% Initialize figure and plot
hFig = figure('color', [.2 .2 .2], 'units', 'pixels');
centerfig(hFig, hWin)
set(hFig, 'name', 'Spiky! Event Triggered Average', 'NumberTitle', 'off');
hAx = axes();
[hMeanLine, hMeanErr] = mean_error_plot(vMean, vErrMean, [0 0 1], vTime);
hold on
[hMedianLine, hMedianErr] = mean_error_plot(vMedian, vErrMedian, [0 1 0], vTime);
set([hMedianLine, hMedianErr], 'visible', 'off', 'tag', 'ETAMedian');
set([hMeanLine, hMeanErr], 'visible', 'on', 'tag', 'ETAMean');

% Plot all individual traces
hold on
hTrials = plot(repmat(vTime, size(mAll, 1), 1)', mAll', 'w');
set([hTrials], 'visible', 'off', 'tag', 'AllTrials');

hCheck = uicontrol(hFig, 'Style', 'checkbox', 'Position', [1 1 100 20], 'String', 'Mean', ...
    'HorizontalAlignment', 'left', 'backgroundcolor', [.2 .2 .2], 'foregroundcolor', 'w', 'value', 1);
set(hCheck, 'Callback', 'if(get(gcbo,''value'')),V=''on'';else,V=''off'';end;set(findobj(''tag'',''ETAMean''),''visible'',V);')

hCheck = uicontrol(hFig, 'Style', 'checkbox', 'Position', [101 1 100 20], 'String', 'Median', ...
    'HorizontalAlignment', 'left', 'backgroundcolor', [.2 .2 .2], 'foregroundcolor', 'w', 'value', 0);
set(hCheck, 'Callback', 'if(get(gcbo,''value'')),V=''on'';else,V=''off'';end;set(findobj(''tag'',''ETAMedian''),''visible'',V);')

hCheck = uicontrol(hFig, 'Style', 'checkbox', 'Position', [201 1 100 20], 'String', 'Trials', ...
    'HorizontalAlignment', 'left', 'backgroundcolor', [.2 .2 .2], 'foregroundcolor', 'w', 'value', 0);
set(hCheck, 'Callback', 'if(get(gcbo,''value'')),V=''on'';else,V=''off'';end;set(findobj(''tag'',''AllTrials''),''visible'',V);')

% Axes properties
axis tight
set(hAx, 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'xlim', [-p_nPre p_nPost])
xlabel('Time (s)')
ylabel(sprintf('Intensity/s^%d', p_nDerivative))
title(sprintf('%s  N=%d events (range %d-%d)', sContCh, size(mAll,1), p_nFirstPulse, p_nLastPulse), ...
    'interpreter', 'none', 'color', [.7 .7 .7], 'FontWeight', 'bold')
zoom on
grid on
set(hFig, 'ToolBar', 'figure')

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotHistogramContinuous(varargin)
[FV,hWin] = GetStruct;
[sCh, bResult] = SelectChannelNumber(FV.csChannels);
hFig = figure('color', [.2 .2 .2]);
set(hFig, 'name', 'Spiky! Histogram', 'NumberTitle', 'off');
vData = FV.tData.(sCh);
% Remove extreme outliers (0.025% off either edge)
vData = sort(vData);
vData = vData((length(vData)/400):(length(vData)-(length(vData)/400)));
[vN, vX] = hist(vData, 100);
hAx = axes;
nCh = strcmpi({FV.tChannelDescriptions.sChannel}, sCh);
vCol = FV.mColors(nCh,:);
if isempty(vCol)
    vCol = [.8 .8 .8]; % default color, gray
end
hBar = bar(vX, vN);
set(hBar, 'faceColor', vCol, 'edgeColor', vCol)
set(hAx, 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'fontsize', 7)
xlabel('V')
ylabel('Number of samples')
hTit = title(sprintf('%s (%s)', FV.tChannelDescriptions(nCh).sDescription, sCh));
set(hTit, 'FontSize', 8, 'FontWeight', 'bold', 'color', vCol, 'backgroundcolor', [.1 .1 .1], 'interpreter', 'none')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotPSTH(varargin)
[FV,hWin] = GetStruct;
[sCh, bResult] = SelectChannelNumber(fieldnames(FV.tSpikes)');
if isempty(sCh), warndlg('No channels were found that contained spiking data.'); return; end

if isfield(FV.tSpikes.(sCh), 'hierarchy')
    % Analyze sorted units
    vUnits = unique(FV.tSpikes.(sCh).hierarchy.assigns); % unit names
else
    % Analyze un-sorted units
    vUnits = NaN;
end

% Find event indices
vEventIndx = [];
cFields = fieldnames(FV.tData);
% Ignore DAQ_Start, _Stop and _Trigger fields (default fields generated by DAQ Toolbox)
cIgnoreFields = {'DAQ_Start_Up', 'DAQ_Stop_Up', 'DAQ_Trigger_Up'};
for e = 1:length(cFields) % iterate over fieldnames
    if strfind(cFields{e}, '_Up') & ~ismember(cFields{e}, cIgnoreFields)
        vEventIndx(end+1) = e;
    end
end
if isempty(vEventIndx)
    waitfor(warndlg('No events to trigger PSTHs on were detected.'));
    return
end

% Create figure
hFig = figure('color', [.2 .2 .2]);
drawnow
set(hFig, 'name', 'Spiky! Peristimulus Time Histograms (PSTH)', 'NumberTitle', 'off');
vXLim = [-.1 .1];
nTrialLen = 100;

% Ask for PSTH parameters
persistent p_nStimDel p_nPreStimDur p_nPostStimDur p_nBinRes
if isempty(p_nStimDel), p_nStimDel = 0; end
if isempty(p_nPreStimDur), p_nPreStimDur = 0.05; end
if isempty(p_nPostStimDur), p_nPostStimDur = 0.15; end
if isempty(p_nBinRes), p_nBinRes = 0.001; end
cAnswer = inputdlg({'Stimulus delay (ms)','Pre-stim period (s)','Post-stim period (s)','Bin resolution (ms)'},...
    'Stimulus Delay', 1, {num2str(p_nStimDel), num2str(p_nPreStimDur), num2str(p_nPostStimDur), num2str(p_nBinRes)});
if isempty(cAnswer), return, end
p_nStimDel = str2num(cAnswer{1});
p_nPreStimDur = str2num(cAnswer{2});
p_nPostStimDur = str2num(cAnswer{3});
p_nBinRes = str2num(cAnswer{4});

% Check that events have enough data
vRemIndx = [];
for e = 1:length(vEventIndx) % iterate over fieldnames
    if isempty(strfind(cFields{vEventIndx(e)}, '_Up'))
        vRemIndx(end+1) = e;
        continue
    end
    vUpTimes = FV.tData.(cFields{vEventIndx(e)}); % sec, abs time
    if (length(vUpTimes) < 2)
        vRemIndx(end+1) = e;
    end
end
vEventIndx(vRemIndx) = [];

% Iterate over units
nRow = 1;
for u = 1:length(vUnits)
    % Iterate over events
    nCol = 1;
    for e = 1:length(vEventIndx) % iterate over fieldnames
        % Get event up times
        %nFs = FV.tData.([cFields{vEventIndx(e)}(1:end-2) 'KHz']) * 1000;
        nFs = FV.tSpikes.(sCh).Fs;
        vUpTimes = FV.tData.(cFields{vEventIndx(e)}); % sec, abs time
        
        % Skip channels that have more than 10,000 events
        if length(vUpTimes) > 10000
            uiwait(warndlg(sprintf('Channel %s has more than 10,000 events and will be skipped.', cFields{vEventIndx(e)}), 'Spiky!'))
            continue
        end
        
        % Get spiketimes
        if isnan(vUnits(u))
            % Un-sorted unit
            vSpiketimes = FV.tSpikes.(sCh).spiketimes(:) ./ nFs; % sec
        else
            % Sorted unit
            vIndx = FV.tSpikes.(sCh).hierarchy.assigns == vUnits(u);
            vSpiketimes = FV.tSpikes.(sCh).spiketimes(vIndx) ./ nFs; % sec
        end

        % Subtract stimulus delay from spiketimes
        vSpiketimes = vSpiketimes - (p_nStimDel/1000); % sec
        
        % Iterate over event times
        vPSTH = [];
        for et = 1:length(vUpTimes)-1
            vIndx = vSpiketimes >= (vUpTimes(et)-p_nPreStimDur) & vSpiketimes < vUpTimes(et+1);
            vRelTimes = vSpiketimes(vIndx) - vUpTimes(et);
            vPSTH = [vPSTH vRelTimes'];
        end

        if ~ishandle(hFig) return; end
        
        % Plot PSTH
        nX = .08; nY = .1; nW = .88; nH = .80;
        nNN = length(vEventIndx); nNNY = length(vUnits);
        hAx = axes('position', [nX+(nCol-1)*(nW/nNN) (1-nY)-u*(nH/nNNY) nW/nNN nH/nNNY]);

        if vUnits(u) == 0, vCol = [.5 .5 .5]; % outliers
        else vCol = FV.mColors(u,:); end

        nYMax = median(diff(vUpTimes));

        vPSTH(vPSTH>nYMax) = [];
        [vC, vT] = hist(vPSTH, -p_nPreStimDur:p_nBinRes:nYMax);
        vC = (vC./et) * (1/p_nBinRes); % normalize to spikes/sec

        % Convolve with a right-angled triangular window
        if ~isempty(vC)
            nLen = 5; % base width 5 ms
            vWin = [linspace(1,0,nLen)];
            vWin = vWin/sum(vWin); % normalize to area 1
            vC = conv(vC, vWin);
        end

        % Plot as line
        %plot(vT, vC(1:end-(nLen-1)), 'color', vCol)
        
        % Plot as bars
        if ~isempty(vT)
            hBar = bar(vT, vC(1:end-(nLen-1)));
            set(hBar, 'facecolor', vCol, 'edgecolor', vCol)
        end
        
        if isnan(nYMax) || nYMax < 0, nYMax = 0.1; end
        set(hAx, 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'fontsize', 7, 'xlim', [-p_nPreStimDur p_nPostStimDur])
        box on

        if nCol == 1
            ylabel('Spikes/s')
            if isnan(vUnits(u))
                % Un-sorted unit
                hTxt = text(0,0,[' Unit UN-SORTED']);
            else
                % Sorted unit
                hTxt = text(0,0,sprintf(' Unit %d', vUnits(u)));
            end
            set(hTxt, 'color', FV.mColors(u,:), 'fontsize', 8, 'fontweight', 'bold', 'backgroundcolor', [.1 .1 .1], ...
                'interpreter', 'none', 'HorizontalAlignment', 'center', 'units', 'normalized', ...
                'backgroundcolor', [.1 .1 .1], 'position', [-.1 .5 0], 'Rotation', 90)
        end

        if nRow == 1
            sChName = cFields{vEventIndx(e)}(1:end-3);
            if isfield(FV, 'tChannelDescriptions')
                nIndx = find(strcmpi({FV.tChannelDescriptions.sChannel}, sChName));
            else nIndx = []; end
            if isempty(nIndx), hTit = title(sChName);
            else
                hTit = title(sprintf('%s (%s)', FV.tChannelDescriptions(nIndx).sDescription, sChName));
            end
            set(hTit, 'FontSize', 8, 'FontWeight', 'bold', 'color', [1 1 0], 'backgroundcolor', [.1 .1 .1], 'interpreter', 'none')
        elseif u == length(vUnits)
            xlabel('Time (s)');
        end
        
        nCol = nCol + 1; % increment column counter
        drawnow
    end
nRow = nRow + 1; % increment row counter
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = ChannelProperties(varargin)
if ~CheckDataLoaded, return, end
[FV, hWin] = GetStruct;
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
else nTimeBegin = NaN; end

% Descriptive string
nIndx = find(strcmp({FV.tChannelDescriptions.sChannel}, sCh));
if ~isempty(nIndx)
    sDescr = FV.tChannelDescriptions(nIndx).sDescription;
else sDescr = ''; end

cAboutText{1} = sprintf('Name: %s', sCh);
cAboutText{2} = sprintf('Description: %s', sDescr);
cAboutText{3} = sprintf('Sampling rate: %.2f kHz', nFs);
cAboutText{4} = sprintf('Length: %d samples', nLen);
cAboutText{5} = sprintf('Gain: %d', nGain);
cAboutText{6} = sprintf('Begin time: %.4f s', nTimeBegin);

hBox = msgbox(cAboutText, 'Properties');

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = DigitizeChannelAuto(hObject, varargin)
DigitizeChannel(hObject, 'auto')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = DigitizeChannel(varargin)

% Check if data and threshold was supplied with function inputs
if length(varargin) == 5 && isnumeric(varargin{1})
    bInteractive = 0;
elseif length(varargin) == 2 && strcmp('auto', varargin{2})
    bInteractive = 1;
    nY = .5; % threshold for automatic digitization set at 0.5 V
else
    bInteractive = 1; % interactive mode; user clicks to set threshold
end

if bInteractive
    % Get data and threshold interactively
    [FV, hWin] = GetStruct;
    if ~CheckDataLoaded, return, end

    % Select channel
    if ~exist('nY')
        zoom off;
        [nX nY] = ginput(1);
        zoom xon
    end
    sTag = get(gca, 'Tag');
    vContData = FV.tData.(sTag);
    nFs = FV.tData.([sTag '_KHz']) * 1000; % sampling frequency (Hz)
    nBeginTime = FV.tData.([sTag '_TimeBegin']); % start of sampling (sec)
    nEndTime = FV.tData.([sTag '_TimeEnd']); % start of sampling (sec)

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
    vUpTimes = vUpTimes - ((1/nFs)/2);
    vDownTimes = vDownTimes - ((1/nFs)/2);
else
    vUpTimes = [];
    vDownTimes = [];
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
    else nIndx = length(FV.tEventThresholds) + 1; end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DuplicateChannel(varargin)
% Duplicate a channel with options for resampling and filtering
[FV, hWin] = GetStruct;
if ~CheckDataLoaded, return, end

[sCh, bResult] = SelectChannelNumber(FV.csChannels);

if isempty(sCh), return, end

% Select channel
vCont = FV.tData.(sCh);
nFs = FV.tData.([sCh '_KHz']) * 1000;

% Decimate and filter
cPrompt = {'High-pass (Hz)','Low-pass (Hz)', 'Rectify (1=yes, 0=no)'};
cAnswer = inputdlg(cPrompt,'Resampling options',1, {'0', num2str(nFs), num2str(FV.nEMGRectify)});
if isempty(cAnswer), return, end

nHiPass = str2num(cAnswer{1});
nLoPass = str2num(cAnswer{2});
bRectify = str2num(cAnswer{3});

if ~(nHiPass == 0 && nLoPass == nFs)
    [vContNew, vTime, nFsNew] = FilterChannel(vCont, [], nFs, nLoPass, nHiPass, bRectify, 'decimate');
else vContNew = vCont; end

FV.tData.([sCh '_Copy_Imported']) = 1;
FV.tData.([sCh '_Copy']) = vContNew;
FV.tData.([sCh '_Copy_KHz']) = nFsNew / 1000;
FV.tData.([sCh '_Copy_KHz_Orig']) = nFsNew / 1000;
FV.tData.([sCh '_Copy_TimeBegin']) = FV.tData.([sCh '_TimeBegin']);
FV.tData.([sCh '_Copy_TimeEnd']) = FV.tData.([sCh '_TimeEnd']);

FV.csDisplayChannels = unique([FV.csDisplayChannels [sCh '_Copy']]);
FV.csChannels = unique([FV.csChannels [sCh '_Copy']]);

SetStruct(FV)
ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DeleteChannel(varargin)
% Delete a channel
[FV, hWin] = GetStruct;
if ~CheckDataLoaded, return, end

csAll = [FV.csDisplayChannels FV.csChannels];
[sCh, bResult] = SelectChannelNumber(csAll);
%[sCh, bResult] = SelectChannelNumber(FV.csChannels);
if isempty(sCh), return, end % no channel was selected (i.e. dialog was closed)
sAns = questdlg(['Are you sure you want to delete channel ' sCh '?'], 'Spiky!', 'Yes', 'No', 'No');
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

SetStruct(FV)
ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UndigitizeChannel(varargin)
% Undigitizes channels that were previously digitized through the GUI
[FV, hWin] = GetStruct;
if ~CheckDataLoaded, return, end

[sCh, bResult] = SelectChannelNumber(FV.csDigitalChannels);
if isempty(sCh), return, end % no channel was selected (i.e. dialog was closed)

if ~isfield(FV.tData, sCh)
    % Check whether channel has an analog equivalent. If not, we cannot 'undigitize'it
    warndlg(['Channel ' sCh ' cannot be un-digitized as it has no analog equivalent.'], 'Spiky!')
    return
end

% 'Un-digitize' channel
FV.csDigitalChannels(strcmpi(FV.csDigitalChannels, sCh)) = [];
SetStruct(FV)
ViewTrialData
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sFilename = CheckFilename(sFilename)
if ispc
    sFilename = strrep(sFilename, '/', '\');
    sFilename = strrep(sFilename, '\\', '\');
elseif isunix
    sFilename = strrep(sFilename, '\', '/');
    sFilename = strrep(sFilename, '//', '/');
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bResult = SpikyWaitbar(nStatus, nLen)
% Usage: SpikyWaitbar(0, 20)    Initialize waitbar with length 20
%        SpikyWaitbar(10, 20)   Waitbar reaches 10 of 20
%        SpikyWaitbar(20, 20)   Waitbar reaches end and is destroyed
%        SpikyWaitbar(0, 0)     Function returns current position (0 if no bar exists)
% Function returns false if the Cancel button is pressed
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
drawnow
hCancelButton = findobj(g_hSpike_Viewer, 'Tag', 'Spiky_Waitbar_Cancel');
if ~isempty(hCancelButton)
    if strcmp(get(hCancelButton(1), 'State'), 'on')
        bResult = false; % process running waitbar was cancelled
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChannelDescriptions(varargin)
[FV, hWin] = GetStruct;
if ~CheckDataLoaded, return, end

hFig = figure('closeRequestFcn', 'set(gcbf,''userdata'',1)');
cColumnNames = {'Channel', 'Description'};
cColumnFormat = {{'numeric' 'Adjustable'}, {'numeric' 'Adjustable'}};
cColumnEditable =  [false true];

% create list of selectable channels
cFields = fieldnames(FV.tData);
csDisplayChannels = {};
cData = {};
for i = 1:length(cFields)
    if ~isempty(strfind(cFields{i}, '_TimeBegin'))
        cData{end+1, 1} = cFields{i}(1:end-10);
        sDescr = '';
        if isfield(FV, 'tChannelDescriptions')
            if ~isfield(FV.tChannelDescriptions, 'sChannel'), nIndx = [];
            else
                nIndx = find(strcmp({FV.tChannelDescriptions.sChannel}, cData{end, 1}));
            end
            % If current channel does not appear in the tChannelDescriptions structure, create it
            if isempty(nIndx)
                FV.tChannelDescriptions(end+1).sChannel = cData(end, 1);
                FV.tChannelDescriptions(end+1).sDescription = '';
                nIndx = length(FV.tChannelDescriptions);
            end
            sDescr = FV.tChannelDescriptions(nIndx).sDescription;
        end
        cData{end, 2} = sDescr;
    end
end

set(hFig, 'color', [.2 .2 .2], 'Name', 'Spiky! Channel Descriptions', 'NumberTitle', 'off', 'ToolBar', 'none', 'menuBar','none')
hTable = uitable('Units', 'normalized','Position', [0 0 1 1], 'Data', cData, 'ColumnName', ...
    cColumnNames, 'ColumnEditable', cColumnEditable, 'ColumnWidth',{150, 377});
waitfor(hFig, 'userdata') % wait for figure to be closed
cData = get(hTable, 'data'); % modified data

if ~isfield(FV, 'tChannelDescriptions')
    FV.tChannelDescriptions = struct([]);
end
for c = 1:size(cData, 1)
    FV.tChannelDescriptions(c).sChannel = cData{c, 1};
    FV.tChannelDescriptions(c).sDescription = cData{c, 2};
end

delete(hFig)
SetStruct(FV)
ViewTrialData()
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExperimentDescriptions(varargin)
[FV, hWin] = GetStruct;
if ~CheckDataLoaded, return, end

hFig = figure('closeRequestFcn', 'set(gcbf,''userdata'',''closed'')');
cColumnNames = {'Variable', 'Value'};
cColumnFormat = {{'char' 'Adjustable'}, {'char' 'Adjustable'}};
cColumnEditable =  [true true];

% create list of selectable channels
if ~isfield(FV, 'tExperimentVariables')
    FV_temp = SetFVDefaults();
    FV.tExperimentVariables = FV_temp.tExperimentVariables;
end
cData = repmat({''},200,2);;
for i = 1:length(FV.tExperimentVariables)
    cData{i, 1} = FV.tExperimentVariables(i).sVariable;
    cData{i, 2} = FV.tExperimentVariables(i).sValue;
end

set(hFig, 'color', [.2 .2 .2], 'Name', 'Spiky! Experiment Variables', 'NumberTitle', 'off', 'ToolBar', 'none', 'menuBar','none')
hTable = uitable('Units', 'normalized','Position', [0 0 1 1], 'Data', cData, 'ColumnName', ...
    cColumnNames, 'ColumnEditable', cColumnEditable, 'ColumnWidth',{150, 352});
waitfor(hFig, 'userdata', 'closed') % wait for figure to be closed
cData = get(hTable, 'data'); % modified data

FV.tExperimentVariables = struct([]); % clear variables
for c = 1:size(cData, 1)
    if isempty(cData{c,1}) | isempty(cData{c,2}), continue, end
    FV.tExperimentVariables(end+1).sVariable = cData{c, 1};
    FV.tExperimentVariables(end).sValue = cData{c, 2};
end

delete(hFig)
SetStruct(FV)
ViewTrialData()
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DeleteSortingData(varargin)
[FV, hWin] = GetStruct;
if ~CheckDataLoaded, return, end

sAns = questdlg('All sorting data, including assignment of spikes as outliers, will be removed. Waveforms and dejittering is not affected. Continue?', ...
    'Spiky!', 'Yes', 'No', 'No');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust channel amplitude by specified gain. Returned gain is in units of mV
function [vCont, FV] = AdjustChannelGain(FV, vCont, sCh) % mV

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
    uiwait(warndlg(['No gain has been specified for ' GetChannelDescription(sCh) '. Adjust gain values in the next dialog window.']))
    bResetGain = 0;
    if isempty(FV.tGain), bResetGain = 1;
    elseif ~isfield(FV.tGain, sCh), bResetGain = 1; end
    if bResetGain, FV.tGain(1).(sCh) = 1; end % default gain of 1
    SetStruct(FV)
    SetGain('noplot')
    [FV,hWin] = GetStruct; % get updated values
    nGain = FV.tGain.(sCh);
end

vCont = vCont ./ nGain; % divide by gain

% Multiply with 1000 so units to from V to mV
if nGain ~= 1
    vCont = vCont .* 1000; % mV
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get description of passed channel name
function sDescr = GetChannelDescription(sCh)
[FV,hWin] = GetStruct;
sDescr = sCh; % default; description same as name
if isfield(FV, 'tChannelDescriptions')
    if isfield(FV.tChannelDescriptions, 'sChannel')
        nIndx = find(strcmp({FV.tChannelDescriptions.sChannel}, sCh));
        if ~isempty(FV.tChannelDescriptions(nIndx))
            if ~isempty(FV.tChannelDescriptions(nIndx).sDescription)
                sDescr = FV.tChannelDescriptions(nIndx).sDescription;
            end
        end
    end
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get current SVN build number
function sBuild = GetBuildNumber

% Set paths
sWT_dir = which('spiky');
sSVNPath = CheckFilename([sWT_dir(1:findstr(sWT_dir, 'spiky.m')-1) '.svn\' 'entries']);

% Get build number from SVN entries files
sBuild = 'Unknown';
if exist(sSVNPath, 'file')
    fid = fopen(sSVNPath);
    for i = 1:5
        sLine = fgetl(fid);
        if strcmp(sLine, 'dir'), sBuild = fgetl(fid); end
    end
    fclose(fid);
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Invert a selected channel
function InvertChannel(varargin)
[FV, hWin] = GetStruct;
if ~CheckDataLoaded, return, end

% Select channel
[sCh, bResult] = SelectChannelNumber(FV.csChannels);

% Invert channel data
FV.tData.(sCh) = FV.tData.(sCh) .* -1;

% Force inverted channel to be displayed
FV.csDisplayChannels = unique([FV.csDisplayChannels sCh]);

SetStruct(FV)
ViewTrialData

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Undo
% TODO
function Undo(varargin)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Redo
% TODO
function Redo(varargin)
keyboard
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show normal time (i.e. start of trial is zero)
function NormalTime(varargin)
[FV, hWin] = GetStruct;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turn on/off x axis grid
function ShowGrid(varargin)
[FV, hWin] = GetStruct;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shift BeginTime time of selected channel
function ShiftBeginTime(varargin)
[FV, hWin] = GetStruct;

% Check that tag of handle references a known channel
sCh = get(varargin{1}, 'Tag');
if ~(isempty(sCh) | ~any(strcmp(fieldnames(FV.tData), sCh)))
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show info stored in DAQ file
function ShowDAQInfo(varargin)
[FV, hWin] = GetStruct;
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
warndlg('No DAQ file loaded or loaded file is not a DAQ (.daq) file.', 'Spiky!')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MeasureLine(varargin)
hObj = varargin{1}; % get calling object

% Find axis with same tag and make that current
hAx = findobj('Type', 'axes', 'Tag', get(hObj, 'Tag'));
if isempty(hAx), return
else axes(hAx(1)), end

mData = get(get(hObj,'parent'),'userdata');
vCont = mData(:,2);
vTime = mData(:,1);

hold on

% Point A
[nXa nYa] = ginput(1);
% Snap to nearest data point on time axis
[nMinVal nMinIndx] = min(abs(vTime - nXa));
nXa_samples = nMinIndx;
nXa = vTime(nMinIndx);
nYa = vCont(nMinIndx);
hPnta = plot(nXa, nYa, 'ro');

% Point B
[nXb nYb] = ginput(1);
% Snap to nearest data point on time axis
[nMinVal nMinIndx] = min(abs(vTime - nXb));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the number of UP/DOWN events on each channel
function ShowEventStatistics(varargin)
[FV, hWin] = GetStruct;
% iterate over digital channels
csChannels = FV.csDigitalChannels;
hFig = figure();
cColumnNames = {'Name', 'kHz', '#Up', '#Down', 'Description'};
cColumnFormat = {{'char' 'Fixed'}, {'numeric' 'Fixed'}, {'numeric' 'Fixed'}, {'numeric' 'Fixed'}, {'char' 'Fixed'}};
cColumnEditable =  [false false false false false];
for nCh = 1:length(csChannels)
    cData{nCh, 1} = csChannels{nCh};
    cData{nCh, 2} = FV.tData.([csChannels{nCh} '_KHz']);
    cData{nCh, 3} = length(FV.tData.([csChannels{nCh} '_Up']));
    cData{nCh, 4} = length(FV.tData.([csChannels{nCh} '_Down']));
    % Get average event duration
    if cData{nCh, 3} == cData{nCh, 4}
        FV.tData.([csChannels{nCh} '_Down']) - FV.tData.([csChannels{nCh} '_Up']);
    end
    nIndx = strcmp(csChannels{nCh}, {FV.tChannelDescriptions.sChannel});
    if ~isempty(FV.tChannelDescriptions(nIndx))
        cData{nCh, 5} = FV.tChannelDescriptions(nIndx).sDescription;
    end
end
set(hFig, 'color', [.2 .2 .2], 'Name', 'Spiky! Event Statistics', 'NumberTitle', 'off', 'ToolBar', 'none', 'menuBar','none')
hTable = uitable('Units', 'normalized','Position', [0 0 1 1], 'Data', cData, 'ColumnName', ...
    cColumnNames, 'ColumnEditable', cColumnEditable, 'ColumnWidth',{100, 100, 70, 70, 180});
centerfig(hFig, hWin)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete section of selected channel. This function is only called by the context
% menu which appears when right-clicking channel traces
function DeleteSection(varargin)
[FV, hWin] = GetStruct;

sCh = get(gcbo, 'tag'); % channel name
hAx = findobj(gcf, 'type', 'axes', 'tag', sCh); % find corresponding axes
hAx = hAx(1);
if isempty(hAx)
    warndlg('Cannot find corresponding axes for selected channel. Aborting...', 'Spiky!')
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
        'Spiky!', 'Delete', 'Cancel', 'Cancel')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClassifyUnit(varargin)
% Classify a selected unit. This function assumes that the calling object
% (gcbo) identifies the channel and unit number in the UserData property
tUserData = get(gcbo, 'UserData');
[FV, hWin] = GetStruct;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetElectrodePosition(varargin)
% Set and store the electrode position used during the experiment.
sCh = get(gcbo, 'tag');
[FV, hWin] = GetStruct;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetReceptiveField(varargin)
% Set and store the electrode position used during the experiment.
[FV, hWin] = GetStruct;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BatchRedo(varargin)
% Set and store the electrode position used during the experiment.
global g_hSpike_Viewer
sAction = varargin{2};
persistent sLastAction
global g_bBatchMode

% List of files
hFileList = findobj(g_hSpike_Viewer, 'Tag', 'MenuFileList');
csFiles = flipud(get(get(hFileList, 'Children'), 'UserData'));
nFiles = length(csFiles);

if strcmp(sAction, 'redo')
    % Redo last recorded BATCH action (B)
    if isempty(sLastAction) % Do nothing
        warndlg('There is no available action to run in batch mode. Note that only actions with the B suffix (e.g. in menus or buttons) can run as batch jobs. To start a batch job, you must first run one of the supported actions and then re-select Batch Redo from the menu.')
    else % Repeat last supported redo action
        sAns = questdlg(['The last supported action was ' sLastAction '. Do you want to repeat this function on (' num2str(nFiles) ') currently loaded files?'], 'Spiky!');
        sAns2 = questdlg('Should the batch operation be repeated on the currently loaded file?', 'Spiky!', 'Yes', 'No', 'No');

        switch sAns
            case 'Yes' % Load each movie, redo last action and save result
                bWaitResult = SpikyWaitbar(0,nFiles+1);
                g_bBatchMode = true;
                % Get index of current file
                [FV, hWin] = GetStruct;
                sLoadedTrial = FV.sLoadedTrial;
                for m = 1:nFiles
                    % Progress bar
                    bWaitResult = SpikyWaitbar(m,nFiles+1);
                    if ~bWaitResult break, end % cancel button pressed
                    % If applicable, skip current file
                    if strcmp(csFiles{m}, sLoadedTrial) && strcmp(sAns2, 'No') continue; end
                    OpenFile([], m)     % Load movie
                    eval(sLastAction)   % Redo action
                    SaveResults();      % Save results
                end
                if bWaitResult % Batch complete notification
                    sStr = sprintf('Batch job %s completed (%d files processed)', sLastAction, nFiles+1);
                else
                    sStr = sprintf('Batch job %s cancelled after file %d/%d', sLastAction, m, nFiles+1);
                end
                bWaitResult = SpikyWaitbar(nFiles+1,nFiles+1);
                uiwait(msgbox(sStr, 'Spiky! Batch Redo', 'modal'));
                sp_disp(sStr);
            case 'No'       % Do nothing
            case 'Cancel'   % Do nothing
        end
    end
else % Record which action was last performed
    sLastAction = sAction;
end
g_bBatchMode = false;
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sp_disp(sStr)
disp(sprintf('%s Spiky says: %s', datestr(now, 'HH:mm:ss'), sStr))
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RunBatchScript(varargin)
RunScript('batch')
SaveResults();
BatchRedo([], 'redo')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RunScript(varargin)
[FV, hWin] = GetStruct;
global g_bBatchMode

% Get default path of ./scripts
sPath = which('spiky');
sPath = checkfilename([sPath(1:end-7) 'scripts\']);

% If the function to be run was not specified as in input, select file manually
if length(varargin) ~= 1
    [sFile, sPath] = uigetfile('*.m', 'Pick an M-file', sPath);
elseif strcmpi(varargin{1}, 'batch')
    % Run script as batch job

    % Get script file
    [sFile, sPath] = uigetfile('*.m', 'Pick an M-file', sPath);
    
    % Save script action in Batch Redo function
    if ~g_bBatchMode BatchRedo([], sprintf('RunScript(''%s'')', sFile)); end
    return
else
    sFile = varargin{1};
end

% The remainder of this function runs if a script is run for a single movie

% Save action for later global Redo
if ~g_bBatchMode BatchRedo([], sprintf('RunScript(''%s'')', sFile)); end

if sFile == 0, return, end
sCurrDir = pwd;
cd(sPath)
FV.ScriptError = '';
%try
    eval(['FV = ' sFile(1:end-2) '(FV);'])
%catch
%    uiwait(warndlg(sprintf('Script execution failed. There appears to be a problem with the script you are attempting to run:\n\n%s', lasterr), 'Spiky!'))
%    return
%end
if ~isempty(FV.ScriptError)
    sp_disp(sprintf('In %s: %s', sFile(1:end-2), FV.ScriptError))
end
cd(sCurrDir)
SetStruct(FV)

% Execute termination commands from script
if isfield(FV, 'ScriptExitCommand')
    eval(FV.ScriptExitCommand)
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%