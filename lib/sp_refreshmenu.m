function sp_refreshmenu()
% Create or refresh menu in Spiky GUI
% 
% Usage:
%   sp_refreshmenu()
%

global Spiky;

% Get Spiky path
sSpikyPath = which('spiky');
sPath = [sSpikyPath(1:end-7) 'icons/'];

% Remove existing menu
hMenu = findobj(Spiky.main.GetGUIHandle(), 'type', 'uimenu');
if ~isempty(hMenu)
    delete(hMenu)
end

hGUI = Spiky.main.GetGUIHandle();

% File menu
hFile  = uimenu(hGUI, 'Label', '&File');
uimenu(hGUI, 'Parent', hFile, 'Label', '&Open...', 'Callback', Spiky.main.OpenFile, 'Accelerator', 'O');
uimenu(hGUI, 'Parent', hFile, 'Label', 'Open &Directory...', 'Callback', Spiky.main.SetDirectory, 'Accelerator', 'D');
uimenu(hGUI, 'Parent', hFile, 'Label', 'Open Directory Tree...', 'Callback', Spiky.main.SetDirectoryTree, 'Accelerator', 'T');
uimenu(hGUI, 'Parent', hFile, 'Label', 'Open Settings...', 'Callback', Spiky.main.OpenSettings);
uimenu(hGUI, 'Parent', hFile, 'Label', 'Import...', 'Callback', Spiky.main.ImportFile, 'Accelerator', 'I');
uimenu(hGUI, 'Parent', hFile, 'Label', '&Save', 'Callback', Spiky.main.SaveResults, 'separator', 'on', 'Accelerator', 'S');
uimenu(hGUI, 'Parent', hFile, 'Label', '&Export...', 'Callback', Spiky.main.ExportData);
hFileList = uimenu(hGUI, 'Parent', hFile, 'Label', 'Files', 'Tag', 'MenuFileList', 'separator', 'on'); % filelist
uimenu(hGUI, 'Parent', hFile, 'Label', 'Previous...', 'Callback', 'Spiky.main.OpenFile([],-1);'); % previous
uimenu(hGUI, 'Parent', hFile, 'Label', '&Next...', 'Callback', 'Spiky.main.OpenFile([],0);', 'Accelerator', 'N'); % next
hMerge  = uimenu(hGUI, 'Parent', hFile, 'Label', '&Merge');
uimenu(hGUI, 'Parent', hMerge, 'Label', '&Distribute Settings', 'Callback', Spiky.main.DistributeSettings);
uimenu(hGUI, 'Parent', hMerge, 'Label', '&Merge Files', 'Callback', Spiky.main.CreateMergeFile, 'Accelerator', 'M');
uimenu(hGUI, 'Parent', hMerge, 'Label', '&Load Merge File...', 'Callback', Spiky.main.LoadMergeFile);
uimenu(hGUI, 'Parent', hMerge, 'Label', '&Delete All Settings', 'Callback', Spiky.main.RemoveSettings, 'Separator', 'on');
uimenu(hGUI, 'Parent', hFile, 'Label', 'Pa&ge Setup...', 'Callback', 'pagesetupdlg(gcbf)', 'separator', 'on');
uimenu(hGUI, 'Parent', hFile, 'Label', 'Print Pre&view...', 'Callback', 'printpreview(gcbf)');
uimenu(hGUI, 'Parent', hFile, 'Label', '&Print...', 'Callback', 'printdlg(gcbf)', 'Accelerator', 'P');
uimenu(hGUI, 'Parent', hFile, 'Label', '&Invert Colors', 'Callback', Spiky.main.InvertColors, 'Checked', 'off');
uimenu(hGUI, 'Parent', hFile, 'Label', 'Export &Figure...', 'Callback', Spiky.main.FigureExport);
uimenu(hGUI, 'Parent', hFile, 'Label', 'E&xit Spiky', 'Callback', Spiky.main.ExitSpiky, 'separator', 'on', 'Accelerator', 'Q');

% Edit menu
hEdit  = uimenu(hGUI, 'Label', '&Edit');
uimenu(hGUI, 'Parent', hEdit, 'Label', '&DAQ Info', 'Callback', Spiky.main.ShowDAQInfo, 'separator', 'on');

% View menu
hView  = uimenu(hGUI, 'Label', '&View');
uimenu(hGUI, 'Parent', hView, 'Label', 'Show Channels...', 'Callback', Spiky.main.SelectChannels);
uimenu(hGUI, 'Parent', hView, 'Label', 'Show Events', 'Callback', Spiky.main.ToggleStatus, 'Checked', 'off', 'Accelerator', 'E');
uimenu(hGUI, 'Parent', hView, 'Label', 'Zoom In', 'Callback', Spiky.main.ZoomIn, 'separator', 'on');
uimenu(hGUI, 'Parent', hView, 'Label', 'Zoom Out', 'Callback', Spiky.main.ZoomOut);
uimenu(hGUI, 'Parent', hView, 'Label', 'Zoom Range', 'Callback', Spiky.main.ZoomRange, 'Accelerator', 'Z');
uimenu(hGUI, 'Parent', hView, 'Label', 'Zoom Reset', 'Callback', Spiky.main.ZoomReset, 'Accelerator', 'X');
uimenu(hGUI, 'Parent', hView, 'Label', 'Pan', 'Callback', Spiky.main.PanWindow, 'Accelerator', 'L');
uimenu(hGUI, 'Parent', hView, 'Label', 'Zoom Amplitude', 'Callback', Spiky.main.ZoomAmplitude, 'Accelerator', 'Y', 'separator', 'on');
uimenu(hGUI, 'Parent', hView, 'Label', 'Zoom Tight', 'Callback', Spiky.main.ZoomTight, 'Accelerator', 'T');
hAmpUnit = uimenu(hGUI, 'Parent', hView, 'Label', 'Amplitude Unit');
uimenu(hGUI, 'Parent', hAmpUnit, 'Label', 'Volts (V)',  'callback', {Spiky.main.SetAmplitudeUnit, 'V'});
uimenu(hGUI, 'Parent', hAmpUnit, 'Label', 'Millivolts (mV)', 'callback', {Spiky.main.SetAmplitudeUnit, 'mV'});
uimenu(hGUI, 'Parent', hAmpUnit, 'Label', 'Microvolts (µV)', 'callback', {Spiky.main.SetAmplitudeUnit, 'uV'});
uimenu(hGUI, 'Parent', hView, 'Label', 'Normal Time', 'Callback', Spiky.main.ToggleStatus, 'separator', 'on', 'Tag', 'ShowNormalTime', 'checked', 'on');

% List of themes
hThemes = uimenu(hGUI, 'Parent', hView, 'Label', '&Themes', 'separator', 'on');
sDir = [sSpikyPath(1:end-7) 'themes' filesep];
tDir = dir([sDir '*.m']);
for th = 1:length(tDir)
    uimenu(hGUI, 'Parent', hThemes, 'Tag', 'Spiky_Menu_ShowTheme', 'Label', [tDir(th).name(1:end-2)], 'Callback', Spiky.main.SetTheme);
end
uimenu(hGUI, 'Parent', hView, 'Label', '&Refresh', 'Accelerator', 'R', ...
    'Callback', Spiky.main.GUIRefresh);

% Channels menu
hChannels = uimenu(hGUI, 'Label', '&Channels');
uimenu(hGUI, 'Parent', hChannels, 'Label', 'Gains...', 'Callback', Spiky.main.SetGain, 'Accelerator', 'G');
uimenu(hGUI, 'Parent', hChannels, 'Label', 'Descriptions...', 'Callback', Spiky.main.ChannelDescriptions);
uimenu(hGUI, 'Parent', hChannels, 'Label', 'Digitize', 'Callback', Spiky.main.DigitizeChannel, 'separator', 'on');
uimenu(hGUI, 'Parent', hChannels, 'Label', 'Duplicate', 'Callback', Spiky.main.DuplicateChannel);
uimenu(hGUI, 'Parent', hChannels, 'Label', 'Delete', 'Callback', Spiky.main.DeleteChannel);
uimenu(hGUI, 'Parent', hChannels, 'Label', 'Delete All Hidden', 'Callback', Spiky.main.DeleteHiddenChannels);
uimenu(hGUI, 'Parent', hChannels, 'Label', 'Noise Reduction... (B)', 'Callback', Spiky.main.PCACleaning);
uimenu(hGUI, 'Parent', hChannels, 'Label', 'Merge...', 'Callback', Spiky.main.MergeChannels);
uimenu(hGUI, 'Parent', hChannels, 'Label', 'Select Filters...', 'Callback', Spiky.main.SetFilterOptions, 'separator', 'on', 'Accelerator', 'F');
hFilters = uimenu(hGUI, 'Parent', hChannels, 'Label', 'Configure Filters');
uimenu(hGUI, 'Parent', hChannels, 'Label', 'Custom Filter...', 'Callback', Spiky.main.SetChannelCalculator, 'Accelerator', 'C');

% Get list of filters
GetScriptList('spikyfilter_', hGUI, hFilters, 'Spiky.main.ConfigureFilter', 0)

uimenu(hGUI, 'Parent', hChannels, 'Label', 'Export settings...', 'Callback', Spiky.main.ExportChannelSettings, 'separator', 'on');
uimenu(hGUI, 'Parent', hChannels, 'Label', 'Import settings...', 'Callback', Spiky.main.ImportChannelSettings);

% Spikes menu
hWaveforms  = uimenu(hGUI, 'Label', '&Spikes');
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Set Threshold...', 'Callback', Spiky.main.SetSpikeThreshold );
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Auto-Detect Thresholds', 'Callback', Spiky.main.AutoDetectThresholds );
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Remove Thresholds', 'Callback', Spiky.main.RemoveThresholds );
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Run Spike Detection (B)', 'Callback', Spiky.main.DetectSpikes, 'separator', 'on');
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Dejitter Spikes (B)', 'Callback', Spiky.main.DejitterSpikes);
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Remove Outliers (B)', 'Callback', Spiky.main.RemoveOutlierSpikes);
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Add Multitrode...', 'Callback', Spiky.main.AddTetrode, 'separator', 'on');
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Remove Multitrode...', 'Callback', Spiky.main.RemoveTetrode );
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Set Max Jitter...', 'Callback', Spiky.main.SetMaxSpikeJitter );
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Show Waveforms', 'Callback', Spiky.main.ViewWaveforms, 'separator', 'on');
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Channel Statistics', 'Callback', Spiky.main.ChannelStatistics);
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Set Pre-Trigger...', 'Callback', Spiky.main.EditPreTriggerTime, 'separator', 'on');
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Set Post-Trigger...', 'Callback', Spiky.main.EditPostTriggerTime);
uimenu(hGUI, 'Parent', hWaveforms, 'Label', 'Set Spike Deadtime...', 'Callback', Spiky.main.SetDeadtime);

% Sorting menu
hSorting  = uimenu(hGUI, 'Label', '&Sort');
uimenu(hGUI, 'Parent', hSorting, 'Label', 'Sort Channel...', 'Callback', Spiky.main.SortSpikes);
uimenu(hGUI, 'Parent', hSorting, 'Label', 'Sort All Channels', 'Callback', Spiky.main.SortAllChannels);
uimenu(hGUI, 'Parent', hSorting, 'Label', 'QuickSort!', 'Callback', Spiky.main.QuickSort);
uimenu(hGUI, 'Parent', hSorting, 'Label', 'Show Sorted Spikes', 'Callback', Spiky.main.ShowSpikeClusters, 'separator', 'on');
uimenu(hGUI, 'Parent', hSorting, 'Label', 'Cluster Control', 'Callback', Spiky.main.ShowOverlappingSpikeClusters, 'Accelerator', 'W');
uimenu(hGUI, 'Parent', hSorting, 'Label', 'Show Aggregration Tree', 'Callback', Spiky.main.ShowAggregationTree, 'Accelerator', 'A');
uimenu(hGUI, 'Parent', hSorting, 'Label', 'View Principal Components', 'Callback', Spiky.main.ViewWaveformPCs);
uimenu(hGUI, 'Parent', hSorting, 'Label', 'Plot Rasters', 'Callback', Spiky.main.PlotRasters);
uimenu(hGUI, 'Parent', hSorting, 'Label', 'Plot Instantaneous Rate', 'Callback', Spiky.main.PlotInstSpikeRate);
uimenu(hGUI, 'Parent', hSorting, 'Label', 'Delete sorting data...', 'Callback', Spiky.main.DeleteSortingData, 'separator', 'on');

% Analysis menu
hAnalysis = uimenu(hGUI, 'Label', '&Analysis');
Spiky.main.CreateAnalysisMenu(hAnalysis, 'continuous')
Spiky.main.CreateAnalysisMenu(hAnalysis, 'discrete')

uimenu(hGUI, 'Parent', hAnalysis, 'Label', 'Experiment &Variables...', 'Callback', Spiky.main.ExperimentDescriptions, 'separator', 'on', 'Accelerator', 'V');

% Tools menu
hTools  = uimenu(hGUI, 'Label', '&Tools');
hExtensions = uimenu(hGUI, 'Parent', hTools, 'Label', 'Extensions');
uimenu(hGUI, 'Parent', hExtensions, 'Label', '&Run Script... (B)', 'Callback', Spiky.main.RunScript);
uimenu(hGUI, 'Parent', hExtensions, 'Label', 'Run &Batch Script...', 'Callback', Spiky.main.RunBatchScript);
uimenu(hGUI, 'Parent', hExtensions, 'Label', 'Redo &Script', 'Callback', Spiky.main.RedoScript, 'Accelerator', 'Y');
uimenu(hGUI, 'Parent', hExtensions, 'Label', 'Get Script Help...', 'Callback', Spiky.main.GetScriptHelp);

% Get list of extensions/scripts
GetScriptList('spiky_', hGUI, hExtensions, 'Spiky.main.RunScript', 1)

uimenu(hGUI, 'Parent', hTools, 'Label', '&Batch Redo', 'Callback', 'Spiky.main.BatchRedo([],''redo'');', 'Accelerator', 'B');
uimenu(hGUI, 'Parent', hTools, 'Label', '&Autoload New Files...', 'Callback', Spiky.main.AutoloadNewFiles, 'Separator', 'on');
uimenu(hGUI, 'Parent', hTools, 'Label', '&Keyboard Mode...', 'Callback', Spiky.main.KeyboardMode, 'Accelerator', 'K');
uimenu(hGUI, 'Parent', hTools, 'Label', '&Clear Persistent Variables', 'Callback', 'clear functions');

% Window menu
hWindow  = uimenu(hGUI, 'Label', '&Window');
uimenu(hGUI, 'Parent', hWindow, 'Label', 'Close All', 'Callback', Spiky.main.CloseSpikyWindows);
uimenu(hGUI, 'Parent', hWindow, 'Label', 'Cascade', 'Callback', Spiky.main.CascadeSpikyWindows);
uimenu(hGUI, 'Parent', hWindow, 'Label', 'All to Front', 'Callback', Spiky.main.BringSpikyFiguresToFront, 'Accelerator', 'F');

% Help menu
hHelp  = uimenu(hGUI, 'Label', '&Help');
hWebResSub = uimenu(hGUI, 'Parent', hHelp, 'Label', 'Web &Resources');
uimenu(hGUI, 'Parent', hHelp, 'Label', '&Version', 'Callback', Spiky.main.ShowVersion, 'Separator', 'on');
uimenu(hGUI, 'Parent', hHelp, 'Label', '&Licence', 'Callback', Spiky.main.ShowLicense);
uimenu(hGUI, 'Parent', hHelp, 'Label', '&Keyboard Shortcuts', 'Callback', Spiky.main.KeyboardShortcuts);
uimenu(hGUI, 'Parent', hHelp, 'Label', '&Check for Updates', 'Callback', Spiky.main.CheckUpdate);
uimenu(hGUI, 'Parent', hHelp, 'Label', '&About Spiky', 'Callback', Spiky.main.AboutSpiky);

% - Web Resources menu
uimenu(hGUI, 'Parent', hWebResSub , 'Label', '&GitHub Repository', 'Callback', 'web(''https://github.com/pmknutsen/spiky'', ''-browser'')');
uimenu(hGUI, 'Parent', hWebResSub , 'Label', '&Wiki', 'Callback', 'web(''https://github.com/pmknutsen/spiky/wiki'', ''-browser'')');
uimenu(hGUI, 'Parent', hWebResSub , 'Label', '&Bug Tracker', 'Callback', 'web(''https://github.com/pmknutsen/spiky/issues'', ''-browser'')');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GetScriptList(sPrefix, hMenu, hParentItem, sCallbackFcn, bSep)
% Get a list of matching filenames and insert them, with callbacks, into a
% UI menu. Function will search for mathing filenames anywhere below the 
% home directory of spiky
sPath = which('spiky');
[sPath, ~] = fileparts(sPath);

% Get a list of all files and folders in this folder.
csFiles = dir(sPath);
csDirs = {csFiles([csFiles.isdir]).name};

csScripts = {};
for p = 3:length(csDirs) % iterate over all directories in path
    sFullPath = sprintf('%s%s%s', sPath, filesep, csDirs{p});
    if exist(sFullPath, 'dir')
        tDir = dir([sFullPath filesep sPrefix '*.m']);
        if ~isempty(tDir)
            csScripts = [csScripts {tDir.name}];
        end
    end
end

% Create links to scripts in the Tools->Scripts menu
for s = 1:length(csScripts)
    sName = strrep(csScripts{s}((length(sPrefix)+1):end-2), '_', ' ');
    iSpace = strfind(sName, ' ');
    sName([1 iSpace+1]) = upper(sName([1 iSpace+1])); % capitalize all words
    hMenu = uimenu(hMenu, 'Label', sName, 'Parent', hParentItem, 'Callback', [sprintf('%s(''%s'')', sCallbackFcn, csScripts{s})]);
    if s == 1 && bSep
        set(hMenu, 'Separator', 'on')
    end
end
return