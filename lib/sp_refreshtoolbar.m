function hToolbar = sp_refreshtoolbar()
% Create or refresh toolbar in Spiky GUI
% 
% Usage:
%   sp_refreshtoolbar()
%

global Spiky;

sSpikyPath = which('spiky');
sPath = [sSpikyPath(1:end-7) 'icons/'];

% Remove existing toolbar
hToolbar = findobj(Spiky.main.GetGUIHandle(), 'type', 'uitoolbar');
if ~isempty(hToolbar)
    delete(hToolbar)
end

hToolbar = uitoolbar('Parent', Spiky.main.GetGUIHandle(), 'Tag', 'Spiky_Waitbar');

mCData = im2double(imread([sPath 'tool_open.png'])); mCData(mCData == 0) = NaN; % open
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_Open', 'TooltipString', 'Open file', 'ClickedCallback', Spiky.main.OpenFile);
mCData = im2double(imread([sPath 'tool_save.png'])); mCData(mCData == 0) = NaN; % save
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_Save', 'TooltipString', 'Save', 'ClickedCallback', Spiky.main.SaveResults);
mCData = im2double(imread([sPath 'tool_zoom_in.png'])); mCData(mCData == 0) = NaN; % zoom in
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_ZoomOut', 'TooltipString', 'Zoom in', 'ClickedCallback', Spiky.main.ZoomIn, 'separator', 'on');
mCData = im2double(imread([sPath 'tool_zoom_out.png'])); mCData(mCData == 0) = NaN; % zoom out
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_ZoomOut', 'TooltipString', 'Zoom out', 'ClickedCallback', Spiky.main.ZoomOut);
mCData = im2double(imread([sPath 'tool_zoom_reset.png'])); mCData(mCData == 1) = NaN; % zoom reset
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_ZoomReset', 'TooltipString', 'Zoom reset', 'ClickedCallback', Spiky.main.ZoomReset);
mCData = im2double(imread([sPath 'tool_hand.png'])); mCData(mCData == 0) = NaN; % pan tool
uitoggletool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_Pan', 'TooltipString', 'Pan', 'ClickedCallback', Spiky.main.PanWindow);
[mCData, mCM] = imread([sPath 'left.gif']); mCData = ind2rgb(mCData, mCM); mCData(mCData == 1) = NaN; % pan left
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_PanLeft', 'TooltipString', 'Pan left', 'ClickedCallback', Spiky.main.PanLeft);
[mCData, mCM] = imread([sPath 'right.gif']); mCData = ind2rgb(mCData, mCM); mCData(mCData == 1) = NaN; % pan right
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_PanRight', 'TooltipString', 'Pan right', 'ClickedCallback', Spiky.main.PanRight);
tZoomX = load([sPath 'zoomx.mat']); % zoom X
uipushtool('Parent', hToolbar, 'cdata', tZoomX.cdata, 'Tag', 'Spiky_WaitbarAction_ZoomX', 'TooltipString', 'Zoom X-scale', 'ClickedCallback', Spiky.main.ZoomRange);
tZoomY = load([sPath 'zoomy.mat']); % zoom Y
uipushtool('Parent', hToolbar, 'cdata', tZoomY.cdata, 'Tag', 'Spiky_WaitbarAction_ZoomY', 'TooltipString', 'Zoom Y-scale', 'ClickedCallback', Spiky.main.ZoomAmplitude);
mCData = im2double(imread([sPath 'tool_text.png'])); mCData(mCData == 0) = NaN; % set threshold
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_Threshold', 'TooltipString', 'Set threshold', 'ClickedCallback', Spiky.main.SetSpikeThreshold, 'separator', 'on');
[mCData, mCM] = imread([sPath 'tools_table.gif']); mCData = ind2rgb(mCData, mCM); mCData(mCData == 1) = NaN; % experiment variables
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_ExperimentVariables', 'TooltipString', 'Edit experiment variables', 'ClickedCallback', Spiky.main.ExperimentDescriptions);
mCData = im2double(imread([sPath 'tool_waveforms.png'])); mCData(mCData == 0) = NaN; % view waveforms
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_Waveforms', 'TooltipString', 'View sorted waveforms', 'ClickedCallback', Spiky.main.ShowOverlappingSpikeClusters, 'separator', 'on');
mCData = im2double(imread([sPath 'tool_2waveforms.png'])); mCData(mCData == 0) = NaN; % view waveforms
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_Waveforms', 'TooltipString', 'View sorted waveforms', 'ClickedCallback', Spiky.main.ShowSpikeClusters);
mCData = im2double(imread([sPath 'tool_psth.png'])); mCData(mCData == 0) = NaN; % plot psth
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_PSTH', 'TooltipString', 'Plot PSTH', 'ClickedCallback', Spiky.analysis.discrete.Peristimulus_Time_Histograms);
mCData = im2double(imread([sPath 'tool_pc.png'])); mCData(mCData == 0) = NaN; % principal components
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_PC', 'TooltipString', 'Show principal components', 'ClickedCallback', Spiky.main.ViewWaveformPCs);
mCData = im2double(imread([sPath 'tool_batch.png'])); mCData(mCData == 1) = NaN; % run batch job
uipushtool('Parent', hToolbar, 'cdata', mCData, 'Tag', 'Spiky_WaitbarAction_BatchRedo', ...
    'TooltipString', 'No command available to batch', 'ClickedCallback', 'global Spiky; Spiky.main.BatchRedo('''',''redo'');', 'enable', 'off');

return