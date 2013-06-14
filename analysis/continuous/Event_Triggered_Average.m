function varargout = Event_Triggered_Average(FV)
% Outputs: [mean, error, time]
%
%

global Spiky g_bBatchMode;

% Select trigger event
persistent p_sEventCh;
if isempty(p_sEventCh) || (~g_bBatchMode && nargout == 0)
    [p_sEventCh, bResult] = Spiky.SelectChannelNumber(FV.csDigitalChannels', 'Select trigger event', p_sEventCh);
    if ~bResult, return, end
end

% Select continuous channel
vIndx = [];
for ch = 1:length(FV.csChannels)
    if isempty(find(strcmp(FV.csChannels(ch), FV.csDigitalChannels), 1))
        vIndx(end+1) = ch;
    end
end
persistent p_sContCh;
if isempty(p_sContCh) || (~g_bBatchMode && nargout == 0)
    [p_sContCh, bResult] = Spiky.SelectChannelNumber(FV.csChannels(vIndx)', 'Select continuous signal', p_sContCh);
    if ~bResult, return, end
end

% Get event up times
vEventTimes = FV.tData.([p_sEventCh '_Up']); % sec, abs time

% Error handling
% i) Check that we have any events
if isempty(vEventTimes)
    waitfor(warndlg(sprintf('No event triggers were detected.\nCannot display event triggered average.')));
    return
end
% ii)Check that we have data
if ~isfield(FV.tData, [p_sContCh '_KHz'])
    waitfor(warndlg(sprintf('Missing data field: %s. Analysis aborted.', [p_sContCh '_KHz'])));
    return
end

% Get continuous signal
nContFs = FV.tData.([p_sContCh '_KHz']) * 1000; % Hz
vCont = Spiky.ChannelCalculator(FV.tData.(p_sContCh), p_sContCh);
vContBegin = FV.tData.([p_sContCh '_TimeBegin']);

% Get parameters interactively (pre/post times and stimulus delay)
% We don't collect parameters when function is called with outputs or if we
% are in batch-mode (assuming parameters are known, which they should be).
persistent p_nStimDel p_nPre p_nPost p_nDetrend p_nFirstPulse p_nLastPulse p_nDerivative
persistent p_bAbsolute p_bLowPassHz p_bHiPassHz p_bHilbert
if isempty(p_nStimDel) || (~g_bBatchMode && nargout == 0)
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
    p_nLastPulse = str2num(cAnswer{6}); %#ok<*ST2NM>
    p_nDerivative = str2num(cAnswer{7});
    p_bAbsolute = str2num(cAnswer{8});
    p_bHiPassHz = str2num(cAnswer{9});
    p_bLowPassHz = str2num(cAnswer{10});
    p_bHilbert = str2num(cAnswer{11});
end

% Limit last pulse to max number of events
nLastPulse = min([p_nLastPulse length(vEventTimes)]);

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
[vCont, ~, nContFs] = Spiky.FilterChannel(vCont, vTime, nContFs, p_bLowPassHz, p_bHiPassHz, p_bAbsolute, 'none');

% Convert unit to Intensity/s^-d
vCont = vCont * (nContFs^p_nDerivative);

% Compute event triggered average
SpikyWaitbar(0, 20);
nLen = length(vCont);
mAll = [];
mAllBaseline = [];
for o = p_nFirstPulse:nLastPulse
    SpikyWaitbar((o/length(vOnsets))*20, 20);
    nStart = vOnsets(o) - round(p_nPre * nContFs);
    nEnd   = vOnsets(o) + round(p_nPost * nContFs);
    if nStart < 1 || nEnd > nLen, continue; end
    vThis = vCont(nStart:nEnd);
    vThisBaseline = vCont((vOnsets(o) - round(.1 * nContFs)):vOnsets(o)); % 100 ms before pulse
    
    % Detrend
    if p_nDetrend
        vThis = detrend(vThis);
        vThisBaseline = detrend(vThisBaseline);
    end

    % Hilbert transform
    if p_bHilbert
        % Subtract set-point (< 2 Hz)
        vSetPoint = filter_series(double(vThis(:)), nContFs, 2);
        vThis = vThis - vSetPoint';
        
        % Replace NaNs with 0
        vNaNIndx = isnan(vThis);
        vThis(vNaNIndx) = 0;
        
        % Hilbert transform
        vHilb = hilbert(vThis);
        
        % Hilbert magnitude and phase
        vThis = abs(vHilb);
    end
    
    mAll(end+1, :) = vThis;
    mAllBaseline(end+1, :) = vThisBaseline;
end
SpikyWaitbar(20, 20);

% Compute statistics
vMean = nanmean(mAll);
vMedian = nanmedian(mAll);
vErrMean = nanstd(mAll) ./ sqrt(size(mAll, 1));
vStdMean = nanstd(mAll);
vErrMedian = vErrMean * 1.25;
vTime = linspace(-p_nPre, p_nPost, size(mAll, 2)); % sec

% Initialize figure and plot
hFig = figure('color', [.2 .2 .2], 'units', 'pixels', 'name', 'Spiky - Event Triggered Average');
if exist('centerfig') centerfig(hFig, hWin); end
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

hCheck = uicontrol(hFig, 'Style', 'checkbox', 'Position', [1 1 70 20], 'String', 'Mean', ...
    'HorizontalAlignment', 'left', 'backgroundcolor', [.2 .2 .2], 'foregroundcolor', 'w', 'value', 1);
set(hCheck, 'Callback', 'if(get(gcbo,''value'')),V=''on'';else,V=''off'';end;set(findobj(''tag'',''ETAMean''),''visible'',V);')

hCheck = uicontrol(hFig, 'Style', 'checkbox', 'Position', [71 1 70 20], 'String', 'Median', ...
    'HorizontalAlignment', 'left', 'backgroundcolor', [.2 .2 .2], 'foregroundcolor', 'w', 'value', 0);
set(hCheck, 'Callback', 'if(get(gcbo,''value'')),V=''on'';else,V=''off'';end;set(findobj(''tag'',''ETAMedian''),''visible'',V);')

hCheck = uicontrol(hFig, 'Style', 'checkbox', 'Position', [141 1 70 20], 'String', 'Trials', ...
    'HorizontalAlignment', 'left', 'backgroundcolor', [.2 .2 .2], 'foregroundcolor', 'w', 'value', 0);
set(hCheck, 'Callback', 'if(get(gcbo,''value'')),V=''on'';else,V=''off'';end;set(findobj(''tag'',''AllTrials''),''visible'',V);')

% Estimate latency
vMeanPostStim = vMean(vTime >= 0);

% baseline: If p_nPre window is negative (i.e. AFTER stimulus), estimate
% baseline from the 0.1 s pre-stim window by default.
if p_nPre <= 0, nBaseL = mean(mAllBaseline(:));
else, nBaseL = mean(vMean(vTime < 0)); end

% find peak
[nMax, nMaxIndx] = max(vMeanPostStim);

% find first sample where signal > max/2
vIndx = find(vMeanPostStim(1:nMaxIndx) >= nBaseL + (nMax - nBaseL)/2);
if isempty(vIndx), nLatS = NaN;
else
    nLatency = vIndx(1); % index
    [nY, nI] = min(abs(vTime - abs(min([0 p_nPre]))));
    %find( vTime == abs(min([0 p_nPre])) )
    nLatS = vTime(nI + nLatency); %s
end
hH(1) = plot([nLatS nLatS], [nBaseL max(vMean)], 'r--', 'linewidth', 2);
hH(2) = plot([nLatS nLatS], [min(vMean) nBaseL], 'g--', 'linewidth', 2);
legend(flipud(hH), {num2str(max(vMean)-nBaseL) num2str(min(vMean)-nBaseL)}, 'textcolor', 'w')
legend boxoff

% Axes properties
axis tight
set(hAx, 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'xlim', [-p_nPre p_nPost])
xlabel('Time (s)')
if p_nDerivative > 0
    ylabel(sprintf('Intensity/s^-%d', p_nDerivative))
else
    ylabel('Intensity')
end
vIndx = strfind(FV.sLoadedTrial, filesep);
if isempty(vIndx)
    % Figure title
    title(sprintf('%s\nCh=%s  N=%d events  Ampl=%.3f  SD=%.3f  Lat=%.1f ms', FV.sLoadedTrial, ...
        p_sContCh, size(mAll,1), max(vMean)-min(vMean), mean(vStdMean), nLatS*1000), ...
        'interpreter', 'none', 'color', [.7 .7 .7], 'FontWeight', 'bold')
else
    sFolder = FV.sLoadedTrial(vIndx(end-1)+1:vIndx(end)-1);
    sFile = FV.sLoadedTrial(vIndx(end)+1:end);
    % Figure title
    title(sprintf('%s\n%s\nCh=%s  N=%d events  Ampl=%.3f  SD=%.3f  Lat=%.1f ms', sFolder, sFile, ...
        p_sContCh, size(mAll,1), max(vMean)-min(vMean), mean(vStdMean), nLatS*1000), ...
        'interpreter', 'none', 'color', [.7 .7 .7], 'FontWeight', 'bold')
end

plot([vTime(1) vTime(end)], [max(vMean) max(vMean)], 'r:') % max value indicator
plot([vTime(1) vTime(end)], [min(vMean) min(vMean)], 'g:') % min value indicator

if nargout > 0, varargout(1) = {vMean}; end
if nargout > 1, varargout(2) = {vErrMean}; end
if nargout > 2, varargout(3) = {vTime}; end
set(hFig, 'toolbar', 'figure')
if ~g_bBatchMode Spiky.BatchRedo([], 'ShowEventTriggeredAverage'); end

return
