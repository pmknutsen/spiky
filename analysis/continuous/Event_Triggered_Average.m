function Event_Triggered_Average(FV)
% Outputs: [mean, error, time]
%
% Plot an event triggered average of a selected channel.
%
% TODO
%   Allow to filter out 'noisy' trials, eg. to those with a high variance
%   (use a relative measure)
%
%   Allow to use a filter (e.g. to exclude or include certain triggers)
%

global Spiky g_bBatchMode
persistent p_iSelected p_csContCh

% Select continuous channels (1D or 2D)
% Check if function was called from a context menu. If so, use the Tag to
% detect the channel
%%
csChannels = FV.csChannels;
csDescr = csChannels;

if isfield(FV, 'tChannelDescriptions')
    csChans = {FV.tChannelDescriptions(:).sChannel};
    for i = 1:length(csChannels)
        iCh = find(strcmp(csChannels{i}, csChans));
        if ~isempty(iCh)
            csDescr{i} = FV.tChannelDescriptions(iCh).sDescription;
        end
    end
end
[csChannels, iOrder] = sort(csChannels);
csDescr = csDescr(iOrder);
p_iSelected = Spiky.main.SelectChannelsUI(csChannels, csDescr, p_iSelected);
if ~any(p_iSelected)
    return
end
p_csContCh = csChannels(p_iSelected);
%%

% Select trigger event
persistent p_sEventCh;
if isempty(p_sEventCh) || (~g_bBatchMode && nargout == 0)
    [p_sEventCh, bResult] = Spiky.main.SelectChannelNumber(FV.csDigitalChannels', 'Select trigger event', p_sEventCh);
    if ~bResult, return, end
end

% Select filter event
persistent p_sFilterCh;
if isempty(p_sFilterCh) || (~g_bBatchMode && nargout == 0)
    [p_sFilterCh, bResult] = Spiky.main.SelectChannelNumber(['None'; FV.csDigitalChannels'], 'Select filter event', p_sFilterCh);
    if ~bResult, return, end
end

% Get event up times
vEventTimes = FV.tData.([p_sEventCh '_Up']); % sec, abs time

% Check that we have events
if isempty(vEventTimes)
    waitfor(warndlg(sprintf('No event triggers were detected.\nCannot display event triggered average.')));
    return
end

% Get parameters interactively (pre/post times and stimulus delay)
% We don't collect parameters when function is called with outputs or if we
% are in batch-mode (assuming parameters are known, which they should be).
persistent p_nStimDel p_nPre p_nPost p_nDetrend p_nFirstPulse p_nLastPulse p_nDerivative
persistent p_bAbsolute p_bHilbert p_nFilterTime
if isempty(p_nStimDel) || (~g_bBatchMode && nargout == 0)
    if isempty(p_nStimDel), p_nStimDel = 0; end
    if isempty(p_nPre), p_nPre = 1; end
    if isempty(p_nPost), p_nPost = 2; end
    if isempty(p_nDetrend), p_nDetrend = 0; end
    if isempty(p_nFirstPulse), p_nFirstPulse = 1; end
    if isempty(p_nLastPulse), p_nLastPulse = length(vEventTimes); end
    if isempty(p_nDerivative), p_nDerivative = 0; end
    if isempty(p_bAbsolute), p_bAbsolute = 0; end
    if isempty(p_bHilbert), p_bHilbert = 0; end
    if isempty(p_nFilterTime), p_nFilterTime = 0; end
    cPrompt = {'Pre-event duration (s)','Post-event duration (s)', 'Detrend (1=yes, 0=no)', ...
        'Stimulus delay (ms)', 'First pulse', sprintf('Last pulse (max %d)', length(vEventTimes)), ...
        'Derivative', 'Use absolute values (1=yes, 0=no)', ...
        'Average Hilbert amplitude (1=yes, 0=no; 1D only)', ...
        'Exclude if filter is ON at time relative to trigger (s)'};
    cAnswer = inputdlg(cPrompt,'Averaging options', 1, ...
        {num2str(p_nPre), num2str(p_nPost), num2str(p_nDetrend), num2str(p_nStimDel), ...
        num2str(p_nFirstPulse), num2str(min([p_nLastPulse length(vEventTimes)])), num2str(p_nDerivative), ...
        num2str(p_bAbsolute), num2str(p_bHilbert), ...
        num2str(p_nFilterTime)});
    if isempty(cAnswer), return, end
    p_nPre  = str2num(cAnswer{1}); % sec
    p_nPost  = str2num(cAnswer{2}); % sec
    p_nDetrend  = str2num(cAnswer{3});
    p_nStimDel = str2num(cAnswer{4}); % ms
    p_nFirstPulse = str2num(cAnswer{5});
    p_nLastPulse = str2num(cAnswer{6});
    p_nDerivative = str2num(cAnswer{7});
    p_bAbsolute = str2num(cAnswer{8});
    p_bHilbert = str2num(cAnswer{9});
    p_nFilterTime = str2num(cAnswer{10});
end

% Initialize figure and plot
hFig = figure('units', 'pixels', 'name', 'Spiky - Event Triggered Average', 'visible', 'off');
centerfig(hFig, Spiky.main.GetGUIHandle())
set(hFig, 'visible', 'on')
Spiky.main.ThemeObject(hFig);

%% Iterate over channels
for ch = 1:length(p_csContCh)
    sContCh = p_csContCh{ch};
    
    % Check that we have data
    if ~isfield(FV.tData, [sContCh '_KHz'])
        waitfor(warndlg(sprintf('Missing data field: %s. Analysis aborted.', [sContCh '_KHz'])));
        return
    end

    % Get continuous, filtered signal
    vCont = Spiky.main.GetFilteredChannel(sContCh, FV.tData.(sContCh));
    
    nContFs = FV.tData.([sContCh '_KHz']) * 1000; % Hz
    vContBegin = FV.tData.([sContCh '_TimeBegin']);
    
    % Check if signal is 1D or 2D
    if all(size(vCont) > 1); bIs2D = true;
    else bIs2D = false; end
    
    % Remove event times where the filter event was UP at the selected relative
    % time (s)
    if ~strcmp('None', p_sFilterCh)
        % Get filter UP and DOWN times
        vFilterUPTimes = FV.tData.([p_sFilterCh '_Up']); % s, abs time
        vFilterDOWNTimes = FV.tData.([p_sFilterCh '_Down']); % s, abs time
        
        vRemInd = [];
        for e = 1:length(vEventTimes)
            % Get absolute time
            nAbsTime = vEventTimes(e) + p_nFilterTime;
            
            % Check if filter event was UP at time nAbsTime
            
            % Get last filter UP event before trigger
            nFiltStart = vFilterUPTimes(nAbsTime > vFilterUPTimes);
            if isempty(nFiltStart), continue; end
            nFiltStart = nFiltStart(end);
            
            nFiltEnd = vFilterDOWNTimes(vFilterDOWNTimes > nFiltStart);
            nFiltEnd = nFiltEnd(1); % first filter DOWN event after UP
            
            % Check if trigger occurred before filter DOWN event
            if nAbsTime < nFiltEnd
                vRemInd = [vRemInd e];
            end
        end
        vEventTimes(vRemInd) = [];
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
    if bIs2D
        vTime = linspace(0, (1/nContFs)*length(vCont), size(vCont, 2));
    else
        vTime = linspace(0, (1/nContFs)*length(vCont), length(vCont));
    end
    
    % Convert unit to Intensity/s^-d
    vCont = vCont * (nContFs^p_nDerivative);
    
    % Compute event triggered average
    Spiky.main.SpikyWaitbar(0, 20);
    nLen = length(vCont);
    mAll = [];
    mAllBaseline = [];
    for o = p_nFirstPulse:nLastPulse
        Spiky.main.SpikyWaitbar((o/length(vOnsets))*20, 20);
        nStart = vOnsets(o) - (p_nPre * nContFs);
        nEnd   = vOnsets(o) + (p_nPost * nContFs);
        
        nStartR = vOnsets(o) - round(p_nPre * nContFs);
        nEndR  = vOnsets(o) + round(p_nPost * nContFs);
        nSegLen = nEndR - nStartR;
        
        % Extract interpolated trace from start to end
        if nStart < 1 || nEnd > nLen, continue; end
        if bIs2D
            vThis = vCont(:, nStartR:nEndR);
            [M, N] = size(vThis); % Assuming square matrix.
            [xx, yy] = meshgrid(1:N,1:M); % xx,yy are both outputs of meshgrid (so called plaid matrices).
            sy = 0;
            sx = nStartR - nStart;
            try
                vThis = interp2(xx, yy, vThis, xx+sx, yy+sy, 'spline', nan);
            catch
            end
            
            vThisBaseline = vCont(:, (vOnsets(o) - round(.1 * nContFs)):vOnsets(o)); % 100 ms before pulse
        else
            if nSegLen < 100
                vThis = interp1(1:length(vCont), vCont, linspace(nStart, nEnd, nSegLen), 'spline');
            else
                vThis = vCont(nStartR:nEndR);
            end
            vThisBaseline = vCont((vOnsets(o) - round(.1 * nContFs)):vOnsets(o)); % 100 ms before pulse
        end
        
        % Detrend
        if p_nDetrend
            vThis = detrend(vThis);
            vThisBaseline = detrend(vThisBaseline);
        end
        
        % Hilbert transform
        if p_bHilbert && ~bIs2D
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
        
        if bIs2D
            mAll(end+1, :, :) = vThis;
            mAllBaseline(end+1, :, :) = vThisBaseline;
        else
            mAll(end+1, :) = vThis;
            mAllBaseline(end+1, :) = vThisBaseline;
        end
    end
    Spiky.main.SpikyWaitbar(20, 20);
    
    % Compute statistics
    vMean = squeeze(nanmean(mAll, 1));
    vMedian = squeeze(nanmedian(mAll, 1));
    vErrMean = squeeze(nanstd(mAll, 1)) ./ sqrt(size(mAll, 1));
    vStdMean = squeeze(nanstd(mAll, 1));
    vErrMedian = vErrMean * 1.25;
    if bIs2D
        vTime = linspace(-p_nPre, p_nPost, size(mAll, 3)); % sec
    else
        vTime = linspace(-p_nPre, p_nPost, size(mAll, 2)); % sec
    end
    
    % Plot
    nCols = ceil(sqrt(length(p_csContCh)));
    nRows = floor(sqrt(length(p_csContCh)));
    if ~ishandle(hFig), return; end
    hAx = subplot(nRows, nCols, ch, 'parent', hFig);
    if bIs2D
        hold on
        vY = FV.tData.([sContCh '_Scale']);
        hMedianLine = imagesc(vTime, vY, vMedian);
        hMeanLine = imagesc(vTime, vY, vMean);
        set(hAx, 'ydir', 'normal')
        set(hMedianLine, 'visible', 'off', 'tag', 'ETAMedian');
        set(hMeanLine, 'visible', 'on', 'tag', 'ETAMean');
    else
        [hMeanLine, hMeanErr] = Spiky.main.PlotMeanError(vMean, vErrMean, [0 0.5 1], vTime);
        hold on
        [hMedianLine, hMedianErr] = Spiky.main.PlotMeanError(vMedian, vErrMedian, [0.5 1 0], vTime);
        set([hMedianLine, hMedianErr], 'visible', 'off', 'tag', 'ETAMedian');
        set([hMeanLine, hMeanErr], 'visible', 'on', 'tag', 'ETAMean');
        
        % Plot all individual traces
        hold on
        hTrials = plot(repmat(vTime, size(mAll, 1), 1)', mAll');
        Spiky.main.ThemeObject(hTrials);
        set([hTrials], 'visible', 'off', 'tag', 'AllTrials');
    end
    
    % Axes properties
    Spiky.main.ThemeObject(hAx);
    vYLim(1) = min(vMean) - max(vErrMean);
    vYLim(2) = max(vMean) + max(vErrMean);
    set(hAx, 'xlim', [-p_nPre p_nPost], 'ylim', vYLim * 1.1)
    xlabel('Time (s)')
    if bIs2D sYStr = FV.tData.([sContCh '_Unit']);
    else sYStr = 'Intensity'; end
    if p_nDerivative > 0
        ylabel(sprintf('%s/s^-%d', sYStr, p_nDerivative))
    else
        ylabel(sYStr)
    end
    vIndx = strfind(FV.sLoadedTrial, filesep);
    if isempty(vIndx)
        % Figure title
        if ~bIs2D
            hTit = title(sprintf('%s  N=%d  Amp=%.2f  SD=%.2f', ...
                sContCh, size(mAll,1), max(vMean)-min(vMean), mean(vStdMean)), ...
                'interpreter', 'none' );
        else
            hTit = title(sprintf('%s  N=%d', ...
                sContCh, size(mAll,1)), 'interpreter', 'none' );
        end
    else
        sFolder = FV.sLoadedTrial(vIndx(end-1)+1:vIndx(end)-1);
        sFile = FV.sLoadedTrial(vIndx(end)+1:end);
        % Figure title
        if ~bIs2D
            hTit = title(sprintf('%s  N=%d  Amp=%.3f  SD=%.3f', ...
                sContCh, size(mAll,1), max(vMean)-min(vMean), mean(vStdMean)), ...
                'interpreter', 'none' );
        else
            hTit = title(sprintf('%s  N=%d', ...
                sContCh, size(mAll,1) ), 'interpreter', 'none');
        end
    end
    Spiky.main.ThemeObject(hTit);
    
    if ~bIs2D
        grid on
    else
        plot([0 0], [0 max(vY)], 'w--')
    end
end
%%

% Check buttons
hCheck = uicontrol(hFig, 'Style', 'checkbox', 'Position', [1 1 70 20], 'String', 'Mean', ...
    'HorizontalAlignment', 'left', 'value', 1);
Spiky.main.ThemeObject(hCheck);
set(hCheck, 'Callback', 'if(get(gcbo,''value'')),V=''on'';else,V=''off'';end;set(findobj(''tag'',''ETAMean''),''visible'',V);')

hCheck = uicontrol(hFig, 'Style', 'checkbox', 'Position', [71 1 70 20], 'String', 'Median', ...
    'HorizontalAlignment', 'left', 'value', 0);
Spiky.main.ThemeObject(hCheck);
set(hCheck, 'Callback', 'if(get(gcbo,''value'')),V=''on'';else,V=''off'';end;set(findobj(''tag'',''ETAMedian''),''visible'',V);')

if ~bIs2D
    hCheck = uicontrol(hFig, 'Style', 'checkbox', 'Position', [141 1 70 20], 'String', 'Trials', ...
        'HorizontalAlignment', 'left', 'value', 0);
    Spiky.main.ThemeObject(hCheck);
    set(hCheck, 'Callback', 'if(get(gcbo,''value'')),V=''on'';else,V=''off'';end;set(findobj(''tag'',''AllTrials''),''visible'',V);')
end


if nargout > 0, varargout(1) = {vMean}; end
if nargout > 1, varargout(2) = {vErrMean}; end
if nargout > 2, varargout(3) = {vTime}; end
set(hFig, 'toolbar', 'figure')
if ~g_bBatchMode
    Spiky.main.BatchRedo([], 'ShowEventTriggeredAverage');
end

% Add an export button to the menu
ToolbarExportVectorsBtn(findall(hFig, 'type', 'uitoolbar'), hTrials, 'Export trial vectors')
ToolbarExportVectorsBtn(findall(hFig, 'type', 'uitoolbar'), hMeanLine, 'Export mean vectors')

%%


return
