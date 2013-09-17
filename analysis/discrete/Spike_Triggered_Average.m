function Spike_Triggered_Average(FV)
% Plot spike-triggered averages of all units from one selected channel
%

global Spiky g_bBatchMode;

% Select spiking channel
persistent p_sSpikeCh;
if isempty(p_sSpikeCh) || (~g_bBatchMode && nargout == 0)
    [p_sSpikeCh, bResult] = Spiky.SelectChannelNumber(fieldnames(FV.tSpikes)', 'Select spiking channel', p_sSpikeCh);
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
    if ~bResult return, end
end

% Get parameters interactively (pre/post times)
% We don't collect parameters when function is called with outputs or if we
% are in batch-mode (assuming parameters are known, which they should be).
persistent p_nPre p_nPost
if isempty(p_nPre) || ~g_bBatchMode
    if isempty(p_nPre), p_nPre = 0.5; end %s
    if isempty(p_nPost), p_nPost = 0.5; end %s
    cPrompt = {'Pre-spike time (s)', 'Post-spike time (s)'};
    cAnswer = inputdlg(cPrompt,'Options', 1, ...
        {num2str(p_nPre), num2str(p_nPost)});
    if isempty(cAnswer), return, end
    p_nPre = str2num(cAnswer{1}); % s
    p_nPost = str2num(cAnswer{2}); % s
end

% Get continuous data
vCont = FV.tData.(p_sContCh);
nContFs = FV.tData.([p_sContCh '_KHz']); % kHz
nContTimeBegin = FV.tData.([p_sContCh '_TimeBegin']); % s
nContTimeEnd = FV.tData.([p_sContCh '_TimeEnd']); % s
vContTime = linspace(nContTimeBegin, nContTimeEnd, length(vCont));

hFig = figure;
set(hFig, 'color', [.2 .2 .2], 'Name', 'Spike Triggered Average', 'NumberTitle', 'off')

% Get unit IDs
if isfield(FV.tSpikes.(p_sSpikeCh), 'hierarchy')
    vUnits = unique(FV.tSpikes.(p_sSpikeCh).hierarchy.assigns); % unit names
else vUnits = NaN; end

% Iterate over units
Spiky.SpikyWaitbar(0, length(vUnits));
for u = 1:length(vUnits)
    nFs = FV.tSpikes.(p_sSpikeCh).Fs;
    
    % Get spiketimes
    if isnan(vUnits(u))
        vSpiketimes = FV.tSpikes.(p_sSpikeCh).spiketimes(:) ./ nFs; % unsorted unit, sec
    else
        vIndx = FV.tSpikes.(p_sSpikeCh).hierarchy.assigns == vUnits(u);
        vSpiketimes = FV.tSpikes.(p_sSpikeCh).spiketimes(vIndx) ./ nFs; % sorted unit, sec
    end

    % Iterate over spiketimes
    nPreLen = round(p_nPre / (1/(nContFs*1000)));
    nPostLen = round(p_nPost / (1/(nContFs*1000)));

    % Plot spike triggered average
    if ~ishandle(hFig) return; end
    nW = .8/length(vUnits);
    hAx = axes('position', [(nW+.15/length(vUnits))*(u-1)+.05 .1 nW .8] );

    vTime = (-nPreLen:nPostLen) .* (1/(nContFs));
    vMean = zeros(size(vTime));
    vStdErr = zeros(size(vTime));

    hold on
    vCol = FV.mColors(u, :);
    hFill = fill([vTime fliplr(vTime)], [vMean+vStdErr fliplr(vMean-vStdErr)], vCol); % error fill
    
    hAvg = plot(vTime, vMean, 'color', vCol);
    set(hAx, 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'fontsize', 7, ...
        'xlim', [vTime(1) vTime(end)])
    
    if u == 1 ylabel('Spikes/s'); end
    xlabel('Time (ms)')
    box on; grid on
    hTit = title('');
    Spiky.ThemeObject(hTit)
    
    % Compute average
    mTrials = [];
    i = 0;
    for s = 1:length(vSpiketimes)
        [~, nMinIndx] = min(abs(vContTime - vSpiketimes(s)));
        vSpan = nMinIndx + (-nPreLen:nPostLen);
        if any(vSpan < 1) || any(vSpan > length(vCont)) continue; end
        mTrials = [mTrials; vCont(vSpan)];
        i = i + 1;
        if i > 100 % update plot every 500 spikes
            vMean = nanmean(mTrials);
            vStdErr = nanstd(mTrials) ./ sqrt(size(mTrials, 1));
            i = 0;
            if ~ishandle(hFig) return; end
            set(hAvg, 'xdata', vTime, 'ydata', vMean)
            delete(hFill)
            hFill = fill([vTime fliplr(vTime)], [vMean+vStdErr fliplr(vMean-vStdErr)], vCol); % error fill
            set(hFill, 'edgeColor', 'none', 'faceAlpha', 0.5)
            uistack(hAvg)
            set(hTit, 'string', sprintf('%d / %d spikes averaged', s, length(vSpiketimes)))
            drawnow
        end
    end
    vMean = nanmean(mTrials);
    vStdErr = nanstd(mTrials) ./ sqrt(size(mTrials, 1));
    set(hAvg, 'xdata', vTime, 'ydata', vMean)
    delete(hFill)
    hFill = fill([vTime fliplr(vTime)], [vMean+vStdErr fliplr(vMean-vStdErr)], vCol); % error fill
    set(hFill, 'edgeColor', 'none', 'faceAlpha', 0.5)
    uistack(hAvg)
    
    % Title
    if isnan(vUnits(u)) sID = ' Unit UN-SORTED';
    else sID = sprintf(' Unit %d', vUnits(u)); end
    set(hTit, 'string', sprintf('%s  n=%d', sID, length(vSpiketimes)));
    
    Spiky.SpikyWaitbar(u, length(vUnits));
    drawnow
end

return
