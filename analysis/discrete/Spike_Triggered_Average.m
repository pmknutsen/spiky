function Spike_Triggered_Average(FV)
% Plot spike-triggered averages of a continuous signal
%
% Usage:
%   Spike_Triggered_Average(FV)
%
% Both 1D and 2D inputs are supported. See README for description of the FV
% input structure. Outlier groups (if present) are ignored.
%
% TODO:
%

global Spiky g_bBatchMode

% Select spiking channel
persistent p_sSpikeCh;
if isempty(p_sSpikeCh) || (~g_bBatchMode && nargout == 0)
    [p_sSpikeCh, bResult] = Spiky.main.SelectChannelNumber(fieldnames(FV.tSpikes)', 'Select spiking channel', p_sSpikeCh);
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
    [p_sContCh, bResult] = Spiky.main.SelectChannelNumber(FV.csChannels(vIndx)', 'Select continuous signal', p_sContCh);
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

% Determine if continuous channel is 1D or 2D
if length(find(size(vCont)>1)) == 1
    bIs2D = false;
else bIs2D = true; end

% Initialize figure
hFig = figure;
set(hFig, 'Name', 'Spike Triggered Average', 'NumberTitle', 'off')
Spiky.main.ThemeObject(hFig);

% Get unit IDs
if isfield(FV.tSpikes.(p_sSpikeCh), 'hierarchy')
    vUnits = unique(FV.tSpikes.(p_sSpikeCh).hierarchy.assigns); % unit names
else vUnits = NaN; end

% Remove outlier group from vUnits
vUnits(vUnits == 0) = [];

% Iterate over units
Spiky.main.SpikyWaitbar(0, length(vUnits));
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
    figure(hFig)
    hAx = subplot(1, length(vUnits), u); %axes('position', [(nW+.5/length(vUnits))*(u-1)+.05 .1 nW .8] );
    Spiky.main.ThemeObject(hAx);

    vTime = (-nPreLen:nPostLen) .* (1/(nContFs));
    vMean = zeros(size(vTime));
    vStdErr = zeros(size(vTime));

    hold on
    vCol = FV.mColors(u, :);
    hFill = fill([vTime fliplr(vTime)], [vMean+vStdErr fliplr(vMean-vStdErr)], vCol); % error fill
    
    if bIs2D
        vScale = FV.tData.([p_sContCh '_Scale']);
        hAvg = imagesc(vTime./1000, vScale, zeros(length(vTime), length(vScale)));
        set(hAx, 'ylim', [min(vScale) max(vScale)]);
        hRefLine = plot([0 0], [min(vScale) max(vScale)], 'w--');
    else
        hAvg = plot(vTime./1000, vMean, 'color', vCol);
    end
    set(hAx, 'fontsize', 7, 'xlim', [vTime(1) vTime(end)]./1000)
    
    if u == 1 ylabel('Spikes/s'); end % not sure if this belongs here?
    if isfield(FV.tData, [p_sContCh '_Unit'])
        ylabel(hAx, FV.tData.([p_sContCh '_Unit']));
    end
    
    xlabel('Time (s)')
    box on; grid on
    hTit = title('');
    Spiky.main.ThemeObject(hTit)
    
    % Compute average
    mTrials = [];
    i = 0;
    for s = 1:length(vSpiketimes)
        [~, nMinIndx] = min(abs(vContTime - vSpiketimes(s)));
        vSpan = nMinIndx + (-nPreLen:nPostLen);
        if any(vSpan < 1) || any(vSpan > length(vCont)) continue; end

        if bIs2D % 2D
            mThis = vCont(:,vSpan);
            mTrials(end+1, :, :) = mThis;
        else % 1D
            mTrials = [mTrials; vCont(vSpan)];
        end
        
        i = i + 1;
        if i > 100 % update plot every 500 spikes
            vMean = squeeze(nanmean(mTrials, 1));
            vStdErr = squeeze(nanstd(mTrials, 0, 1)) ./ sqrt(size(mTrials, 1));
            i = 0;
            if ~ishandle(hFig) return; end

            if bIs2D
                % Update average
                set(hAvg, 'cdata', vMean)
            else
                % Update average
                set(hAvg, 'xdata', vTime, 'ydata', vMean)
                % Update standard error
                delete(hFill)
                hFill = fill([vTime fliplr(vTime)], [vMean+vStdErr fliplr(vMean-vStdErr)], vCol); % error fill
                set(hFill, 'edgeColor', 'none', 'faceAlpha', 0.5)
                uistack(hAvg) % bring mean line to front
            end
            % Update title
            set(hTit, 'string', sprintf('%d / %d spikes averaged', s, length(vSpiketimes)))
            drawnow
        end
    end
    vMean = squeeze(nanmean(mTrials, 1));
    vStdErr = squeeze(nanstd(mTrials, 0, 1)) ./ sqrt(size(mTrials, 1));
    if bIs2D
        % Update average
        set(hAvg, 'cdata', vMean)
    else
        set(hAvg, 'xdata', vTime, 'ydata', vMean)
        delete(hFill)
        hFill = fill([vTime fliplr(vTime)], [vMean+vStdErr fliplr(vMean-vStdErr)], vCol); % error fill
        set(hFill, 'edgeColor', 'none', 'faceAlpha', 0.5)
        uistack(hAvg)
    end
    Spiky.main.ThemeObject(hAx);
    
    % Title
    if isnan(vUnits(u)) sID = ' Unit UN-SORTED';
    else sID = sprintf(' Unit %d', vUnits(u)); end
    set(hTit, 'string', sprintf('%s  n=%d', sID, length(vSpiketimes)));
    
    Spiky.main.SpikyWaitbar(u, length(vUnits));
    drawnow
end
drawnow

return
