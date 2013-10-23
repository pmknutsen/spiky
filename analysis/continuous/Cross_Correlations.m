function Cross_Correlations(FV)
% Usage:
%   Cross_Correlations(FV)
%
% Plot cross correlations between all continuous channels.
%
% TODO
%   - Plot CCs only of channels currently displayed in GUI
%   - Plot cross correlations between all continuous data traces; rows vs axes
%   - See additional TODOs below

global Spiky
csChannels = FV.csDisplayChannels;

hFig = figure;
set(hFig, 'name', 'Spiky Continuous Cross Correlations', 'NumberTitle', 'off')
Spiky.ThemeObject(hFig)

% Iterate over channels and generate subplot
for ch1 = 1:length(csChannels)
    vTrace1 = FV.tData.(csChannels{ch1});
    if isempty(vTrace1) || length(find(size(vTrace1)>1)) > 1
        continue
    end
    nFs1 = FV.tData.([csChannels{ch1} '_KHz']) * 1000; % sampling frequency (Hz)
    nBeginTime1 = FV.tData.([csChannels{ch1} '_TimeBegin']); % start of sampling (sec)
    %vTime1 = (nBeginTime1+1/nFs1):(1/nFs1):(nBeginTime1+length(vTrace1)/nFs1); % time vector (sec)

    for ch2 = 1:length(csChannels)
        vTrace2 = FV.tData.(csChannels{ch2});
        if isempty(vTrace2) || length(find(size(vTrace2)>1)) > 1
            continue
        end
        nFs2 = FV.tData.([csChannels{ch2} '_KHz']) * 1000; % sampling frequency (Hz)
        nBeginTime2 = FV.tData.([csChannels{ch2} '_TimeBegin']); % start of sampling (sec)
        %vTime2 = (nBeginTime2+1/nFs2):(1/nFs2):(nBeginTime2+length(vTrace2)/nFs2); % time vector (sec)

        % Abort if the two traces dont have same length, Fs or BeginTime
        % TODO: Allow this, by resampling and interpolating signals to
        % match each other in time...
        if (length(vTrace1) ~= length(vTrace2)) ... % length
                || (nFs1 ~= nFs2) ... % Fs
                || (nBeginTime1 ~= nBeginTime2) % BeginTime
            warndlg('Two traces dont start at same time or have different sampling rates. Cannot cross correlate these.')
            continue
        end

        % axes
        if ~ishandle(hFig) return; end
        nX = .08; nY = .1; nW = .88; nH = .85;
        nNN = length(csChannels);
        hAx = axes('position', [nX+(ch2-1)*(nW/nNN) (1-nY)-ch1*(nH/nNN) nW/nNN nH/nNN], 'Color', [.1 .1 .1]); %#ok<*LAXES>

        % Cross correlate
        nMaxLagSec = .1; % sec
        nMaxLagSamp = nFs1/(1/nMaxLagSec); % samples
        
        % Scaled signals and insert zeros where there are NaNs
        vTrace1Scaled = (vTrace1 - nanmean(vTrace1)) / nanstd(vTrace1);
        vTrace1Scaled(isnan(vTrace1)) = 0;
        vTrace2Scaled = (vTrace2 - nanmean(vTrace2)) / nanstd(vTrace2);
        vTrace2Scaled(isnan(vTrace2)) = 0;

        [vC, ~] = xcorr(vTrace1Scaled, vTrace2Scaled, nMaxLagSamp, 'coeff');

        % X axis values
        vX = -(nMaxLagSamp/nFs1):(1/nFs1):(nMaxLagSamp/nFs1);

        % Plot
        if ch1 == ch2
            plot(vX, vC, 'w', 'linewidth', 1.5)
            set(hAx, 'Color', [.25 .25 .25])
        else
            plot(vX, vC, 'y')
        end
        set(hAx, 'fontsize', 7, 'ylim', [-.5 1] )
        Spiky.ThemeObject(hAx)

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

hChild = get(hFig,'children');
if ~isempty(hChild)
    linkaxes(hChild, 'xy')
    zoom xon
end

hHeader = header('Channel Cross Correlations', 12);
set(hHeader, 'color', 'w', 'interpreter', 'none')

return
