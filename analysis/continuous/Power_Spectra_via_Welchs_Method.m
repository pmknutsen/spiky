function Power_Spectra_via_Welchs_Method(FV)
% Compute power spectral density estimates via Welch's method (pwelch)
%
% Spectra of all analog channels are computed and displayed.
%
% Requirements:
%   Signal Processing Toolbox
%

global Spiky

hFig = figure;
set(hFig, 'color', [.2 .2 .2], 'name', 'Spiky Spectral Densities', 'NumberTitle', 'off')
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
    set(hCheckbox, 'Tag', char(nCh), 'callback', Spiky.ToggleWaveforms);
end
xlabel('Frequency (Hz)');
ylabel('Power (arb.)')

hHeader = header('Power Spectral Densities', 12);
set(hHeader, 'color', 'w', 'interpreter', 'none')
set(hAx, 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'xscale', 'log', 'yscale', 'li', 'xtick', [.1 1 10 100 1000 10000], 'xticklabel', [.1 1 10 100 1000 10000])
axes(hAx); axis tight

return

