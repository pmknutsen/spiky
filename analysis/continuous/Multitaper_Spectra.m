function Multitaper_Spectra(FV)
% Compute time-frequency spectra via a multitaper method.
%
% Usage:
%   Multitaper_Spectra(FV)
%
% Spectra of all currently displayed analog channels are computed and then
% plotted with jack-knife error-bar estimates.
%
% Requirements:
%   Chronux spectral analysis library (re-distributed with Spiky)
%
%

global Spiky g_bBatchMode

% Get parameters interactively (pre/post times and stimulus delay)
% We don't collect parameters when function is called with outputs or if we
% are in batch-mode (assuming parameters are known, which they should be).
persistent p_nMinF p_nMaxF p_nTW p_nD
if isempty(p_nMinF) || ~g_bBatchMode
    if isempty(p_nMinF), p_nMinF = 1; end
    if isempty(p_nMaxF), p_nMaxF = 10000; end
    if isempty(p_nTW), p_nTW = 5; end
    if isempty(p_nD), p_nD = 1; end
    
    cPrompt = {'Min frequency (Hz)', 'Max frequency (Hz)', 'Time-bandwidth product (TW; s*Hz)', ...
        'Signal derivative'};
    cAnswer = inputdlg(cPrompt,'Options', 1, ...
        {num2str(p_nMinF), num2str(p_nMaxF), num2str(p_nTW), num2str(p_nD)});
    if isempty(cAnswer), return, end
    p_nMinF = str2num(cAnswer{1}); % hz
    p_nMaxF = str2num(cAnswer{2}); % hz
    p_nTW = str2num(cAnswer{3});
    p_nD = str2num(cAnswer{4});
end

% mtspec parameters structure
tParams.tapers   = [p_nTW p_nTW*2-1]; % [NW K]
tParams.pad      = 2;
tParams.err      = [2 .1];
tParams.fpass    = [p_nMinF p_nMaxF];
tParams.trialave = 1;

hFig = figure;
set(hFig, 'color', [.2 .2 .2], 'name', 'Multitaper Spectra', 'NumberTitle', 'off')
hAx = axes('position', [.1 .1 .74 .8], 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], ...
    'ycolor', [.6 .6 .6], 'xscale', 'li', 'yscale', 'log');
hold on
sCh = FV.csDisplayChannels;

for nCh = 1:length(sCh)
    vCont = double(FV.tData.(FV.csDisplayChannels{nCh}));    
    if all(size(vCont) > 1) continue; end
    if isempty(vCont) continue; end

    vCol = FV.mColors(nCh,:);
    vColLo = vCol - [.3 .3 .3];
    vColLo(vColLo < 0) = 0;
    vColHi = vCol + [.3 .3 .3];
    vColHi(vColHi > 1) = 1;
        
    nFs = FV.tData.([FV.csDisplayChannels{nCh} '_KHz']) * 1000;
    
    % Decimate signal to match user-defined frequency range (speeds up spectral analysis)
    nR = floor(nFs / (p_nMaxF*2.5));
    vCont = decimate(vCont, nR);
    nFs = nFs/nR;
    tParams.Fs = nFs;
    
    nSegLen = (1 / p_nMinF) * nFs * 2;
    %nSegLen = nFs * 3;
    nSegs = floor(length(vCont) / nSegLen);
    mS = [];
    mCont = [];
    for i = 1:nSegs
        vRange = [(i-2)*nSegLen + nSegLen + 1 (i-1)*nSegLen + nSegLen];
        mCont(:, i) = vCont(vRange(1):vRange(2));
    end
    if p_nD > 0 mCont = diff(mCont, p_nD, 1); end
    [S, f, Serr] = mtspectrumc(mCont, tParams);
    
    % Error bar estimate
    axes(hAx);
    hFill = fill([f fliplr(f)], [Serr(1,:) fliplr(Serr(2,:))], vColHi, 'edgecolor', 'none');
    hLine = plot(hAx, f, S, 'color', vColLo);

    % Toggle view check box
    set([hFill hLine], 'tag', char(nCh));
    hCheckbox = uicontrol(hFig, 'Style', 'checkbox', 'units', 'normalized', ...
        'Position', [.85 .85-([nCh-1]*.06) .15 .05], 'String', sCh(nCh), ...
        'HorizontalAlignment', 'left', 'backgroundcolor', [.2 .2 .2], ...
        'value', 1, 'foregroundColor', vCol);
    set(hCheckbox, 'Tag', char(nCh), 'callback', Spiky.ToggleWaveforms);
end

xlabel('Frequency (Hz)')
ylabel('Spectral Power')
title('Multitaper Spectra', 'color', 'w', 'interpreter', 'none');

return

