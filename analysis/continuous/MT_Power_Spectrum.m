function MT_Spectrogram(FV)
% Multi taper power spectrum using Chronux library functions
%
% Usage:
%   S = Multitaper_Power_Spectrum(FV)
%
% Where S is a structure containing the spectrogram of a selected channel,
% its sampling rate (kHz) and start/end times (s). Segments that contain
% NaNs are either filled with mirrors of adjacent data (first) or
% interpolate (second) before the spectrogram is computed (to avoid edge
% effects).
%
%

global Spiky g_bBatchMode

persistent p_sContCh
if isempty(p_sContCh) || ~g_bBatchMode
    [p_sContCh, bResult] = Spiky.main.SelectChannelNumber(FV.csDisplayChannels, 'Select continuous signal', p_sContCh);
    if ~bResult; return, end
end
drawnow

% Fetch data
vCont = double(FV.tData.(p_sContCh)');
if all(size(vCont) > 1); return; end
nFs = FV.tData.([p_sContCh '_KHz']) * 1000;
nDur = length(vCont) / nFs; % signal duration, s

% Get channel descriptive string
sDescr = Spiky.main.GetChannelDescription(p_sContCh);
if isempty(sDescr); sDescr = p_sContCh; end

% Get parameters interactively
persistent p_nMinF p_nMaxF p_nTW p_nD p_nP
if isempty(p_nMinF) || ~g_bBatchMode
    if isempty(p_nMinF), p_nMinF = 0.1; end
    if isempty(p_nMaxF), p_nMaxF = 500; end
    if isempty(p_nTW), p_nTW = 5; end
    if isempty(p_nD), p_nD = 1; end
    if isempty(p_nP), p_nP = 0.05; end
    
    cPrompt = {'Min frequency (Hz):', 'Max frequency (Hz):', 'Time-bandwidth product (TW; s*Hz):', ...
        'Signal derivative:', 'Significance level:'};
    cAnswer = inputdlg(cPrompt, 'Options', 1, ...
        {num2str(p_nMinF), num2str(p_nMaxF), num2str(p_nTW), num2str(p_nD), num2str(p_nP)});
    if isempty(cAnswer), return, end
    p_nMinF = str2num(cAnswer{1}); % hz
    p_nMaxF = str2num(cAnswer{2}); % hz
    p_nTW = str2num(cAnswer{3});
    p_nD = str2num(cAnswer{4});
    p_nP = str2num(cAnswer{5});
end

p_nP = max(0.001, min(1, p_nP));
p_nMinF = max(0.01, p_nMinF);

% mtspec parameters structure
tParams.tapers   = [p_nTW p_nTW*2-1]; % [NW K]
tParams.pad      = 2;
tParams.err      = [2 p_nP];
tParams.fpass    = [p_nMinF p_nMaxF];
tParams.trialave = 1;

% Initialize waitbar
hMsg = waitbar(.2, 'Computing multitaper spectrogram...');

%%
% Decimate signal to match user-defined frequency range (speeds up spectral analysis)
if length(vCont) > 10000
    nR = floor(nFs / (p_nMaxF*2.5));
    vCont = decimate(vCont, nR);
    nFs = nFs/nR;
    tParams.Fs = nFs;
end
waitbar(.4, hMsg)

% Split signal into signals (for averaging)
nSegLen = (1 / p_nMinF) * nFs * 2;
nSegLen = min(length(vCont), nSegLen);
nSegs = floor(length(vCont) / nSegLen);
mS = [];
mCont = [];
for i = 1:nSegs
    vRange = [(i-2)*nSegLen + nSegLen + 1 (i-1)*nSegLen + nSegLen];
    mCont(:, i) = vCont(vRange(1):vRange(2));
end
waitbar(.6, hMsg)

% Compute derivatives
if p_nD > 0
    mCont = diff(mCont, p_nD, 1);
end
[S, f, Serr] = mtspectrumc(mCont, tParams);
waitbar(1, hMsg)
close(hMsg);

%% Plot
hFig = figure;
set(hFig, 'name', 'Multitaper Spectra')
Spiky.main.ThemeObject(hFig);
hAx = axes('position', [.1 .125 .85 .8], 'xscale', 'li', 'yscale', 'log');
Spiky.main.ThemeObject(hAx);
hold(hAx, 'on')
axes(hAx);
box(hAx, 'on')
grid(hAx, 'on')

hPw = plot(hAx, f, S);
vCol = get(hPw, 'color');
vColErr = min([1 1 1], vCol + .2);
vColPw = max([0 0 0], vCol - .2);
delete(hPw);
hPw = plot(hAx, f, Serr, f, S);
set(hPw(1:2), 'color',  vColErr)
set(hPw(3), 'color',  vColPw)

xlabel('Frequency (Hz)')
ylabel('Spectral Power')
hTit = title(sprintf('Multitaper Power Spectrum:  %s', p_sContCh), 'interpreter', 'none');
Spiky.main.ThemeObject(hTit);

%%

return
