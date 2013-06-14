function tSig = Multitaper_Spectrogram(FV)
% Multi taper spectrogram using Chronux library functions
%
% Usage:
%   S = Multitaper_Spectrogram(FV)
%
% Where S is a structure containing the spectrogram of a selected channel,
% its sampling rate (kHz) and start/end times (s).
%
% Requirements:
%   Chronux spectral analysis library (re-distributed with Spiky)
%
% To-do:
%   Collect parameters from user:

global Spiky g_bBatchMode
tSig = struct([]);

persistent p_sContCh
if isempty(p_sContCh) || ~g_bBatchMode
    [p_sContCh, bResult] = Spiky.SelectChannelNumber(FV.csDisplayChannels, 'Select continuous signal', p_sContCh);
    if ~bResult return, end
end
drawnow

% Fetch data
vCont = double(FV.tData.(p_sContCh)');
if all(size(vCont) > 1) return; end
nFs = FV.tData.([p_sContCh '_KHz']) * 1000;

% Get parameters interactively (pre/post times and stimulus delay)
% We don't collect parameters when function is called with outputs or if we
% are in batch-mode (assuming parameters are known, which they should be).
persistent p_nMinF p_nMaxF p_nWinSize p_nWinStep p_nTW p_nD
if isempty(p_nMinF) || ~g_bBatchMode
    if isempty(p_nMinF), p_nMinF = 1; end
    if isempty(p_nMaxF), p_nMaxF = floor(nFs/2.5); end
    if isempty(p_nWinSize), p_nWinSize = 0.5; end % todo: guess better...
    if isempty(p_nWinStep), p_nWinStep = 0.1; end % todo: guess better...
    if isempty(p_nTW), p_nTW = 3; end % todo: guess better...
    if isempty(p_nD), p_nD = 1; end % todo: guess better...
    
    cPrompt = {'Min frequency (Hz)', 'Max frequency (Hz)', 'Window size (s)', 'Window step (s)', ...
        'Time-bandwidth product (TW; s*Hz)', 'Signal derivative'};
    cAnswer = inputdlg(cPrompt,'Options', 1, ...
        {num2str(p_nMinF), num2str(p_nMaxF), num2str(p_nWinSize), num2str(p_nWinStep), ...
        num2str(p_nTW), num2str(p_nD)});
    if isempty(cAnswer), return, end
    p_nMinF = str2num(cAnswer{1}); % hz
    p_nMaxF = str2num(cAnswer{2}); % hz
    p_nWinSize = str2num(cAnswer{3});
    p_nWinStep = str2num(cAnswer{4});
    p_nTW = str2num(cAnswer{5});
    p_nD = str2num(cAnswer{6});
end

% Decimate signal to match user-defined frequency range (speeds up spectral analysis)
nR = floor(nFs / (p_nMaxF*2.5));
vCont = decimate(vCont, nR);
nFs = nFs/nR;

% mtspecgramc parameters structure
% TW - time-bandwidth product (i.e. sec * Hz),
% K  - number of tapers
tParams.tapers = [p_nTW p_nTW*2-1]; % [NW K]
tParams.pad    = 2;
tParams.Fs     = nFs;
tParams.err    = 0; % make sure we don't compute error bars. too slow
tParams.fpass  = [p_nMinF p_nMaxF];

% Compute the multitaper spectrogram
hMsg = waitbar(0.5, 'Computing multitaper spectrogram...');
if p_nD > 0 vCont = diff(vCont, p_nD); end
[S, ~, f] = mtspecgramc(vCont, [p_nWinSize p_nWinStep], tParams);
close(hMsg);

% Create output structure
sPreFix = [p_sContCh '_MtSpc'];
tSig(1).(sPreFix) = log10(S)';
tSig.([sPreFix '_KHz']) = 1/p_nWinStep/1000; % FIX
tSig.([sPreFix '_TimeBegin']) = FV.tData.([p_sContCh '_TimeBegin']) + p_nWinSize / 2;
tSig.([sPreFix '_TimeEnd']) = FV.tData.([p_sContCh '_TimeBegin']) + size(S, 1) * p_nWinStep;
tSig.([sPreFix '_Unit']) = 'Hz';
tSig.([sPreFix '_Scale']) = f;

return


% Some other useful processing we'll recycle later....
% Gamma vs time
tFreqRange(1).hz = [2 4];
tFreqRange(1).name = 'delta (2-4 Hz)';
tFreqRange(2).hz = [6 10];
tFreqRange(2).name = 'theta (6-10 Hz)';
tFreqRange(3).hz = [8 12];
tFreqRange(3).name = 'alpha (8-12 Hz)';
tFreqRange(4).hz = [12 30];
tFreqRange(4).name = 'beta (12-30 Hz)';
tFreqRange(5).hz = [40 90];
tFreqRange(5).name = 'gamma (40-90 Hz)';

for k = 1:length(tFreqRange)
    vFreqRange = tFreqRange(k).hz;
    [~, vFreqRangeI(1)] = min(abs(f-vFreqRange(1)));
    [~, vFreqRangeI(2)] = min(abs(f-vFreqRange(2)));
    vGamma = [];
    for i = 1:size(S, 1)
        mGamma(i,k) = sum(S(i, vFreqRangeI(1):vFreqRangeI(2)));
    end
end

figure
hAx(1) = subplot(length(tFreqRange)+1,1,1);
plot(linspace(0, (length(vCont)/nFs), length(vCont)), vCont)
axis tight; grid on
for k = 1:length(tFreqRange)
    hAx(end+1) = subplot(length(tFreqRange)+1,1, k+1);
    plot(t, log10(mGamma(:, k)))
    legend(tFreqRange(k).name)
    legend boxoff
    axis tight; grid on
end
linkaxes(hAx, 'x')
zoom xon

