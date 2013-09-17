function tSig = Multitaper_Spectrogram(FV)
% Multi taper spectrogram using Chronux library functions
%
% Usage:
%   S = Multitaper_Spectrogram(FV)
%
% Where S is a structure containing the spectrogram of a selected channel,
% its sampling rate (kHz) and start/end times (s). Note that segments
% containing NaNs are either filled with mirrors of adjacent data (first)
% or interpolated (second). Thus, returned spectrograms are without gaps.
%
% Requirements:
%   Chronux spectral analysis library (re-distributed with Spiky)
%
% To-do:
%   Implement code in 'recycling area' at bottom
%   See additional TODO's below

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

% Get channel descriptive string
sDescr = Spiky.GetChannelDescription(p_sContCh);
if isempty(sDescr) sDescr = p_sContCh; end

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
    p_nMinF = str2double(cAnswer{1}); % hz
    p_nMaxF = str2double(cAnswer{2}); % hz
    p_nWinSize = str2double(cAnswer{3});
    p_nWinStep = str2double(cAnswer{4});
    p_nTW = str2double(cAnswer{5});
    p_nD = str2double(cAnswer{6});
end

% Derivative
if p_nD > 0
    vCont = diff(vCont, p_nD);
end

% Interpolate (nearest) indices that contain NaNs
% Step 1 - NaN segments are filled with mirror images adjacent data
% Step 2 - Remaining NaNs (e.g. lone occurences, NaNs at start/end of
%          vector) are interpolated
iNaN = isnan(vCont);
iTurnsNaN = find([0;diff(iNaN)] == 1) - 1; % start indices of nan segments
iNaNDone = find([0;diff(iNaN)] == -1); % end indices of nan segments
iNaNDone(iNaNDone <= iTurnsNaN(1)) = [];
iNaNDone = iNaNDone(1:length(iTurnsNaN));
for i = 1:length(iTurnsNaN)
    nL = ceil((iNaNDone(i)-iTurnsNaN(i)+1)/2);
    if (iTurnsNaN(i)+nL <= length(vCont)) && (iTurnsNaN(i)-nL > 0)
        vCont(iTurnsNaN(i):(iTurnsNaN(i)+nL)) = fliplr(vCont((iTurnsNaN(i)-nL):iTurnsNaN(i))');
    end
    if ((iNaNDone(i)-nL) > 0) && (iNaNDone(i)+nL <= length(vCont))
        vCont((iNaNDone(i)-nL):iNaNDone(i)) = fliplr(vCont(iNaNDone(i):(iNaNDone(i)+nL))');
    end
end
% Interpolate remaining NaNs
iNaNb = isnan(vCont);
if any(iNaNb)
    vCont(iNaNb) = interp1(find(~iNaNb), vCont(~iNaNb), find(iNaNb), 'linear');
end

% Decimate signal to match user-defined frequency range (speeds up spectral analysis)
nR = floor(nFs / (p_nMaxF*2.5));
vCont = decimate(vCont, nR);
nFsO = nFs;
nFs = nFs/nR;

% mtspecgramc parameters structure
% TW - time-bandwidth product (i.e. sec * Hz),
% K  - number of tapers
tParams.tapers = [p_nTW p_nTW*2-1]; % [NW K]
tParams.pad    = 1;
tParams.Fs     = nFs;
tParams.err    = 0; % make sure we don't compute error bars. too slow
tParams.fpass  = [p_nMinF p_nMaxF];

% Compute the multitaper spectrogram
hMsg = waitbar(0.5, 'Computing multitaper spectrogram...');
[S, t, f] = mtspecgramc(vCont, [p_nWinSize p_nWinStep], tParams);
close(hMsg);

% Remove values from spectrogram where data was a NaN
%if any(iNaN)
%    iSNaN = interp1(t, 1:length(t), find(iNaN)./nFsO, 'linear');
%    iSNaN = unique(round(iSNaN));
%    S(iSNaN, :) = NaN;
%end

% Create output structure
sPreFix = [sDescr '_MtSpc'];
tSig(1).(sPreFix) = log10(S)';
nInt = unique(round(diff(t)*1000)/1000);
tSig.([sPreFix '_KHz']) = (1/nInt(1)/1000);
tSig.([sPreFix '_TimeBegin']) = FV.tData.([p_sContCh '_TimeBegin']) + t(1);
tSig.([sPreFix '_TimeEnd']) = FV.tData.([p_sContCh '_TimeBegin']) + ((size(S, 1)+1) * nInt(1));
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

