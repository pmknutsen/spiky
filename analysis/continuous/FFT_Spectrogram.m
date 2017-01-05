function tSig = FFT_Spectrogram(FV)
% Spectrogram using short-time Fourier transform
%
% Usage:
%   S = Spectrogram(FV)
%
% Where S is a structure containing the spectrogram of a selected channel,
% its sampling rate (kHz) and start/end times (s). Segments that contain
% NaNs are either filled with mirrors of adjacent data (first) or
% interpolate (second) before the spectrogram is computed (to avoid edge
% effects).
%
%

global Spiky g_bBatchMode
tSig = struct([]);

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
% We don't collect parameters when function in batch-mode (when known)
persistent p_nMinF p_nMaxF p_nWinSize p_nWinOverlap p_nD
if isempty(p_nMinF) || ~g_bBatchMode
    if isempty(p_nMinF), p_nMinF = 0; end
    if isempty(p_nMaxF), p_nMaxF = round((nFs / 2.5) * 100) / 100; end
    if isempty(p_nWinSize), p_nWinSize = 1; end % s
    if isempty(p_nWinOverlap), p_nWinOverlap = 0.5; end % s
    if isempty(p_nD), p_nD = 1; end    
    cPrompt = {'Min frequency (Hz)', 'Max frequency (Hz)', 'Window size (s)', ...
        'Window step (s)', 'Signal derivative' };
    cAnswer = inputdlg(cPrompt,'Options', 1, ...
        {num2str(p_nMinF), num2str(p_nMaxF), num2str(p_nWinSize), num2str(p_nWinOverlap), ...
        num2str(p_nD) });
    if isempty(cAnswer), return, end
    p_nMinF = str2double(cAnswer{1}); % hz
    p_nMaxF = str2double(cAnswer{2}); % hz
    p_nWinSize = str2double(cAnswer{3}); % s
    p_nWinOverlap = str2double(cAnswer{4}); % s
    p_nD = str2double(cAnswer{5});
end

% Derivative
if p_nD > 0
    vCont = diff(vCont, p_nD);
end

% Initialize waitbar
hMsg = waitbar(.2, 'Computing short-time Fourier transform spectrogram...');
centerfig(hMsg, Spiky.main.GetGUIHandle());

% Fill in missing NaN values in data
vCont = intnans(vCont);

%%
% Low-pass filter and decimate signal to match user-defined frequency range
% This speeds up the spectral analysis significantly, as out-of-band
% frequencies are ignored.
waitbar(.4, hMsg)
%nBegin = FV.tData.([p_sContCh '_TimeBegin']); % sampling start, sec
%nFs = FV.tData.([p_sContCh '_KHz']) * 1000; % sampling frequency Hz
%vTime = (nBegin+1/nFs):(1/nFs):(nBegin+length(vCont)/nFs); % absolute time, sec
%[vCont, ~, nFs] = Spiky.main.FilterChannel(vCont, vTime, nFs, (p_nMaxF*5), 0, 0, 'decimate');

%% Compute spectrogram
waitbar(.6, hMsg)
nWinOverlap = round(p_nWinOverlap * nFs);
nWinSize = round(p_nWinSize * nFs);
[~, f, t, S] = spectrogram(vCont, nWinSize, nWinOverlap, p_nMinF:1:p_nMaxF, nFs);
S = 10*log10(abs(S)+eps);

%figure
%imagesc(t,f,S)
%power/frequency (dB/Hz)

% Create output structure
sPreFix = [sDescr '_Spec'];
tSig(1).(sPreFix) = S;
nInt = unique(round(diff(t)*1000)/1000);
tSig.([sPreFix '_KHz']) = (1/nInt(1)/1000);
tSig.([sPreFix '_TimeBegin']) = FV.tData.([p_sContCh '_TimeBegin']) + t(1);
tSig.([sPreFix '_TimeEnd']) = FV.tData.([p_sContCh '_TimeBegin']) + ((size(S, 1)+1) * nInt(1));
tSig.([sPreFix '_Unit']) = 'Hz';
tSig.([sPreFix '_Scale']) = f;

% Create a Properties field to record spectrogram variables
cVars = {'p_nMinF' 'p_nMaxF' 'p_nWinSize' 'p_nWinOverlap' 'p_nD'};
tProps = struct('');
for i = 1:length(cVars)
    tProps(i).Var = cVars{i};
    tProps(i).Descr = cPrompt{i};
    tProps(i).Value = eval(cVars{i});
end
tSig.([sPreFix '_Properties']) = tProps;
close(hMsg);

return

