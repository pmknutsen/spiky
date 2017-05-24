function tSig = Wavelet_Spectrogram(FV)
% Wavelet spectrogram
%
% Usage:
%   S = Wavelet_Spectrogram(FV)
%
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
% We don't collect parameters when function is in batch-mode (when known)
persistent p_nMinF p_nMaxF p_nFSteps p_nD
if isempty(p_nMinF) || ~g_bBatchMode
    if isempty(p_nMinF)
        % Estimate the lowest detectable frequency
        p_nMinF = round((2.5 / nDur) * 1000)  / 1000;
    end
    if isempty(p_nMaxF)
        % Estimate the highest detectable frequency
        p_nMaxF = round((nFs / 2.5) * 100) / 100;
    end
    if isempty(p_nFSteps), p_nFSteps = 30; end
    if isempty(p_nD), p_nD = 1; end
    
    cPrompt = {'Min frequency (Hz)', 'Max frequency (Hz)', 'Frequency steps', 'Signal derivative' };
    cAnswer = inputdlg(cPrompt,'Options', 1, ...
        {num2str(p_nMinF), num2str(p_nMaxF), num2str(p_nFSteps), num2str(p_nD) });
    if isempty(cAnswer), return, end
    p_nMinF = str2double(cAnswer{1}); % hz
    p_nMaxF = str2double(cAnswer{2}); % hz
    p_nFSteps = str2double(cAnswer{3});
    p_nD = str2double(cAnswer{4});
end

% Derivative
if p_nD > 0
    vCont = diff(vCont, p_nD);
end

% Initialize waitbar
hWait = waitbar(0, 'Computing wavelet spectrogram...');
centerfig(hWait, Spiky.main.GetGUIHandle());

% Interpolate NaNs
vCont = intnans(vCont);

% Frequencies in logarithmic space (makes sense if the spectra is wide)
p_nFSteps = 30; % how many bins between p_nMinF and p_nMaxF
vFreq = logspace(log10(p_nMinF), log10(p_nMaxF), p_nFSteps);

% Wavelet time vector (can be extended to -2:2 s)
vWaveTime = -1 : 1 / nFs : 1; % s

nHalfWaveLen = (length(vWaveTime) - 1) / 2; % half wavelet length (to correct edge artifacts)
nConvLen = length(vWaveTime) + length(vCont) - 1; % FFT length (wavelet + data)
nConvPow2 = pow2(nextpow2(nConvLen)); % pow2 is faster and more accurate

% Create wavelet cycles in log space
% Many cycles increase frequency but decrease temporal resolution... The lowest cycles
% should not be lower than 3 (you want to have at least 3 cycles for the calculation)
% e.g. 3 cycles for lowest and 16 for highest frequency.
% This is how many cycles are used to compute the specific frequency spectra.
vWaveCycs = logspace(log10(3), log10(16), p_nFSteps);

% FFT of the data
vFFT = fft(vCont, nConvPow2);

%%
for fi = 1:length(vFreq)
    waitbar(fi/length(vFreq), hWait)

    % Create complex wavelet and get its FFT (multiply gaussian with complex wavelet)
    vWavelet = (pi*vFreq(fi)*sqrt(pi)) ^ -.5 * exp(2 * 1i * pi * vFreq(fi) .* vWaveTime) .* exp(-vWaveTime .^ 2 ./ (2 * ( vWaveCycs(fi) /(2 * pi * vFreq(fi))) ^ 2)) / vFreq(fi);
    vFFTWave = fft(vWavelet, nConvPow2);
    
    % Multiply vFFT with vFFTWave and compute the inverse fft
    vInvDTF = ifft(vFFTWave .* vFFT, nConvPow2);
    
    % Now you have your frequency spectra but the vector is too long. You
    % have to cut the second half of the transformation (after nConvLen)
    % and then the half of the wavelet at the beginning and end.
    vInvDTF = vInvDTF(1:nConvLen);
    vInvDTF = vInvDTF(nHalfWaveLen+1:end-nHalfWaveLen);
    
    % Compute power
    mPower(fi, :) = abs(vInvDTF) .^ 2;
end

%%
close(hWait)

%% Plot time-frequency-power spectra over time for mean
figure
imagesc(1:size(vCont, 2), vFreq, mPower);

% Plot power spectral density for mean (I am not 100% sure if this is the
% right implementation for power spectral density...)
%concat_data = reshape(mPower, 1,numel(mPower));
%transform = fft(concat_data, nFs/2);
%PSD = transform.*conj(transform) / (nFs / 2);
%f = 1000 / (nFs/2) * (0:127); % this is hard coded... change for different frequency axis
%figure;
%plot(f, PSD(1:128)) % this is hard coded... change for different frequency axis
%title('Power spectral density over all channels')
%xlabel('Frequency (Hz)')

%%

return



% Low-pass filter and decimate signal to match user-defined frequency range
% This speeds up the spectral analysis significantly, as out-of-band
% frequencies are ignored.

% Remove values from spectrogram where data was a NaN
waitbar(.8, hWait)
if any(iNaN)
    iSNaN = interp1(t, 1:length(t), find(iNaN)./nFs, 'linear', 'extrap');
    iSNaN = unique(round(iSNaN));
    S(iSNaN, :) = NaN;
end

% Remove columns with only zeros, and the columns immediately adjacent
iRemCols = find(all(S == 0, 2));
iRemCols = [min(iRemCols) - 1; iRemCols; max(iRemCols) + 1];
iRemCols(iRemCols < 1 | iRemCols > size(S, 1)) = [];
S(iRemCols, :) = NaN;

% Insert median value where spectrogram is zero
S(S == 0) = median(S(:));

% Create output structure
sPreFix = [sDescr '_MtSpc'];
tSig(1).(sPreFix) = log10(S)';
nInt = unique(round(diff(t)*1000)/1000);
tSig.([sPreFix '_KHz']) = (1/nInt(1)/1000);
tSig.([sPreFix '_TimeBegin']) = FV.tData.([p_sContCh '_TimeBegin']) + t(1);
tSig.([sPreFix '_TimeEnd']) = FV.tData.([p_sContCh '_TimeBegin']) + ((size(S, 1)+1) * nInt(1));
tSig.([sPreFix '_Unit']) = 'Hz';
tSig.([sPreFix '_Scale']) = f;

% Create a Properties field to record spectrogram variables
cVars = {'p_nMinF' 'p_nMaxF' 'p_nWinSize' 'p_nWinStep' 'p_nD'};
tProps = struct('');
for i = 1:length(cVars)
    tProps(i).Var = cVars{i};
    tProps(i).Descr = cPrompt{i};
    tProps(i).Value = eval(cVars{i});
end
tSig.([sPreFix '_Properties']) = tProps;
close(hWait);

return

