function tSig = Hilbert_Transform(FV)
% Analyze a continuous signal via the Hilbert transform.
%
% Usage:
%   S = Hilbert_Transform(FV)
%
% Where S is a structure containing the spectrogram of a selected channel,
% its sampling rate (kHz) and start/end times (s).
%

% TODO:
%   Add low-pass parameter (of raw signal)
%   Remove phase information when signal is not fluctuating
%       for this, add a parameter which is a relative value, e.g. std of
%       hilbert amplitude

global Spiky g_bBatchMode
persistent p_sContCh
tSig = struct([]);

% Select channel
if isempty(p_sContCh) || ~g_bBatchMode
    [p_sContCh, bResult] = Spiky.SelectChannelNumber(FV.csDisplayChannels, 'Select continuous signal', p_sContCh);
    if ~bResult return, end
end

% Get data
vCont = double(FV.tData.(p_sContCh)');
if all(size(vCont) > 1) return; end
nFs = FV.tData.([p_sContCh '_KHz']) * 1000;
vCont = Spiky.ChannelCalculator(vCont, p_sContCh);

% Check that signal is not 2D
if all(size(vCont) > 1)
    warndlg('Cannot compute the Hilbert transform on a matrix. Signal needs to be a vector', 'Spiky');
    return;
end

% Get channel descriptive string
sDescr = Spiky.GetChannelDescription(p_sContCh);
if isempty(sDescr) sDescr = p_sContCh; end

% Get parameters interactively
% We don't collect parameters when function in batch-mode (when known)
persistent p_nSigLoPass p_nSPLoPass p_nAmplCutOff
if isempty(p_nSigLoPass) || ~g_bBatchMode
    if isempty(p_nSigLoPass), p_nSigLoPass = 50; end
    if isempty(p_nSPLoPass), p_nSPLoPass = 2; end
    if isempty(p_nAmplCutOff), p_nAmplCutOff = 0.1; end
    
    cPrompt = {'Signal low-pass (Hz)', 'Setpoint low-pass (Hz)', ...
        'Minimum amplitude cut-off'};
    cAnswer = inputdlg(cPrompt,'Options', 1, ...
        {num2str(p_nSigLoPass), num2str(p_nSPLoPass), num2str(p_nAmplCutOff)});
    if isempty(cAnswer), return, end
    p_nSigLoPass = str2double(cAnswer{1}); % hz
    p_nSPLoPass = str2double(cAnswer{2}); % hz
    p_nAmplCutOff = str2double(cAnswer{3}); % unit of original signal
end

% Low-pass signal
[vCont, ~, ~] = Spiky.FilterChannel(vCont, 1:length(vCont), nFs, p_nSigLoPass, 0, 0, 'none');

% Subtract set-point
[vSetPoint, ~, ~] = Spiky.FilterChannel(vCont, 1:length(vCont), nFs, p_nSPLoPass, 0, 0, 'none');
vSig = vCont - vSetPoint;

% Replace NaNs with 0
vSig(isnan(vSig)) = 0;

% Compute the Hilbert transform
vH = hilbert(vSig);

% Hilbert amplitude
vA = abs(vH);

% Hilbert phase
vP = angle(vH); % unwrap(angle(vH))

% Remove phase values where amplitude is less than cut-off
vP((vA < p_nAmplCutOff)) = NaN;

if 0
    figure
    hAx(1) = subplot(3,1,1);
    grid
    plot(vSig)
    hAx(2) = subplot(3,1,2);
    grid
    plot(vP)
    hAx(3) = subplot(3,1,3);
    grid
    plot(vA)
    linkaxes(hAx, 'x')
end

% Create output structure
sPreFix = [sDescr '_Phase'];
tSig(1).(sPreFix) = vP';
tSig.([sPreFix '_KHz']) = nFs / 1000;
tSig.([sPreFix '_TimeBegin']) = FV.tData.([p_sContCh '_TimeBegin']);
tSig.([sPreFix '_TimeEnd']) = FV.tData.([p_sContCh '_TimeEnd']);
tSig.([sPreFix '_Unit']) = 'rad';


return
