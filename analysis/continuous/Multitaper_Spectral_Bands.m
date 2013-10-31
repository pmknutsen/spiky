function tSig = Multitaper_Spectral_Bands(FV)
% Average spectral power within a range a frequencies.
%
% Usage:
%   S = Multitaper_Spectral_Bands(FV)
%
% Requirements:
%   Chronux spectral analysis library (re-distributed with Spiky)
%   Spectrogram must be computed first (using Multitaper_Spectrogram)
%
% Default spectral bands:
%       2 -  4 Hz    Delta
%       6 - 10 Hz    Theta
%       8 - 12 Hz    Alpha
%      12 - 30 Hz    Beta
%      40 - 90 Hz    Gamma
%

global Spiky g_bBatchMode
persistent p_sSpecCh p_mSpecBands

% Select spectrogram
if isempty(p_sSpecCh) || ~g_bBatchMode
    % Select spectrogram channel
    iCh = [];
    for i = 1:length(FV.csDisplayChannels)
        if ~isempty(strfind(FV.csDisplayChannels{i}, '_MtSpc'))
            iCh(end+1) = i;
        end
    end
    if isempty(iCh)
        warndlg('No spectrograms were found. You must compute a spectrogram before using this function.', 'Spiky');
        return
    end
    [p_sSpecCh, bResult] = Spiky.main.SelectChannelNumber(FV.csDisplayChannels(iCh), 'Select spectrogram', p_sSpecCh);
    if ~bResult; return, end
end

% Get data
mSpec = FV.tData.(p_sSpecCh);
nFs = FV.tData.([p_sSpecCh '_KHz']);
nTimeBegin = FV.tData.([p_sSpecCh '_TimeBegin']); % s
nTimeEnd = FV.tData.([p_sSpecCh '_TimeEnd']); % s
vScale = FV.tData.([p_sSpecCh '_Scale']); % Hz

% Get parameters interactively
% We don't collect parameters when function in batch-mode (when known)
if isempty(p_mSpecBands) || ~g_bBatchMode
    if isempty(p_mSpecBands)
        p_mSpecBands = [2 4;6 10;8 12; 12 30;40 90];
    end
    cPrompt = {'Spectral bands [min max]. One band per line.'};
    cAns = inputdlg(cPrompt,'Options', 5, {num2str(p_mSpecBands)});
    if isempty(cAns), return, end
    p_mSpecBands = str2double(cAns{1}); % hz
end

% Iterate over frequency bands
for fb = 1:size(p_mSpecBands, 1)
    iBand = vScale >= p_mSpecBands(fb, 1) & vScale <= p_mSpecBands(fb, 2);
    
    % Get band
    mBand = mSpec(iBand, :);
    
    % Average across band in at each time point
    vBandAvg = nanmean(mBand, 1);
    
    % Create new channel that will be returned as result
    iSep = findstr('_', p_sSpecCh);
    sCh = p_sSpecCh(1:(iSep(end)-1));
    sPreFix = sprintf('%s_%d_%dHz', sCh, p_mSpecBands(fb, 1), p_mSpecBands(fb, 2));
    tSig(1).(sPreFix) = vBandAvg;
    tSig.([sPreFix '_KHz']) = nFs;
    tSig.([sPreFix '_TimeBegin']) = nTimeBegin;
    tSig.([sPreFix '_TimeEnd']) = nTimeEnd;
    tSig.([sPreFix '_Unit']) = 'Power';
end

return
