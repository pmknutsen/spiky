function tSig = Multitaper_Spectrogram(FV)
% Multi taper spectrogram using Chronux library functions
%
% Usage:
%   S = Multitaper_Spectrogram(FV)
%
% Where S is a structure containing the spectrogram of a selected channel,
% its sampling rate (kHz) and start/end times (s). Segments that contain
% NaNs are either filled with mirrors of adjacent data (first) or
% interpolate (second) before the spectrogram is computed (to avoid edge
% effects).
%
% To-do:
%   See TODO's below
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
if all(size(vCont) > 1) return; end
nFs = FV.tData.([p_sContCh '_KHz']) * 1000;

% Get channel descriptive string
sDescr = Spiky.main.GetChannelDescription(p_sContCh);
if isempty(sDescr) sDescr = p_sContCh; end

% Get parameters interactively
% We don't collect parameters when function in batch-mode (when known)
persistent p_nMinF p_nMaxF p_nWinSize p_nWinStep p_nTW p_nD p_bRemLines
if isempty(p_nMinF) || ~g_bBatchMode
    if isempty(p_nMinF), p_nMinF = 0; end
    if isempty(p_nMaxF), p_nMaxF = floor(nFs/2.5); end
    if isempty(p_nWinSize), p_nWinSize = 0.5; end % todo: guess better...
    if isempty(p_nWinStep), p_nWinStep = 0.1; end % todo: guess better...
    if isempty(p_nTW), p_nTW = 3; end % todo: guess better...
    if isempty(p_nD), p_nD = 1; end
    if isempty(p_bRemLines), p_bRemLines = 0; end
    
    cPrompt = {'Min frequency (Hz)', 'Max frequency (Hz)', 'Window size (s)', ...
        'Window step (s)', 'Time-bandwidth product (TW; s*Hz)', 'Signal derivative', ...
        'Remove line noise' };
    cAnswer = inputdlg(cPrompt,'Options', 1, ...
        {num2str(p_nMinF), num2str(p_nMaxF), num2str(p_nWinSize), num2str(p_nWinStep), ...
        num2str(p_nTW), num2str(p_nD), num2str(p_bRemLines) });
    if isempty(cAnswer), return, end
    p_nMinF = str2double(cAnswer{1}); % hz
    p_nMaxF = str2double(cAnswer{2}); % hz
    p_nWinSize = str2double(cAnswer{3});
    p_nWinStep = str2double(cAnswer{4});
    p_nTW = str2double(cAnswer{5});
    p_nD = str2double(cAnswer{6});
    p_bRemLines = str2double(cAnswer{7});
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
if ~isempty(iTurnsNaN)
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
end

% Decimate signal to match user-defined frequency range (speeds up spectral analysis)
% Note that nFs after decimation must be an integer as there will otherwise
% be a temporal offset in the output of mtspecgramc (due to internal rounding of nFs)
if 0
    nFsO = nFs;
    nR_max = floor(nFs / (p_nMaxF*2.5));
    nR = nR_max;
    nFs = nFsO/nR;
    while 1
        % Check that both nR and nFs are integers
        if  (ceil(nR) == floor(nR)) && (ceil(nFs) == floor(nFs))
            break
        else
            nR = nR - 1;
            nFs = nFsO/nR;
        end
        if nR == 0 break; end
    end
    vCont = decimate(vCont, nR);
else
    nFsO = nFs;
end

% Ensure window size and step (in samples) are integers, to avoid internal
% rounding in mtspecgramc and temporal offset of spectrogram
nWinSize = (round((p_nWinSize * nFs)/10) * 10) / nFs;
nWinStep = (round((p_nWinStep * nFs)/10) * 10) / nFs;

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
[S, t, f] = mtspecgramc(vCont, [nWinSize nWinStep], tParams);
close(hMsg);

% Remove values from spectrogram where data was a NaN
if any(iNaN)
    iSNaN = interp1(t, 1:length(t), find(iNaN)./nFsO, 'linear');
    iSNaN = unique(round(iSNaN));
    S(iSNaN, :) = NaN;
end

% Detect and remove noise lines
if p_bRemLines
    vMeanPrime = diff(nanmean(S, 1));
    nStd = nanstd(vMeanPrime);
    iI = vMeanPrime < (nanmean(vMeanPrime)-3*nStd) | vMeanPrime > (nanmean(vMeanPrime)+3*nStd);
    plot(vMeanPrime)
    plot(iI)
    iI = find(diff(iI) == 1);
    iIon = iI(1:2:end) - 2;
    iIoff = iI(2:2:end) + 2;
    Sp = S;
    for i = 1:min([length(iIon) length(iIoff)])
        Sp(:, (iIon(i):iIoff(i))+3) = NaN;
    end
    S = Sp;
end


%% Subtract mean from spectrogram to remove lines
% TODO
% In presence of line noise, subtract lines by dividing/subtracting by mean
% power at each frequency (correcting for variance may be good too)
%vMean = median(S, 1);
%mMean = repmat(vMean, size(S, 1), 1);
%mRes = S ./ mMean;
%mRes = mRes - min(min(mRes));
%subplot(2,1,1)
%imagesc(log10(S'))
%%

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

