function Multitaper_Coherence(FV)
% Multi taper coherence of two continuous channels.
%
% Usage:
%   S = Multitaper_Coherence(FV)
%
% Where S is a structure containing the spectrogram of a selected channel,
% its sampling rate (kHz) and start/end times (s). Segments that contain
% NaNs are either filled with mirrors of adjacent data (first) or
% interpolate (second) before the spectrogram is computed (to avoid edge
% effects).
%
% Requirements:
%   Chronux spectral analysis library (re-distributed with Spiky)
%

global Spiky g_bBatchMode
persistent p_sContChA p_sContChB p_nMinF p_nMaxF p_nTW p_nD

% Select channel A
if isempty(p_sContChA) || ~g_bBatchMode
    [p_sContChA, bResult] = Spiky.main.SelectChannelNumber(FV.csDisplayChannels, 'Select continuous signal A', p_sContChA);
    if ~bResult; return, end
end

% Select channel B
if isempty(p_sContChB) || ~g_bBatchMode
    [p_sContChB, bResult] = Spiky.main.SelectChannelNumber(FV.csDisplayChannels, 'Select continuous signal B', p_sContChB);
    if ~bResult; return, end
end

% Fetch signal A
vContA = double(FV.tData.(p_sContChA)');
if all(size(vContA) > 1); return; end
nFsA = FV.tData.([p_sContChA '_KHz']) * 1000; % Hz
nTimeBeginA = FV.tData.([p_sContChA '_TimeBegin']);
nTimeEndA = FV.tData.([p_sContChA '_TimeEnd']);

% Fetch signal B
vContB = double(FV.tData.(p_sContChB)');
if all(size(vContB) > 1); return; end
nFsB = FV.tData.([p_sContChB '_KHz']) * 1000; % Hz
nTimeBeginB = FV.tData.([p_sContChB '_TimeBegin']);
nTimeEndB = FV.tData.([p_sContChB '_TimeEnd']);

% Check that sample rate is same for both channels
if nFsA > nFsB % decimate A
    vContA = [flipud(vContA); vContA; flipud(vContA)]; % mirror to avoid edge effects
    vContA = resample(vContA, round(nFsB), round(nFsA));
    vContA = vContA(round([length(vContA)/3:[length(vContA)] - length(vContA)/3]));
    nFs = nFsB;
elseif nFsA < nFsB % decimate B
    vContB = [flipud(vContB); vContB; flipud(vContB)]; % mirror to avoid edge effects
    vContB = resample(vContB, round(nFsA), round(nFsB));
    vContB = vContB(round([length(vContB)/3:[length(vContB)] - length(vContB)/3]));
    nFs = nFsA;
else
    nFs = nFsA;
end

% Calculate time vectors
vTimeA = Spiky.main.GetTime(nTimeBeginA, nTimeEndA, vContA, nFs);
vTimeB = Spiky.main.GetTime(nTimeBeginB, nTimeEndB, vContB, nFs);

% Force A and B to start at same time and have same length
% NOTE: NaN's remain NaN's and may produce more NaN's
nTimeBegin = max([nTimeBeginA nTimeBeginB]);
nTimeEnd = max([nTimeEndA nTimeEndB]);
vTime = nTimeBegin:(1/nFs):nTimeEnd;
vContA = interp1(vTimeA, vContA, vTime, 'linear');
vContB = interp1(vTimeB, vContB, vTime, 'linear');

% Get parameters interactively
% We don't collect parameters when function in batch-mode (when known)
if isempty(p_nMinF) || ~g_bBatchMode
    if isempty(p_nMinF), p_nMinF = max([round(((1/nFs)*2.5)*100)/100 0.01]); end
    if isempty(p_nMaxF), p_nMaxF = round((nFs/2.5)*100)/100; end
    if isempty(p_nTW), p_nTW = 3; end
    if isempty(p_nD), p_nD = 1; end
    cPrompt = {'Minimum frequency (Hz)', 'Maximum frequency (Hz)', ...
        'Time-bandwidth product (TW; s*Hz)', 'Signal derivative'};
    cAnswer = inputdlg(cPrompt,'Options', 1, ...
        {num2str(p_nMinF), num2str(p_nMaxF), ...
        num2str(p_nTW), num2str(p_nD)});
    if isempty(cAnswer), return, end
    p_nMinF = str2double(cAnswer{1}); % hz
    p_nMaxF = str2double(cAnswer{2}); % hz
    p_nTW = str2double(cAnswer{3});
    p_nD = str2double(cAnswer{4});
end

% Check parameters
if p_nMinF == 0
    warndlg('Minimum frequency must be greater than zero.')
    return;
end

% Derivative (whitening)
if p_nD > 0
    vContA = diff(vContA, p_nD);
    vContB = diff(vContB, p_nD);
end

% Remove NaN's by local mirroring
vContA = intnans(vContA);
vContB = intnans(vContB);

% Decimate signal to match user-defined frequency range (speeds up spectral analysis)
% Note that nFs after decimation must be an integer as there will otherwise
% be a temporal offset in the output of mtspecgramc (due to internal rounding of nFs)
if 1
    nFsO = nFs;
    nR_max = floor(nFs / (p_nMaxF*5));
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
        if nR == 0, break; end
    end
    if nR > 0
        vContA = decimate(vContA, nR);
        vContB = decimate(vContB, nR);
    else
        nFs = nFsO;
    end
end

% Segment signal vector (to enable trial averaging)
nSegLen = (1 / (p_nMinF)) * nFs * 2;
nSegs = floor(length(vContA) / nSegLen);
mS = [];
mContA = [];
mContB = [];
if nSegs == 0
    mContA = vContA;
    mContB = vContB;
else
    for i = 1:nSegs
        vRange = [(i-2)*nSegLen + nSegLen + 1 (i-1)*nSegLen + nSegLen];
        vRange = round(vRange(1):vRange(2));
        mContA(:, i) = vContA(vRange);
        mContB(:, i) = vContB(vRange);
    end
end

% Initialize coherencyc parameters structure
tParams.tapers   = [p_nTW (p_nTW*2)-1]; % [NW K]
tParams.pad      = 2;
tParams.err      = [2 .05]; % 0 = no error bars, [1 p] = theoretical error bars, [2 p] = jackknife error
tParams.fpass    = [p_nMinF p_nMaxF];
tParams.trialave = 1;
tParams.Fs       = nFs;

[C, phi, S12, S1, S2, f, confC, phistd, Cerr] = coherencyc(mContA, mContB, tParams);

% Plot results
vCol = [0 0 1];
vColLo = vCol - [.3 .3 .3];
vColLo(vColLo < 0) = 0;
vColHi = vCol + [.4 .4 .4];
vColHi(vColHi > 1) = 1;

hFig = figure;
set(hFig, 'name', 'Multitaper Coherence')
hAx = axes('position', [.1 .1 .74 .8], 'xscale', 'li', 'yscale', 'li');
Spiky.main.ThemeObject([hFig hAx])
hold(hAx, 'on');

hAx(1) = subplot(2, 1, 1);
hLin(1) = line([f(1) f(end)], confC*[1 1], 'lineStyle', '--'); % confidence level
hold on
hLin(2) = plot(f, C); % hz vs coherence
ylabel('|Coherence|')
xlabel('Frequency (Hz)')
plot(f, Cerr, 'r')

hAx(2) = subplot(2, 1, 2);
hLin(3) = plot(f, (phi));
ylabel('Phase')
xlabel('Frequency (Hz)')

hTit = title(hAx(1), sprintf('Multitaper Coherence\n%s vs %s, n=%d trials, %.0f%% confidence limit', ...
    p_sContChA, p_sContChB, nSegs, tParams.err(2).*100), 'color', 'w', 'interpreter', 'none');
Spiky.main.ThemeObject([hAx hLin hTit])

return



function vCont = intnans(vCont)
% Interpolate (nearest) indices that contain NaNs
% Step 1 - NaN segments are filled with mirror images adjacent data
% Step 2 - Remaining NaNs (e.g. lone occurences, NaNs at start/end of
%          vector) are interpolated
iNaN = isnan(vCont(:));
iTurnsNaN = find([0;diff(iNaN)] == 1) - 1; % start indices of nan segments
if ~isempty(iTurnsNaN)
    iNaNDone = find([0;diff(iNaN)] == -1); % end indices of nan segments
    iNaNDone(iNaNDone <= iTurnsNaN(1)) = [];
    iTurnsNaN(iTurnsNaN > iNaNDone(end)) = [];
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
return


