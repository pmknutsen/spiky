function tSig = Multitaper_Coherence(FV)
% Multi taper coherence of two continuous channels.
%
% Usage:
%   S = Multitaper_Coherence(FV)
%
% Requirements:
%   Chronux spectral analysis library (re-distributed with Spiky)
%
% To-do:
%   Everything!!
%   Allow coherence of channels with different sampling rate;
%       i.e. 

global Spiky g_bBatchMode
persistent p_sContChA p_sContChB p_nMinF p_nMaxF p_nTW p_nD
tSig = struct([]);

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

% Fetch data
vContA = double(FV.tData.(p_sContChA)');
% TODO: Get vTime (for downsampling etc below)
if all(size(vContA) > 1); return; end
nFsA = FV.tData.([p_sContChA '_KHz']) * 1000;
nTimeBeginA = FV.tData.([p_sContChA '_TimeBegin']);

vContB = double(FV.tData.(p_sContChB)');
% TODO: Get vTime (for downsampling etc below)
if all(size(vContB) > 1); return; end
nFsB = FV.tData.([p_sContChB '_KHz']) * 1000;
nTimeBeginB = FV.tData.([p_sContChB '_TimeBegin']);

% Get parameters interactively
% We don't collect parameters when function in batch-mode (when known)
if isempty(p_nMinF) || ~g_bBatchMode
    if isempty(p_nMinF), p_nMinF = 1; end
    if isempty(p_nMaxF), p_nMaxF = min(floor([nFsA nFsB]/2.5)); end
    if isempty(p_nTW), p_nTW = 3; end
    if isempty(p_nD), p_nD = 1; end
    cPrompt = {'Min frequency (Hz)', 'Max frequency (Hz)', ...
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

% Derivative
if p_nD > 0
    vContA = diff(vContA, p_nD);
    vContB = diff(vContB, p_nD);
end

% Check that onset time is the same for both signals
if nTimeBeginA ~= nTimeBeginB
    warndlg('Both channels need to start at the same time.', 'Multitaper Coherence')
    return
end

% Check that sample rate is same for both channels
if nFsA ~= nFsB
    % TODO: resample highest Fs trace to lower resolution HERE
    warndlg('Sampling rate must be same for both channels.', 'Multitaper Coherence')
    return
end

% TODO: Fix traces with NaNs HERE

% Decimate signal to match user-defined frequency range (speeds up spectral analysis)
% Note that nFs after decimation must be an integer as there will otherwise
% be a temporal offset in the output of mtspecgramc (due to internal rounding of nFs)
if 0
    nFsO = nFsA;
    nR_max = floor(nFsA / (p_nMaxF*5)); % scalar must be at least 2.5
    nR = nR_max;
    nFsA = nFsO/nR;
    nFsB = nFsO/nR;
    while 1
        % Check that both nR and nFs are integers
        if  (ceil(nR) == floor(nR)) && (ceil(nFsA) == floor(nFsA))
            break
        else
            nR = nR - 1;
            nFsA = nFsO/nR;
            nFsB = nFsO/nR;
        end
        if nR == 0 break; end
    end
    vContA = decimate(vContA, nR);
    vContB = decimate(vContB, nR);
else
    nFsO = nFsA;
end


%% Segment signal vector (to enable trial averaging)
nSegLen = (1 / (p_nMinF+2)) * nFsA * 2;
nSegs = floor(length(vContA) / nSegLen)
mS = [];
mContA = [];
mContB = [];
if nSegs == 0
    mContA = vContA;
    mContB = vContB;
else
    for i = 1:nSegs
        vRange = [(i-2)*nSegLen + nSegLen + 1 (i-1)*nSegLen + nSegLen];
        mContA(:, i) = vContA(vRange(1):vRange(2));
        mContB(:, i) = vContB(vRange(1):vRange(2));
    end
end

% Initialize coherencyc parameters structure
p_nTW = 5;
tParams.tapers   = [p_nTW (p_nTW*2)-1]; % [NW K]
tParams.pad      = 2;
tParams.err      = [2 .05]; % 0 = no error bars, [1 p] = theoretical error bars, [2 p] = jackknife error
tParams.fpass    = [p_nMinF p_nMaxF];
tParams.trialave = 1;
tParams.Fs       = nFsA;

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

%%


% Error bar estimate
axes(hAx);
hFill = fill([f fliplr(f)], [C'+Cerr(1,:) fliplr(C'-Cerr(2,:))], vColHi, 'edgecolor', 'none');
hLine = plot(hAx, f, C, 'color', vColLo);

Spiky.main.ThemeObject(hTit)

%%

return


