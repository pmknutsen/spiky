function tSig = Multitaper_Coherence(FV)
% Multi taper coherence between two analog channels using Chronux library functions
%
% Usage:
%   S = Multitaper_Coherence(FV)
%
% Requirements:
%   Chronux spectral analysis library (re-distributed with Spiky)
%
% To-do:
%   Everything!!

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
if all(size(vContA) > 1) return; end
nFsA = FV.tData.([p_sContChA '_KHz']) * 1000;

vContB = double(FV.tData.(p_sContChB)');
% TODO: Get vTime (for downsampling etc below)
if all(size(vContB) > 1) return; end
nFsB = FV.tData.([p_sContChB '_KHz']) * 1000;

% Get parameters interactively
% We don't collect parameters when function in batch-mode (when known)
if isempty(p_nMinF) || ~g_bBatchMode
    if isempty(p_nMinF), p_nMinF = 1; end
    if isempty(p_nMaxF), p_nMaxF = min(floor([nFsA nFsB]/2.5)); end
    if isempty(p_nTW), p_nTW = 3; end % todo: guess better...
    if isempty(p_nD), p_nD = 1; end % todo: guess better...
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


if nFsA ~= nFsB
    % TODO: resample highest Fs trace to lower resolution HERE
    % Check that nFsA now is equal to nFsB
end

% TODO: Fix traces with NaNs HERE

% TODO: Both traces need to start and end at the same time!
% If different, interpolate one trace based on vTime of other trace??


% Decimate signal to match user-defined frequency range (speeds up spectral analysis)
nRA = floor(nFsA / (p_nMaxF*2.5));
vContA = decimate(vContA, nRA);
nFsA = nFsA/nRA;

nRB = floor(nFsB / (p_nMaxF*2.5));
vContB = decimate(vContB, nRB);
nFsB = nFsB/nRB;

%% Segment signal vector (to enable trial averaging)
nSegLen = (1 / p_nMinF) * nFsA * 2;
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
        mContA(:, i) = vContA(vRange(1):vRange(2));
        mContB(:, i) = vContB(vRange(1):vRange(2));
    end
end
if p_nD > 0
    mContA = diff(mContA, p_nD, 1);
    mContB = diff(mContB, p_nD, 1);
end

% Initialize coherencyc parameters structure
tParams.tapers   = [p_nTW p_nTW*2-1]; % [NW K]
tParams.pad      = 1;
tParams.err      = [2 .1];
tParams.fpass    = [p_nMinF p_nMaxF];
tParams.trialave = 1;
tParams.Fs = nFsA;

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

% Error bar estimate
axes(hAx);
hFill = fill([f fliplr(f)], [C'+Cerr(1,:) fliplr(C'-Cerr(2,:))], vColHi, 'edgecolor', 'none');
hLine = plot(hAx, f, C, 'color', vColLo);

xlabel('Frequency (Hz)')
ylabel('Coherence')
hTit = title(sprintf('Multitaper Coherence (n=%d trials)', nSegs), 'color', 'w', 'interpreter', 'none');
Spiky.main.ThemeObject(hTit)

%%

return


