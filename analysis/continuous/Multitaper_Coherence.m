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
    [p_sContChA, bResult] = Spiky.SelectChannelNumber(FV.csDisplayChannels, 'Select continuous signal A', p_sContChA);
    if ~bResult return, end
end

% Select channel B
if isempty(p_sContChB) || ~g_bBatchMode
    [p_sContChB, bResult] = Spiky.SelectChannelNumber(FV.csDisplayChannels, 'Select continuous signal B', p_sContChB);
    if ~bResult return, end
end

% Get parameters interactively (pre/post times and stimulus delay)
% We don't collect parameters when function is called with outputs or if we
% are in batch-mode (assuming parameters are known, which they should be).
if isempty(p_nMinF) || ~g_bBatchMode
    if isempty(p_nMinF), p_nMinF = 1; end
    if isempty(p_nMaxF), p_nMaxF = 10000; end
    if isempty(p_nTW), p_nTW = 5; end
    if isempty(p_nD), p_nD = 1; end
    
    cPrompt = {'Min frequency (Hz)', 'Max frequency (Hz)', 'Time-bandwidth product (TW; s*Hz)', ...
        'Signal derivative'};
    cAnswer = inputdlg(cPrompt,'Options', 1, ...
        {num2str(p_nMinF), num2str(p_nMaxF), num2str(p_nTW), num2str(p_nD)});
    if isempty(cAnswer), return, end
    p_nMinF = str2num(cAnswer{1}); % hz
    p_nMaxF = str2num(cAnswer{2}); % hz
    p_nTW = str2num(cAnswer{3});
    p_nD = str2num(cAnswer{4});
end

% Fetch data
vContA = double(FV.tData.(p_sContChA)');
if all(size(vContA) > 1) return; end
nFsA = FV.tData.([p_sContChA '_KHz']) * 1000;

vContB = double(FV.tData.(p_sContChB)');
if all(size(vContB) > 1) return; end
nFsB = FV.tData.([p_sContChB '_KHz']) * 1000;

% Decimate signals to match user-defined frequency range
nRA = floor(nFsA / (p_nMaxF*2.5));
vContA = decimate(vContA, nRA);
nFsA = nFsA/nRA;

nRB = floor(nFsB / (p_nMaxF*2.5));
vContB = decimate(vContB, nRB);
nFsB = nFsB/nRB;

% TODO: Check nFsA and nFsB are equal 
%       If not, resample the longest trace down

% mtspec parameters structure
tParams.tapers   = [p_nTW p_nTW*2-1]; % [NW K]
tParams.pad      = 1;
tParams.err      = [2 .1];
tParams.fpass    = [p_nMinF p_nMaxF];
tParams.trialave = 1;
tParams.Fs = nFsA;

%% Segment signal vector
nSegLen = (1 / p_nMinF) * nFsA * 2;
nSegs = floor(length(vContA) / nSegLen);
mS = [];
mContA = [];
mContB = [];
for i = 1:nSegs
    vRange = [(i-2)*nSegLen + nSegLen + 1 (i-1)*nSegLen + nSegLen];
    mContA(:, i) = vContA(vRange(1):vRange(2));
    mContB(:, i) = vContB(vRange(1):vRange(2));
end
if p_nD > 0
    mContA = diff(mContA, p_nD, 1);
    mContB = diff(mContB, p_nD, 1);
end
%%
[C, phi, S12, S1, S2, f, confC, phistd, Cerr] = coherencyc(mContA, mContB, tParams);










hFig = figure;
set(hFig, 'color', [.2 .2 .2], 'name', 'Multitaper Coherence', 'NumberTitle', 'off')
hAx = axes('position', [.1 .1 .74 .8], 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], ...
    'ycolor', [.6 .6 .6], 'xscale', 'li', 'yscale', 'log');
hold on



% Error bar estimate
axes(hAx);
hFill = fill([f fliplr(f)], [Serr(1,:) fliplr(Serr(2,:))], vColHi, 'edgecolor', 'none');
hLine = plot(hAx, f, S, 'color', vColLo);


xlabel('Frequency (Hz)')
ylabel('Coherence')
title('Multitaper Coherence', 'color', 'w', 'interpreter', 'none');

return


