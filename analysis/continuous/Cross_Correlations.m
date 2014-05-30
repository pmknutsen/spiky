function Cross_Correlations(FV)
% Usage:
%   Cross_Correlations(FV)
%
% Plot cross correlations of two selected continuous channels
%
%

global Spiky g_bBatchMode
persistent p_sContChA p_sContChB p_nLag p_nD

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

% Detrend
%vContA = detrend(vContA);
%vContB = detrend(vContB);

% Get parameters interactively
% We don't collect parameters when function in batch-mode (when known)
if isempty(p_nLag) || ~g_bBatchMode
    if isempty(p_nLag), p_nLag = 1; end % s

    cPrompt = {'Temporal lag (s)'};
    cAnswer = inputdlg(cPrompt,'Options', 1, ...
        {num2str(p_nLag)});
    if isempty(cAnswer), return, end
    p_nLag = str2double(cAnswer{1}); % s
end

% Cross-correlate signals
nMaxLagSamp = round(nFs / (1 / p_nLag)); % samples
vX = -(nMaxLagSamp / nFs):(1 / nFs):(nMaxLagSamp / nFs);

if 0
    % Average segments / trials separated by NaNs
    vNonNan = ~isnan(vContA);
    iNotNan = find(vNonNan);
    vStart = [iNotNan(1) find(diff(vNonNan) == 1)];
    vStop = find(diff(vNonNan) == -1);

    vContAi = intnans(vContA);
    vContBi = intnans(vContB);
    
    mC = [];
    for i = 1:(length(vStart)-1)
        [mC(end+1, :), ~] = xcov(vContAi(vStart(i):vStop(i)), vContBi(vStart(i):vStop(i)), nMaxLagSamp, 'coeff');
    end
    vC = mean(mC);
    vErr = std(mC);
else    
    % Remove NaN's by local mirroring
    vContA = intnans(vContA);
    vContB = intnans(vContB);
    [vC, ~] = xcov(detrend(vContA), detrend(vContB), nMaxLagSamp, 'coeff');
    vErr = zeros(size(vC));
end

%% Plot results
hFig = figure;
set(hFig, 'name', 'Cross Correlation')
hAx = subplot(1,1,1);%axes('position', [.125 .1 .8 .8], 'xscale', 'li', 'yscale', 'li');
hold(hAx, 'on');
hPatch = patch([vX fliplr(vX)], [vC+vErr fliplr(vC-vErr)], 'r', 'facealpha', .6);
plot(vX, vC, 'w')
ylabel('C')
xlabel('Lag (s)')
grid on
hTit = title(sprintf('%s vs %s', p_sContChA, p_sContChB));
Spiky.main.ThemeObject([hFig hAx hTit])
set(hTit, 'interpreter', 'none')
%%
return
