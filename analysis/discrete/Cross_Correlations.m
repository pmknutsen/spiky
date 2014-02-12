function Cross_Correlations(FV)
% Plot cross-correlations between units on same electrode
%
% Usage:
%   Cross_Correlations(FV)
%
% Supports batch processing.
%
% TODO:
%

global Spiky g_bBatchMode

% Select spiking channel
persistent p_sCh;
if isempty(p_sCh) || (~g_bBatchMode && nargout == 0)
    [p_sCh, bResult] = Spiky.main.SelectChannelNumber(fieldnames(FV.tSpikes)', 'Select spiking channel', p_sCh);
    if ~bResult, return, end
end

% Open figure
hFig = figure;
sDescr = Spiky.main.GetChannelDescription(p_sCh);
set(hFig, 'name', sprintf('Spiky - Unit Cross Correlations on %s', sDescr), 'NumberTitle', 'off');
Spiky.main.ThemeObject(hFig)
vXLim = [-.05 .05];

% Iterate over units (plot auto-corrs and cross-corrs)
vUnits = unique(FV.tSpikes.(p_sCh).hierarchy.assigns);
nFs = FV.tData.([p_sCh '_KHz']) * 1000;
for u1 = 1:length(vUnits)
    vIndx1 = FV.tSpikes.(p_sCh).hierarchy.assigns == vUnits(u1);
    vSpiketimes1 = FV.tSpikes.(p_sCh).spiketimes(vIndx1) / nFs; % sec
    % limit number of spikes to 10000
    nMaxLen = 10000;
    if length(vSpiketimes1) > nMaxLen
        vRand = randperm(length(vSpiketimes1));
        vSpiketimes1 = vSpiketimes1(vRand(1:nMaxLen));
    end
    for u2 = 1:length(vUnits)
        if u1 > u2, continue, end
        vIndx2 = FV.tSpikes.(p_sCh).hierarchy.assigns == vUnits(u2);
        vSpiketimes2 = FV.tSpikes.(p_sCh).spiketimes(vIndx2) / nFs; % sec
        % limit number of spikes to 10000
        if length(vSpiketimes2) > nMaxLen
            vRand = randperm(length(vSpiketimes2));
            vSpiketimes2 = vSpiketimes2(vRand(1:nMaxLen));
        end
        nX = .08; nY = .1; nW = .88; nH = .85; nNN = length(vUnits);
        hAx = axes('position', [nX+(u2-1)*(nW/nNN) (1-nY)-u1*(nH/nNN) nW/nNN nH/nNN], 'Color', [.1 .1 .1]);

        % Compute cross-correlation
        try
            % number of spikes must be same
            nNumSpikes = min([length(vSpiketimes1) length(vSpiketimes2)]);
            vRand = randperm(length(vSpiketimes1));
            vSpk1 = vSpiketimes1(vRand(1:nNumSpikes));
            vRand = randperm(length(vSpiketimes2));
            vSpk2 = vSpiketimes2(vRand(1:nNumSpikes));
            [vC, vT] = ccorr(sort(vSpk1)', sort(vSpk2)', vXLim, 'n', [], 1, 1, 200);
        catch
            continue
        end
        vT = vT*1000; % ms

        if all(vUnits([u1 u2]) == 0), vCol = [.5 .5 .5]; % outliers
        elseif u1 == u2, vCol = FV.mColors(u1,:);
        else vCol = [.7 .7 .7]; end
        plot(vT, vC, 'color', vCol)
        hold on
        plot([0 0], [0 max(vC)*2], ':', 'color', [.6 .6 .6])
        set(hAx, 'xlim', vXLim*1000, 'ylim', [0 max(vC)+.0001] )
        Spiky.main.ThemeObject(hAx);
        if u1 > u2, set(hAx, 'xtick', []); end
        if u1 == 1
            hTit = title(sprintf(' Unit %d ', vUnits(u2)));
            set(hTit, 'color', FV.mColors(u2,:), 'fontsize', 8, 'fontweight', 'bold', 'backgroundcolor', [.1 .1 .1])
        end
        xlabel('ms'); ylabel('')
        drawnow
    end
end

return
