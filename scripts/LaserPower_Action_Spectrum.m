function FV = LaserPower_Action_Spectrum(FV)
% Plot laser power vs spike-rate / spike probability
% Specifically, graph:
%   Rasters for all positions
%   Map image (with statistics)
%
% TODO:
%   Convert volts to W/mm2

% Check if any channel has been sorted
csFieldnames = fieldnames(FV.tSpikes)';
vDel = [];
for fn = 1:length(csFieldnames)
    if ~isfield(FV.tSpikes.(csFieldnames{fn}), 'hierarchy'), vDel(end+1) = fn; end
end
csFieldnames(vDel) = [];

% If no channel has been sorted, then pre-process DAQ file automatically
%if isempty(csFieldnames)
%    FV = Process_DAQ_File_Info(FV);
%end

% Get input parameters from user
persistent p_nWinDur p_nWV p_nBeamDiam

if isempty(p_nWinDur), p_nWinDur = 50; end % ms
if isempty(p_nWV), p_nWV = 2.75; end % mW
if isempty(p_nBeamDiam), p_nBeamDiam = 30; end % microns

cPrompt = {'Window duration (ms post-stim)', 'mW/V', 'Beam diameter (microns)' };
cAns = inputdlg(cPrompt, 'Action Spectrum', 1, ...
    {num2str(p_nWinDur), num2str(p_nWV), num2str(p_nBeamDiam)} );
if isempty(cAns) return, end
p_nWinDur = str2num(cAns{1}); % window
p_nWV = str2num(cAns{2}); % window
p_nBeamDiam = str2num(cAns{3}); % window

% Conversion factor form Volts to mW/mm^2
nWV = p_nWV / ((pi*(p_nBeamDiam/2)^2)/1000);

% Select channel
csFields = fieldnames(FV.tSpikes);
sFields = '{';
for f = 1:length(csFields)
    sFields = [sFields '''' csFields{f} ''' '];
end
sFields = [sFields '}'];
[sCh] = spiky(sprintf('SelectChannelNumber(%s)', sFields));
if isempty(sCh), return, end

% Iterate over units
vUnits = unique(FV.tSpikes.(sCh).hierarchy.assigns);
for u = 1:length(vUnits)
    
    % Initialize figure
    hFig = figure;
    %set(hFig, 'color', [.2 .2 .2], 'name', 'Spiky! Laser power vs spiking', 'NumberTitle', 'off')
    centerfig(hFig)
    drawnow
    
    nUnit = vUnits(u);
    
    % Get spiketimes of selected unit
    vSpkIndx = FV.tSpikes.(sCh).hierarchy.assigns == nUnit;
    vSpiketimes = FV.tSpikes.(sCh).spiketimes(vSpkIndx); % samples
    nTimebegin = FV.tData.(sprintf('%s_TimeBegin', sCh));
    nFs = FV.tData.(sprintf('%s_KHz', sCh)); % kHz
    vSpiketimes = vSpiketimes ./ (nFs*1000);% - nTimebegin;
    
    % Get laser ON times
    nIndx = strcmpi({FV.tChannelDescriptions.sDescription}, 'LaserShutter');
    sChLaser = FV.tChannelDescriptions(nIndx).sChannel;
    vUp = FV.tData.([sChLaser '_Up']); % sec
    
    % Get laserpower and time vectors
    sPowCh = FV.tChannelDescriptions(strcmp('LaserPower', {FV.tChannelDescriptions.sDescription})).sChannel;
    vLaserPow = FV.tData.(sPowCh); % V
    nFs = FV.tData.([sPowCh '_KHz']); % KHz
    nTB = FV.tData.([sPowCh '_TimeBegin']);
    vTime = linspace(nTB, nTB + length(vLaserPow)./(nFs.*1000), length(vLaserPow));

    % To get unique laser voltages;
    %  - round to nearest 0.01
    vLaserPow = floor(vLaserPow .* 100) ./ 100;
    %  - histogram in increments of millivolts
    [vN vX] = hist(vLaserPow, 0:.001:max(vLaserPow));
    %  - find modes (values) that occur more frequently than 0.001% of the time
    vUniqPowVals = vX( vN > length(vLaserPow) * .0001 );
    %  
    vX = 1:length(vUniqPowVals);
    vUniqPowVals = interp1([1 length(vUniqPowVals)], [min(vUniqPowVals) max(vUniqPowVals)], vX, 'linear');
    %mode(diff(sort(vUniqPowVals)))
    
    %  - substitute original trace with unique stimulation values
    vLaserPow = nearest(vLaserPow, vUniqPowVals);

    % UP times are times when the laser power is at any of the values above
    vUp = zeros(size(vTime));
    for i = 1:length(vUniqPowVals)
        if vUniqPowVals(i) == 0
            continue
        end
        vUp(vLaserPow == vUniqPowVals(i)) = 1;
    end
    vUp = vTime(find(diff(vUp) == 1));
    
    % Get laserpower at each laser ON time, delayed by N ms
    vLaserPowPulses = interp1(vTime, vLaserPow, vUp+(1/4000), 'nearest');
    
    % Round to nearest 0.05 V
    %vLaserPowPulses = round(vLaserPowPulses * 20) ./ 20;
    
    % Iterate over pulses
    vNumSpikes = zeros(length(vUp), 1);
    vSpikeLat = zeros(length(vUp), 1);
    vSpikeProb = zeros(length(vUp), 1);
    for i = 1:length(vUp)
        nUp = vUp(i);
        if i == length(vUp)
            nUpNext = vUp(i) + median(diff(vUp));
        else
            nUpNext = vUp(i+1);
        end
        vIndx = find(vSpiketimes >= nUp & vSpiketimes < nUp+(p_nWinDur/1000));
        % Spike count
        vNumSpikes(i) = length(vIndx);
        % Latency to 1st spike & Spike probability (0 or 1)
        if ~isempty(vIndx)
            vSpikeLat(i) = vSpiketimes(vIndx(1)) - vUp(i);
            vSpikeProb(i) = 1;
        else
            vSpikeLat(i) = NaN;
            vSpikeProb(i) = 0;
        end
    end

    % Plot
    subplot(2,1,1)
    [B,I,J] = unique(vLaserPowPulses);
    vCumSpikeCount = zeros(length(B), 1);
    vStimNum = zeros(length(B), 1);
    vAvgSpikeLat = zeros(length(B), 1);
    vStdSpikeProb = zeros(length(B), 1);
    for i = 1:length(B)
        vIndx = J == i;
        vCumSpikeCount(i) = sum(vNumSpikes(vIndx));
        vStdSpikeCount(i) = std(vNumSpikes(vIndx));
        vAvgSpikeProb(i) = mean(vSpikeProb(vIndx));
        vStdSpikeProb(i) = std(vSpikeProb(vIndx));
        vAvgSpikeLat(i) = nanmean(vSpikeLat(vIndx));
        vStimNum(i) = length(find(vIndx));
    end
    % Average spike count and probability
    vMean = vCumSpikeCount./vStimNum;
    vErr = vStdSpikeCount'./sqrt(vStimNum);
    errorbar(B, vMean, vErr); % mean +/- SEM
    xlabel('Laser-power (V)')
    ylabel('Spikecount (N)')
    %axis tight
    set(gca, 'xlim', [min(B)-.1 max(B)+.1], 'ylim', [-.1 max(vMean+vErr)+.1])

    % Plot laser-power vs spike latency
    subplot(2,2,3)
    hist(vSpikeLat.*1000, 50)
    xlabel('1st spike latency (ms)')
    ylabel('N')

    
    % Fits
    
    % - fit sigmoid (logistic function; sigmoid on linear scale)
    % [base height loc slope]
    %vSigF = @(p,x) p(1)^2 + p(2) ./ (1 + exp(-(x-p(3))/p(4)));
    %vParams = [0 max(vAvgSpikeCnt) mean(B) .1];
    
    %  - Hill equation (sigmoid on log scale)
    % [I_max n K]
    % [height slope pos]
    vSigF = @(p,x) p(1) .* ( x.^p(2) ./ [p(3)^p(2) + x.^p(2)]);

    % spiking probability does not work when there is spontaneous activity
    subplot(2,2,1)
    errorbar(B, vAvgSpikeProb, vStdSpikeProb./sqrt(vStimNum), 'ko') % mean +/- SEM

    hold on
    [vB, vR, vJ] = nlinfit(B, vAvgSpikeProb, vSigF, [max(vAvgSpikeProb) 1 mean(B) .1]); % fit
    vXX = linspace(0.001, max(B)+10, 1000);
    plot(vXX, vSigF(vB, vXX), 'k--')
    xlabel('Power (mW mm^-^2)')
    ylabel('p(spike)')
    nR2 = 1 - sum(vR.^2) / sum((vAvgSpikeProb - mean(vAvgSpikeProb)).^2); % compute R-square
    
    vXTickLabels = [.01 .05 .1:.1:1 2:5 10:10:100]; % mW
    vXTicks = vXTickLabels ./ nWV; % V
    set(gca, 'ylim', [-.05 1.05], 'xlim', [0.025 max(B)+min(B)], 'xtick', vXTicks, 'xtickLabel', vXTickLabels, 'xscale', 'log')
    
    % Find threshold P=0.75
    vX = linspace(min(B), max(B), 100);
    nThresh = interp1(vSigF(vB, vX), vX, .75, 'cubic');
    plot([nThresh nThresh], [0 1], 'r-')
    title(sprintf('R^2=%.2f, N=%d, p(%.2f)=.75', nR2, round(median(vStimNum)), nThresh))
    
    
    % For cases with high spontaneous acticity, the response function
    % should be estimated by spike counts, rather than spike probability
    subplot(2,2,2)
    vAvgSpikeCnt = vCumSpikeCount'./sqrt(vStimNum');
    errorbar(B, vAvgSpikeCnt, vStdSpikeCount, 'ko') % mean +/- SEM
    hold on

    [vB, vR, vJ] = nlinfit(B, vAvgSpikeCnt, vSigF, [max(vAvgSpikeCnt) 1 mean(B)]); % Hill fit
    plot(vXX, vSigF(vB, vXX), 'k--')
    xlabel('Power (mW mm^-^2)')
    ylabel('Spike count')
    nR2 = 1 - sum(vR.^2) / sum((vAvgSpikeCnt - mean(vAvgSpikeCnt)).^2); % compute R-square
    set(gca, 'ylim', [0 max(vAvgSpikeCnt)+1], 'xlim', [.025 max(B)+min(B)*2], 'xtick', vXTicks, 'xtickLabel', vXTickLabels, 'xscale', 'lo')
    title(sprintf('R^2=%.2f, N=%d', nR2, round(median(vStimNum))))
    
    
    % Plot laser vs latency to 1st spike
    subplot(2,2,4)
    vIndx = ~isnan(vSpikeLat);
    plot(vLaserPowPulses(vIndx), vSpikeLat(vIndx).*1000, 'k.')
    hold on
    plot(B, vAvgSpikeLat.*1000, 'ro')
    xlabel('Laser-power (V)')
    ylabel('1st spike latencym (ms)')
    % - regression
    [B,BINT,R,RINT,STATS] = regress(vSpikeLat(vIndx).*1000, [vLaserPowPulses(vIndx)' ones(length(find(vIndx)),1)]);
    vX = [min(vLaserPowPulses) max(vLaserPowPulses)];
    vY = vX .* B(1) + B(2);
    plot(vX, vY, 'r-')
    axis tight

end

return


function [vC, vInd] = nearest(vA, vB)
% NEAREST Substitute for nearest matrching number
% [C, I] = NEAREST(A,B) substitutes numbers in vector A with the nearest number
% in B and returns the result to C, where the length of C is equal to the
% length of A. The substitution indices of B are returned in I.
%
% Written by Per Magne Knutsen, January 2005

% Reshape vectors
vA = vA(:);
vB = vB(:)';

% Subtract every number in vB from every number in vA
mB = repmat(vB, length(vA), 1);
mA = repmat(vA, 1, length(vB));
mR = mA - mB;
[vMin, vInd] = min(abs(mR),[],2);
vC = vB(vInd);

return