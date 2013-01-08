function FV = Plot_S1_Input_Map(FV)
% Display the S1 input map to SpVi units
% Specifically, graph:
%   Rasters for all positions
%   Map image (with statistics)
%
% TODO
%   Plot maps on top of vessel image
%   Plot maps on top of barrel map (from IOS)
%   Change unit of smoothing parameter to millimeters
%   Add prestim and poststim parameters
%   Add scaling factors to image titles (max spike-rate or max z-score)
%

% Check if any channel has been sorted
csFieldnames = fieldnames(FV.tSpikes)';
vDel = [];
for fn = 1:length(csFieldnames)
    if ~isfield(FV.tSpikes.(csFieldnames{fn}), 'hierarchy'), vDel(end+1) = fn; end
end
csFieldnames(vDel) = [];
% Abort if no channels have been sorted
if isempty(csFieldnames)
    uiwait(warndlg('You need to spike sort channels first.', 'modal'))
    return;
%    FV = Process_DAQ_File_Info(FV); % process automatically TODO
end

persistent p_nSmooth p_nWinStart p_nWinEnd p_nPreStimPeriod p_nPostStimPeriod p_bShowZScore
persistent p_vSearchStrSel p_nZCutOff p_bRemoveInfreqPos
csStr = {'_ScanPath_' '_ScanPathContAir_'};
global g_bMergeMode g_bPlotS1MapQuiet;
if isempty(g_bPlotS1MapQuiet), g_bPlotS1MapQuiet = 0; end

bDebug = 0;

% Default parameters
if isempty(p_nSmooth), p_nSmooth = 0; end
if isempty(p_nWinStart), p_nWinStart = 0.005; end
if isempty(p_nWinEnd), p_nWinEnd = 0.04; end
if isempty(p_nPreStimPeriod), p_nPreStimPeriod = 0.025; end
if isempty(p_nPostStimPeriod), p_nPostStimPeriod = 0.1; end
if isempty(p_bShowZScore), p_bShowZScore = 0; end
if isempty(p_vSearchStrSel), p_vSearchStrSel = 1; end
if isempty(p_nZCutOff), p_nZCutOff = 0; end
if isempty(p_bRemoveInfreqPos), p_bRemoveInfreqPos = 1; end

if ~isfield(FV, 'S1Map_GetROI') && ~g_bPlotS1MapQuiet
    cPrompt = {'Smoothing (sigma)','Evoked activity window start (s)', 'Evoked activity window end (s)', ...
        'Pre-stimulus period (s)', 'Post-stimulus period (s)', 'Show z-scores (1=yes)', ...
        'Z-score cutoff (SDs)', 'Remove infrequent positions (1=yes)' };
    cAns = inputdlg(cPrompt, 'GalvoScan Map', 1, ...
        {num2str(p_nSmooth), num2str(p_nWinStart), num2str(p_nWinEnd), num2str(p_nPreStimPeriod), ...
        num2str(p_nPostStimPeriod), num2str(p_bShowZScore), num2str(p_nZCutOff), num2str(p_bRemoveInfreqPos)} );
    if isempty(cAns) return, end
    p_nSmooth = str2num(cAns{1}); % sigma of gaussian smoothing window
    p_nWinStart = str2num(cAns{2}); % sec
    p_nWinEnd = str2num(cAns{3}); % sec
    p_nPreStimPeriod = str2num(cAns{4}); % sec
    p_nPostStimPeriod = str2num(cAns{5}); % sec
    p_bShowZScore = str2num(cAns{6}); % sec
    p_nZCutOff = str2num(cAns{7}); % sec
    p_bRemoveInfreqPos = str2num(cAns{8}); % sec
    
    % Ask under what conditions maps should be generated
    % Options are: ScanPath, ScanPathContAir or Both
    [p_vSearchStrSel, nOK] = listdlg('PromptString', 'Map conditions:',...
        'ListString', csStr, ...
        'ListSize', [200 50], ...
        'SelectionMode','multiple', ...
        'InitialValue', p_vSearchStrSel);
    if ~nOK, return; end
end

% Select channel
csFields = fieldnames(FV.tSpikes);
if length(csFields) > 1
    sFields = '{';
    for f = 1:length(csFields)
        sFields = [sFields '''' csFields{f} ''' '];
    end
    sFields = [sFields '}'];
    [sCh] = spiky(sprintf('SelectChannelNumber(%s)', sFields));
    if isempty(sCh), return, end
else
    sCh = csFields{1};
end

hFig = figure;
set(hFig, 'color', [.2 .2 .2], 'name', 'Spiky! GalvoScan Map', 'NumberTitle', 'off', ...
    'position', [200 150 950 550])
centerfig(hFig)
drawnow
colormap(fireice(1024))

% Get [start end] timestamps of files that are to be included
if isfield(FV.tData, 'FileStart')
    csSearchStr = csStr(p_vSearchStrSel);
    mInclRanges = [];
    for f = 1:length(FV.tData.FileStart)
        for s = 1:length(csSearchStr)
            sSearchStr = csSearchStr{s};
            if ~isempty(strfind(FV.tData.FileStart(f).File, sSearchStr))
                nStart = FV.tData.FileStart(f).Timestamp;
                nEnd = FV.tData.FileEnd(f).Timestamp;
                mInclRanges = [mInclRanges; nStart nEnd];
            end
        end
    end
    mInclRanges = unique(mInclRanges, 'rows'); % s
else mInclRanges = [0 inf]; end

% Iterate over units
vUnits = unique(FV.tSpikes.(sCh).hierarchy.assigns);
for u = 1:length(vUnits)
    nUnit = vUnits(u);
    
    % Get spiketimes of selected unit
    vSpkIndx = FV.tSpikes.(sCh).hierarchy.assigns == nUnit;
    vSpiketimes = FV.tSpikes.(sCh).spiketimes(vSpkIndx); % samples
    nTimebegin = FV.tData.(sprintf('%s_TimeBegin', sCh)); % samples
    
    nFs = FV.tData.(sprintf('%s_KHz', sCh)); % kHz
    vSpiketimes = vSpiketimes ./ (nFs*1000) - nTimebegin; % s
    
    % Get laser ON times
    nIndx = strcmpi({FV.tChannelDescriptions.sDescription}, 'LaserShutter');
    sChLaser = FV.tChannelDescriptions(nIndx).sChannel;
    vUp = FV.tData.([sChLaser '_Up']); % sec
    vUp = vUp - nTimebegin; % sec
    
    % Get galvo scanner position
    vAPCont = FV.tData.('GalvoScanAP'); % mm
    vMLCont = FV.tData.('GalvoScanML'); % mm
    nFs = FV.tData.('GalvoScanML_KHz'); % KHz
    
    % Drop spikes and event times that are not in the ranges specified by mInclRanges
    vSpkKeepIndx = [];
    vUpKeepIndx = [];
    for r = 1:size(mInclRanges, 1)
        vSpkKeepIndx = unique([vSpkKeepIndx; ...
            find((vSpiketimes >= mInclRanges(r, 1)) & (vSpiketimes <= mInclRanges(r, 2)))]);
        vUpKeepIndx = unique([vUpKeepIndx ...
            find((vUp >= mInclRanges(r, 1)) & (vUp <= mInclRanges(r, 2)))]);
    end
    vSpiketimes = vSpiketimes(vSpkKeepIndx); % s
    vUp = vUp(vUpKeepIndx); % s
        
    % Construct time vector
    if ~g_bMergeMode
        % Single file analysis
        vTime = linspace(vUp(1), length(vMLCont)./(nFs.*1000) + vUp(1), length(vMLCont)); % sec
    else
        % Analysis of merge files
        nBeginTime = FV.tData.('GalvoScanML_TimeBegin'); % start of sampling (sec)
        nEndTime = FV.tData.('GalvoScanML_TimeEnd'); % start of sampling (sec)
        vTime = linspace(nBeginTime, nEndTime, length(vAPCont));
    end
    
    % Get [AP ML] coordinate at each light pulse
    vAPCont(isnan(vAPCont)) = interp1(vTime(~isnan(vAPCont)), vAPCont(~isnan(vAPCont)), ...
        vTime(isnan(vAPCont)), 'nearest');
    vMLCont(isnan(vMLCont)) = interp1(vTime(~isnan(vMLCont)), vMLCont(~isnan(vMLCont)), ...
        vTime(isnan(vMLCont)), 'nearest');
    vAP = interp1(vTime, vAPCont, vUp, 'nearest');
    vML = interp1(vTime, vMLCont, vUp, 'nearest');

    vAP(isnan(vAP)) = [];
    vML(isnan(vML)) = [];
    
    % Check if unique values are regularly spaced (e.g. in a grid). If they
    % are, then find the mode of the grid spacing and round all positions
    % to the nearest vertex. This is a fix that removes rarely occurring values.
    if p_bRemoveInfreqPos
        vUap = unique(vAP);
        vUml = unique(vML);
        % Number of repetitions we assume for each unique position
        vNthresh = length(vAP) / ( length(vUap) * length(vUml) );
        vRemIndx = [];
        vL = [];
        for uap = 1:length(vUap)
            for uml = 1:length(vUml)
                % If number of repetitions for current position is below
                % threshold, then remove this position
                nL = length( find( vAP == vUap(uap) & vML == vUml(uml) ) );
                vL(end+1) = nL;
                if nL < (vNthresh*.75) || nL > (vNthresh*4)
                    % Remove all events related to current position
                    vRemIndx = [vRemIndx find( vAP == vUap(uap) & vML == vUml(uml))];
                end
            end
        end
        vAP(vRemIndx) = [];
        vML(vRemIndx) = [];
        vUp(vRemIndx) = [];
    end
    
    mPos = [vAP' vML'];
    
    if bDebug
        figure; hold on
        plot(vTime, vAPCont,'k-')
        scatter(vUp, vAP, 'co', 'filled')
        plot(vSpiketimes, zeros(size(vSpiketimes)), 'r.')
        set(gca, 'xlim', [min(vUp) max(vUp)])
        xlabel('Time (s)')
        legend('AP', 'On', 'Spikes')
    end
    
    % Iterate over all positions
    cSpikes = cell(size(mPos, 1), 1);
    for i = 1:size(mPos, 1)
        % Stim time
        nUp = vUp(i);
        
        % Find all spikes that occurred between p_nPreStimPeriod and
        % p_nPostStimPeriod (i.e. windows of PSTH)
        vIndx = vSpiketimes >= (nUp - p_nPreStimPeriod) & vSpiketimes <= (nUp + p_nPostStimPeriod);
        vSpikes = vSpiketimes(vIndx) - nUp; % spiketimes relative to stim
        
        % Keep spikes
        cSpikes{i} = vSpikes;
    end
    
    % Get list of all unique stimulus positions
    mPosUniq = unique(mPos, 'rows');

    % Create lookup table for converting between positions and matrix indices
    vAPLookup = unique(mPos(:,1));
    vMLLookup = unique(mPos(:,2));
    
    % Generate spike-rate map (#spikes in window)
    mImg = zeros(sqrt(size(mPosUniq, 1)));
    mImgBaseline = zeros(sqrt(size(mPosUniq, 1)));
    
    vStimSpikesInsideROIAll = [];
    vStimSpikesOutsideROIAll = [];
    for c = 1:size(mPosUniq, 1)
        % Get matrix location
        nI = mPosUniq(c, 1) == vAPLookup;
        nJ = mPosUniq(c, 2) == vMLLookup;
        
        % Get all repetitions at current position
        vIndx = find(all(mPos == repmat(mPosUniq(c, :), size(mPos, 1), 1), 2));
        vStimSpikes = [];
        vPreStimSpikes = [];
        for i = 1:length(vIndx)
            % Evoked spikes
            vSpikes = cSpikes{vIndx(i)}(cSpikes{vIndx(i)} >= p_nWinStart & cSpikes{vIndx(i)} <= p_nWinEnd);
            vStimSpikes = [vStimSpikes(:); vSpikes];
            
            % Pre-stim spikes
            vSpikes = cSpikes{vIndx(i)}(cSpikes{vIndx(i)} <= 0);
            vPreStimSpikes = [vPreStimSpikes(:); vSpikes];
            
            % Post-stim spikes (post-stim is taken as from one stim-window duration
            % after the end of the stim window until the end of the trial)
            vSpikes = cSpikes{vIndx(i)}(cSpikes{vIndx(i)} >= p_nWinEnd);
            vPostStimSpikes = [vPreStimSpikes(:); vSpikes];
        end

        % Evoked firing rate (spikes/sec)
        mImg(nJ, nI) = length(vStimSpikes) / (p_nWinStart+p_nWinEnd) / i;
        
        % Baseline firing rate (spikes/sec)
        nPreRate = length(vPreStimSpikes) / (p_nPreStimPeriod) / i;
        nPostRate = length(vPostStimSpikes) / (p_nPostStimPeriod - (p_nWinEnd*2)) / i;
        mImgBaseline(nJ, nI) = (nPreRate + nPostRate) / 2;
    end
    
    if p_bShowZScore
        % Compute z-score matrix
        mImg = (mImg - mean(mImgBaseline(:))) ./ std(mImgBaseline(:));

        
        %imagesc( (mImg - mean(mImgBaseline(:))) ./ [mean(mImg(:)) / sqrt()] )
        
        % Z-Test
        mImg(mImg(:) < p_nZCutOff & mImg(:) > -p_nZCutOff) = 0;
    else
        % Subtract baseline and normalize to peak (absolute)
        mImg = mImg - mean(mImgBaseline(:));
        nMaxSpikeRate = max(mImg(:));
        mImg = mImg ./ max(abs(mImg(:)));
    end
    
    % Display intensity map
    if ~ishandle(hFig), return; end
    figure(hFig)
    hImAx = subplot(2, length(vUnits), u);
    if p_nSmooth > 0
        hImg = imagesc(filt2(mImg, p_nSmooth, 'lm'));
    else
        hImg = imagesc(mImg);
    end

    % Plot ROI
    if isfield(FV, 'S1Map_GetROI')
        hold on
        plot(FV.S1Map_ROIx, FV.S1Map_ROIy, 'g.:')
    end
    
    % Set axis properties
    axis square
    if p_bShowZScore
        nA = max(abs(mImg(:)));
        if isnan(nA)
            set(hImAx, 'tickdir', 'out', 'fontsize', 7, 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'clim', [-.1 .1])
        else
            set(hImAx, 'tickdir', 'out', 'fontsize', 7, 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'clim', [-nA max([0.01 nA])])
        end
    else
        set(hImAx, 'tickdir', 'out', 'fontsize', 7, 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'clim', [-1 1])
    end
    vAxPos = get(hImAx, 'position');
    
    % Place context menu on axis
    hMenu = uicontextmenu;
    set(hMenu, 'UserData', FV)
    hItems(1) = uimenu(hMenu, 'Label', '&Re-compute PSTH with ROI', 'Callback', @ReRunWithROI, 'Tag', 'Spiky_GetROI');
    set(get(hImAx, 'child'), 'uicontextmenu', hMenu)
    
    % Estimate centroid (location of peak activation)
    % Simple estimate is find pixels with vals > 2* std, threshold and find median XY of points
    mThresh = zeros(size(mImg));
    mThresh(mImg >= 2*std(mImg(:))) = 1;
    [vY, vX] = find(mThresh);
    nX = round(median(vX));
    nY = round(median(vY));

    % Plot 'centroid'
    hold on
    plot(nX, nY, 'wx', 'markersize', 12)

    % Use actual anatomical coordiates as axis labels
    vXlabels = unique(mPos(:,1));
    vYlabels = unique(mPos(:,2));

    % Get anatomical coordinates of centroid
    if ~isnan(nX)
        nCentroidX = vXlabels(nX);
        nCentroidY = vYlabels(nY);
    else
        nCentroidX = NaN;
        nCentroidY = NaN;
    end
    
    vXInterp = -5:.5:5;
    vXInterpVals = interp1(vXlabels, 1:length(vXlabels), vXInterp);
    set(hImAx, 'xtick', vXInterpVals, 'xtickLabel', vXInterp, 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6])
    vYInterp = -5:.5:5;
    vYInterpVals = interp1(vYlabels, 1:length(vYlabels), vYInterp);
    set(hImAx, 'ytick', vYInterpVals, 'ytickLabel', vYInterp, 'ydir', 'normal', 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6])
    grid
    
    % Set title
    if isfield(FV.tSpikes.(sCh), 'quality')
        nUnitIndx = find([FV.tSpikes.(sCh).quality.unit] == nUnit);
        if ~isempty(nUnitIndx)
            nQuality = FV.tSpikes.(sCh).quality(nUnitIndx).score;
        else nQuality = NaN; end
    else nQuality = NaN; end
    
    if isfield(FV.tSpikes.(sCh), 'receptivefield')
        nUnitIndx = find([FV.tSpikes.(sCh).receptivefield.unit] == nUnit);
        if isempty(FV.tSpikes.(sCh).receptivefield(nUnitIndx))
            sRF = 'Unknown';
        else
            sRF = FV.tSpikes.(sCh).receptivefield(nUnitIndx).rf;
        end
    else sRF = 'NA'; end

    [nMax, nMaxI] = max(abs(mImg(:)));
    if p_bShowZScore
        hTit = title(sprintf('U:%d  Q:%d  RF:%s  MaxZ:%.2f', nUnit, nQuality, sRF, mImg(nMaxI)));
    else
        hTit = title(sprintf('U:%d  Q:%d  RF:%s  MaxI:%.2f', nUnit, nQuality, sRF, nMaxSpikeRate));
    end
    set(hTit, 'color', 'w', 'interpreter', 'none', 'fontsize', 8);
    
    % Plot barrels on top of map (if available)
    if isfield(FV.tData, 'BarrelOutlines')
        if isfield(FV.tData.BarrelOutlines, 'vXmm')
            vInternalInt = [];
            vBoundaryInt = [];
            for b = 1:length(FV.tData.BarrelOutlines)
                axes(hImAx); hold on
                vX = FV.tData.BarrelOutlines(b).vXmm;
                vY = FV.tData.BarrelOutlines(b).vYmm;
                % Contour
                vXi = interp1(vXlabels, 1:length(vXlabels), vY, 'linear', 'extrap');
                vYi = interp1(vYlabels, 1:length(vYlabels), vX, 'linear', 'extrap');
                
                % Emphasize contour if barrel is in cell's receptive field
                cIDs = textscan(sRF,'%s','Delimiter',',');
                cIDs = squeeze(cIDs{1});
                % Compare only first two letters
                for c = 1:length(cIDs)
                    cIDs{c} = cIDs{c}(1:2);
                end
                if any(strcmpi(cIDs, FV.tData.BarrelOutlines(b).sID))
                    plot(vXi, vYi, 'g-', 'linewidth', 3)
                else
                    plot(vXi, vYi, 'w-', 'linewidth', 1)
                end
                
                % For each pixel enclosed by barrel, determine its value
                % and distance from edge
                [vYY, vXX] = ind2sub(size(mImg), 1:prod(size(mImg)));
                [vIN, vON] = inpolygon(vXX, vYY, round(vXi), round(vYi));
                mOnMat = zeros(size(mImg));
                mOnMat(vON) = 1;
                mInMat = zeros(size(mImg));
                mInMat(vIN) = 1;
                mInMat = mInMat - mOnMat;
                % Average value inside barrel and on boundary
                vInternalInt = [vInternalInt; mImg(logical(mInMat))];
                vBoundaryInt = [vBoundaryInt; mImg(logical(mOnMat))];
            end
        end
    end
    
    if ~ishandle(hFig), return; end
    figure(hFig)
    
    % Compute PSTH
    if ~isfield(FV, 'S1Map_GetROI') FV.S1Map_Mask = ones(size(mImg)); end
    vStimSpikesInsideROI = [];
    vStimSpikesOutsideROI = [];
    for c = 1:size(mPosUniq, 1)
        % Get matrix location
        nI = mPosUniq(c, 1) == vAPLookup;
        nJ = mPosUniq(c, 2) == vMLLookup;
        
        % Iterate over all repetitions at current position
        vIndx = find(all(mPos == repmat(mPosUniq(c, :), size(mPos, 1), 1), 2));
        for i = 1:length(vIndx)
            % Determine if current position is inside or outside ROI
            if FV.S1Map_Mask(nJ, nI)
                vStimSpikesInsideROI = [vStimSpikesInsideROI(:); cSpikes{vIndx(i)}]; % spikes inside mask
            else
                vStimSpikesOutsideROI = [vStimSpikesOutsideROI(:); cSpikes{vIndx(i)}]; % spikes outside mask
            end
        end
    end
    
    % Plot PSTH
    if ~ishandle(hFig), return; end
    figure(hFig)
    hAx = subplot(2, length(vUnits), u+(length(vUnits))); hold on

    nBin = 0.001;
    [vN, vX] = hist(vStimSpikesInsideROI, (-p_nPreStimPeriod):nBin:p_nPostStimPeriod);
    vN_psth = vN./nBin./length(cSpikes);
    
    % Draw thick bar to denote evoked activity window
    hPatch = patch([p_nWinStart p_nWinEnd p_nWinEnd p_nWinStart p_nWinStart], [0 0 max(vN_psth) max(vN_psth) 0], [.7 0 0]);
    set(hPatch, 'edgeColor', [.7 0 0])
    
    % Plot PSTH
    bar(vX, vN_psth, 'facecolor', [.8 .8 .8], 'edgecolor', [.8 .8 .8]);
    
    % Plot baseline activity on PSTH
    nBaseline = mean(mImgBaseline(:));
    plot([-p_nPreStimPeriod p_nPostStimPeriod], [nBaseline nBaseline], 'g--')
    set(hAx, 'xlim', [-p_nPreStimPeriod p_nPostStimPeriod], 'color', [.1 .1 .1], 'xcolor', [.6 .6 .6], ...
        'ycolor', [.6 .6 .6], 'ylim', [0 max([0.001 max(vN_psth)])])
    xlabel('Post-stimulus time (s)')
    ylabel('Spikes/s')
    hTit = title(sprintf('PSTH N=%d', length(vStimSpikesInsideROI)));
    set(hTit, 'color', 'w', 'interpreter', 'none');
    axis tight square
    drawnow

end % end of unit loop

% Create colorbar
axes(hImAx)
try
    hCol = colorbar;
    vPos = get(hImAx, 'position');
    set(hCol, 'Position', [.93 vPos(2) .02 vPos(4)], 'fontsize', 7, 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6]);
end

% Bring last barrel map to top
if exist('hAxBarrels')
    axes(hAxBarrels)
end

hHeader = header(['Channel ' sCh ' - GalvoScanner Map'], 10);
hFooter = footer(fullfile(FV.sDirectory, FV.sLoadedTrial), 10);
set([hHeader hFooter], 'color', 'w', 'interpreter', 'none')

FV = rmfield(FV, 'S1Map_Mask');

% If g_bPlotS1MapQuiet is true, we are presumably running this function
% through some automated scripting. In such case, close the figure and make
% available a global structure which contains the S1 map and PSTH.
if g_bPlotS1MapQuiet
    global g_S1MapAnalysis
    g_S1MapAnalysis = struct([]);
    
    % Parameters
    g_S1MapAnalysis(1).p_nSmooth = p_nSmooth;
    g_S1MapAnalysis(1).p_nWinStart = p_nWinStart;
    g_S1MapAnalysis(1).p_nWinEnd = p_nWinEnd;
    g_S1MapAnalysis(1).p_nPreStimPeriod = p_nPreStimPeriod;
    g_S1MapAnalysis(1).p_nPostStimPeriod = p_nPostStimPeriod;
    g_S1MapAnalysis(1).p_bShowZScore = p_bShowZScore;
    g_S1MapAnalysis(1).p_vSearchStrSel = p_vSearchStrSel;
    g_S1MapAnalysis(1).p_nZCutOff = p_nZCutOff;
    g_S1MapAnalysis(1).p_bRemoveInfreqPos = p_bRemoveInfreqPos;
    
    % Results
    g_S1MapAnalysis(1).map = mImg;
    g_S1MapAnalysis(1).map_xcoords_mm = vXlabels;
    g_S1MapAnalysis(1).map_ycoords_mm = vYlabels;
    g_S1MapAnalysis(1).psth_spikes = vStimSpikesInsideROI;
    g_S1MapAnalysis(1).psth_bins = (-p_nPreStimPeriod):nBin:p_nPostStimPeriod;
    g_S1MapAnalysis(1).psth_baseline = nBaseline;
    g_S1MapAnalysis(1).psth_rate = vN_psth; % normalized spikes/sec
    g_S1MapAnalysis(1).map_max = nMaxSpikeRate; % max spikes/sec
    
    close(hFig);
end

return


function ReRunWithROI(hObject, handles)
% Draw region of interest (polygon) on current axis
[mMask, vX, vY] = roipoly();
FV = get(get(hObject, 'Parent'), 'UserData');
FV.S1Map_GetROI = 1;
FV.S1Map_ROIx = vX;
FV.S1Map_ROIy = vY;
FV.S1Map_Mask = mMask;
eval(sprintf('FV=%s(FV);', mfilename))
return


function cmap = fireice(m)
% FIREICE LightCyan-Cyan-Blue-Black-Red-Yellow-LightYellow Colormap
%
%  FIREICE(M) Creates a colormap with M colors
%
%   Inputs:
%       M - (optional) an integer between 1 and 256 specifying the number
%           of colors in the colormap. Default is 64.
%
%   Outputs:
%       CMAP - an Mx3 colormap matrix
%
%   Example:
%       imagesc(peaks(500))
%       colormap(fireice), colorbar
%
%   Example:
%       imagesc(interp2(rand(10),2))
%       colormap(fireice); colorbar
%
% See also: hot, jet, hsv, gray, copper, bone, cold, vivid, colors
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
% Release: 1.0
% Date: 07/29/09


% Default Colormap Size
if ~nargin
    m = 64;
end

% LightCyan-Cyan-Blue-Black-Red-Yellow-LightYellow
clrs = [0.75 1 1; 0 1 1; 0 0 1;...
    0 0 0; 1 0 0; 1 1 0; 1 1 0.75];

y = -3:3;
if mod(m,2)
    delta = min(1,6/(m-1));
    half = (m-1)/2;
    yi = delta*(-half:half)';
else
    delta = min(1,6/m);
    half = m/2;
    yi = delta*nonzeros(-half:half);
end
cmap = interp2(1:3,y,clrs,1:3,yi);
return
