function FV = GalvoScanner_Continuous_Response_Map(FV, varargin)
% Display responses on a continuous channel to stimulation of S1/M1 in a
% grid-pattern. This function was adapted from Plot_S1_Input_Map, but plots
% the amplitude in a certain time window of a continuous signal, rather
% than spike responses.
% 
% 
% This script can be run silently, using the default parameters, by setting
% the global variable g_bPlotS1MapQuiet to true (or 1).
% 
% If g_bPlotS1MapQuiet is used, the default parameters can also be changed
% by passing a structure as a second input. This structure should contain
% fields with names that match any of the persistent variables and the
% corresponding values to use. For instance:
% 
%   tParams = struct('p_nSmooth', 0.5, 'p_vSearchStrSel', [1 2])
% 
%       FV = GalvoScanner_Continuous_Response_Map(FV, tParams)
%

%
% TODO
%   Add option to click on pixel and show corresponding curve in new color
%   in lower subplot
%   Add option to show multiple time slices (redundant with above fix?)
%   Add option to display map over an image (eg blue image)
%   Better method to get both negative and positive amplitude
%       Average
%       Two maps, next to each other: min and max
%

global g_bPlotS1MapQuiet Spiky g_S1MapAnalysis g_bMergeMode g_bBatchMode
persistent p_nSmooth p_nTempSmooth p_nWinStart p_nWinEnd p_nPreStimPeriod p_nPostStimPeriod p_vSearchStrSel p_bRemoveInfreqPos p_sContCh

bPlotETA = 0;

if g_bPlotS1MapQuiet
    g_S1MapAnalysis = struct([]);
end

% Check that data is loaded and that we have the required data fields
if ~Spiky.main.CheckDataLoaded(), return; end
if ~isfield(FV.tData, 'GalvoScanAP')
    Spiky.main.sp_disp('GalvoScanner data is missing from this file.');
    return
end

% Select continuous channel to measure responses from
if isempty(p_sContCh) || ~g_bBatchMode
    [p_sContCh, bResult] = Spiky.main.SelectChannelNumber(FV.csDisplayChannels, 'Select continuous signal', p_sContCh);
    if ~bResult; return, end
end

csStr = {'_ScanPath_' '_ScanPathContAir_'};
if isempty(g_bPlotS1MapQuiet), g_bPlotS1MapQuiet = 0; end

if g_bPlotS1MapQuiet
    disp('Plot_S1_Response_Map is in Quiet mode. To disable: global g_bPlotS1MapQuiet; g_bPlotS1MapQuiet = 0;')
end
bDebug = 0;

% Default parameters
if isempty(p_nSmooth), p_nSmooth = 0; end
if isempty(p_nTempSmooth), p_nTempSmooth = 0; end
if isempty(p_nWinStart), p_nWinStart = 0.01; end
if isempty(p_nWinEnd), p_nWinEnd = 0.1; end
if isempty(p_nPreStimPeriod), p_nPreStimPeriod = 0.5; end
if isempty(p_nPostStimPeriod), p_nPostStimPeriod = 0.5; end
if isempty(p_vSearchStrSel), p_vSearchStrSel = 1; end % eg [1] or [1 2] for multiple selections
if isempty(p_bRemoveInfreqPos), p_bRemoveInfreqPos = 1; end

% Change default parameters if a second input parameters is passed
if nargin == 2
    tParams = varargin{1};
    if isstruct(tParams)
        csFieldnames = fieldnames(tParams);
        for i = 1:length(csFieldnames)
            eval(sprintf('%s = tParams.%s;', csFieldnames{i}, csFieldnames{i}));
        end
    end
end

if ~isfield(FV, 'S1Map_GetROI') && ~g_bPlotS1MapQuiet
    cPrompt = {'Spatial smoothing (sigma mm)', 'Temporal smoothing (ms)'...
        'Response window start (s)', 'Response window end (s)', ...
        'Display period pre-stim (s)', 'Display period post-stim (s)', ...
        'Remove infrequent grid locations (1=yes)' };
    cAns = inputdlg(cPrompt, 'GalvoScan Map', 1, ...
        {num2str(p_nSmooth), num2str(p_nTempSmooth), num2str(p_nWinStart), num2str(p_nWinEnd), num2str(p_nPreStimPeriod), ...
        num2str(p_nPostStimPeriod), num2str(p_bRemoveInfreqPos)} );
    if isempty(cAns) return, end
    p_nSmooth = str2num(cAns{1}); % sigma of gaussian smoothing window
    p_nTempSmooth = str2num(cAns{2}); % sec
    p_nWinStart = str2num(cAns{3}); % sec
    p_nWinEnd = str2num(cAns{4}); % sec
    p_nPreStimPeriod = str2num(cAns{5}); % sec
    p_nPostStimPeriod = str2num(cAns{6}); % sec
    p_bRemoveInfreqPos = str2num(cAns{7}); % sec
    
    % Ask under what conditions maps should be generated
    % Options are: ScanPath, ScanPathContAir or Both
    [p_vSearchStrSel, nOK] = listdlg('PromptString', 'Map conditions:',...
        'ListString', csStr, ...
        'ListSize', [200 50], ...
        'SelectionMode','multiple', ...
        'InitialValue', p_vSearchStrSel);
    if ~nOK, return; end
end
sCh = p_sContCh;

% Initialize figure
hFig = figure('visible', 'off');
set(hFig, 'name', 'GalvoScanner Grid Response Map', 'NumberTitle', 'off', 'position', [200 150 420 600])
Spiky.main.ThemeObject(hFig)
centerfig(hFig, Spiky.main.GetGUIHandle())
if ~g_bPlotS1MapQuiet, set(hFig, 'visible', 'on'); end
drawnow

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

% Abort now if there are zero data inclusion ranges (typically means that
% the file inclusion string is not relevant to this dataset.
if isempty(mInclRanges), return; end

% Get the continuous signal and associated parameters
vCont = FV.tData.(sCh);
vCont = Spiky.main.ChannelCalculator(vCont, sCh); % apply math filter

nContBeginTime = FV.tData.([sCh '_TimeBegin']); % samples
nContEndTime = FV.tData.([sCh '_TimeEnd']); % samples
nContFs = FV.tData.([sCh '_KHz']) * 1000; % Hz
vContTime = Spiky.main.GetTime(nContBeginTime, nContEndTime, vCont, nContFs); % s, absolute time

% Get scanning laser event times
nIndx = strcmpi({FV.tChannelDescriptions.sDescription}, 'LaserShutter');
if isempty(nIndx)
    Spiky.main.sp_disp('LaserShutter channel is missing from this file.');
    return
end
sChLaser = FV.tChannelDescriptions(nIndx).sChannel;
vUp = FV.tData.([sChLaser '_Up']); % s, absolute time
vUpFs = FV.tData.([sChLaser '_KHz']) * 1000; % Hz


%% Get a list of unique GalvoScanner locations
% Get GalvoScanner data
vAPCont = FV.tData.('GalvoScanAP'); % mm
vMLCont = FV.tData.('GalvoScanML'); % mm
nGalvScanFs = FV.tData.('GalvoScanML_KHz'); % KHz

% Construct time vector for GalvoScanner data
if ~g_bMergeMode
    % Single file analysis
    vGalvScanTime = linspace(vUp(1), length(vMLCont)./(nGalvScanFs.*1000) + vUp(1), length(vMLCont)); % sec
else
    % Analysis of merge files
    nGalvScanBeginTime = FV.tData.('GalvoScanML_TimeBegin'); % start of sampling (sec)
    nGalvScanEndTime = FV.tData.('GalvoScanML_TimeEnd'); % start of sampling (sec)
    vGalvScanTime = linspace(nGalvScanBeginTime, nGalvScanEndTime, length(vAPCont));
end

% Get [AP ML] coordinate at each light pulse & create continuous signals
vAPCont(isnan(vAPCont)) = interp1(vGalvScanTime(~isnan(vAPCont)), vAPCont(~isnan(vAPCont)), ...
    vGalvScanTime(isnan(vAPCont)), 'nearest');
vMLCont(isnan(vMLCont)) = interp1(vGalvScanTime(~isnan(vMLCont)), vMLCont(~isnan(vMLCont)), ...
    vGalvScanTime(isnan(vMLCont)), 'nearest');
vAP = interp1(vGalvScanTime, vAPCont, vUp, 'nearest');
vML = interp1(vGalvScanTime, vMLCont, vUp, 'nearest');
vAP(isnan(vAP)) = [];
vML(isnan(vML)) = [];

%%

% Remove infrequent locations
% Check if unique values are regularly spaced (e.g. in a grid). If they
% are, then find the mode of the grid spacing and round all positions
% to the nearest vertex.
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
mPos = [vAP' vML']; % all locations after removing infrequent ones
mPosUniq = unique(mPos, 'rows'); % only unique locations

% Create lookup table for converting between positions and matrix indices
vAPLookup = unique(mPos(:,1));
vMLLookup = unique(mPos(:,2));

%% Get continuous signal after each scan event

% Convert vUp to sample times of the continuous signal
vUpTime = round(interp1(vContTime, 1:length(vContTime), vUp, 'linear'));

% Create a cell matrix corresponding to all unique scan locations
cTrigEvents = cell(length(vAPLookup), length(vMLLookup));

% Create vectors with same length and times as continuous signal
vMLContSame = interp1(vGalvScanTime, vMLCont, vContTime, 'linear');
vAPContSame = interp1(vGalvScanTime, vAPCont, vContTime, 'linear');

% Iterate over all scan events
Spiky.main.SpikyWaitbar(0, 20);
nLen = length(vCont);
nPreStimSamples = round(p_nPreStimPeriod * nContFs);
nPostStimSamples = round(p_nPostStimPeriod * nContFs);

vStartR = vUpTime - nPreStimSamples;
vEndR  = vUpTime + nPostStimSamples;

% Loop over unique positions, and exctract all stimuli events for each.
vAPUptimes = vAPContSame(vUpTime);
vMLUptimes = vMLContSame(vUpTime);
cTrigEvents = cell(length(vAPLookup), length(vMLLookup));
for ap = 1:length(vAPLookup)
    Spiky.main.SpikyWaitbar((ap/length(vAPLookup))*20, 20);
    nAP = vAPLookup(ap);
    for ml = 1:length(vMLLookup)
        nML = vMLLookup(ml);
        
        % Extract all events at this position
        vIndx = vAPUptimes == nAP & vMLUptimes == nML;
        if isempty(vIndx), continue; end
        
        % Insert triggered signal into correct location cell
        for i = find(vIndx)
            cTrigEvents{ap, ml}(end+1, :) = vCont(:, vStartR(i):vEndR(i));
        end
    end
end
Spiky.main.SpikyWaitbar(20, 20);

%% Compute triggered averages at each locations

% Initialize figure
if bPlotETA
    hTrigAvgFig = figure('visible', 'off');
    centerfig(hFig, Spiky.main.GetGUIHandle())
    set(hTrigAvgFig, 'visible', 'on')
    drawnow
    Spiky.main.ThemeObject(hTrigAvgFig)
    hAx = [];
end

% Get relative time vector
vSigTime = (1:size(cTrigEvents{1,1}, 2)) ./ nContFs; % s
vSigTime = vSigTime - p_nPreStimPeriod;

% Create variables for holding average data
cTrigEventsAvg = cell(length(vAPLookup), length(vMLLookup));
mTrigAmpl = zeros(length(vAPLookup), length(vMLLookup));
mTrigAmplBaseline = zeros(length(vAPLookup), length(vMLLookup));

% Baseline indices (all indices before event trigger)
vPreTimeIndx = vSigTime <= 0;

% Signal measurement indices
vPostWinTimeIndx = (vSigTime <= p_nWinEnd) & (vSigTime >= p_nWinStart);

for i = 1:size(cTrigEvents, 1)
    for j = 1:size(cTrigEvents, 2)
        nInd = (i-1) * size(cTrigEvents, 2) + j;
        
        % Get all triggered events
        mSig = cTrigEvents{i, j};
        
        % Clean events by removing those with high baseline variance
        for s = 1:size(mSig, 1)
            vSig = mSig(s, :);
            vSTD(s) = std(vSig(vPreTimeIndx));
        end
        
        % Remove events with baseline standard deviation that exceed 50% of
        % the mean standard deviation
        vRemIndx = vSTD >= (mean(vSTD) * 1.5);
        mSig(vRemIndx, :) = [];
        
        % Average events
        vSig = mean(mSig);
        
        % Get baseline (average before trigger event)
        nBaseline = mean(vSig(vPreTimeIndx));
        
        % Subtract baseline
        vSig = vSig' - nBaseline;
        
        % Smooth event triggered average in time
        nSmoothFact = floor(p_nTempSmooth * nContFs / 1000); % indices
        if nSmoothFact > 0
            vSig = smooth(vSig, nSmoothFact);
        end
        vSigAmpl = detrend(vSig);
        
        % Compute evoked amplitude at each scan location
        mTrigAmpl(i, j) = mean(vSig(vPostWinTimeIndx)); % average

        % Compute baseline amplitude at each scan location
        mTrigAmplBaseline(i, j) = mean(vSig(vPreTimeIndx));

        % Assign
        cTrigEventsAvg{i, j} = vSig;
        
        % Plot event triggered average at current location
        if bPlotETA
            hAx(end+1) = subplot(size(cTrigEvents, 1), size(cTrigEvents, 2), nInd, 'parent', hTrigAvgFig);
            Spiky.main.ThemeObject(hAx(end))
            hold(hAx(end), 'on')
            hSig = plot(vSigTime, vSig, 'k-');
            Spiky.main.ThemeObject(hSig)
            plot([0 0], [-10 10], 'r--'); % vertical line at time zero
            hTxt = text(min(vSigTime), .1, sprintf('%d', size(mSig, 1)), 'parent', hAx(end), 'color', 'r', 'fontsize', 8);
            Spiky.main.ThemeObject(hTxt)
            axis(hAx, 'square')
            drawnow
        end
        
    end
end
mAvg = cell2mat(cTrigEventsAvg(:)');
vYLim = [min(mAvg(:)) max(mAvg(:))];

if bPlotETA
    set(hAx, 'ylim', vYLim, 'xlim', [min(vSigTime) max(vSigTime)], 'xticklabel', [], 'yticklabel', [])
    subplotsqueeze(hAx, 1.3)
    axis(hAx, 'off')
end

%%

% Subtract baseline and normalize to peak (absolute)
mImg = mTrigAmpl - mTrigAmplBaseline;

% Get image resolution (pixels/mm)
vImSizePx = max(mPos) - min(mPos); % image size, px
vPixMM = size(mImg) ./ vImSizePx; % pixels/mm ([height x width])

% Display intensity map
if ~ishandle(hFig), return; end
if gcf ~= hFig, figure(hFig); end
hImAx = subplot(2, 1, 1);
nSmoothFact = p_nSmooth * vPixMM(1);
if nSmoothFact > 0
    % Convert p_nSmooth from mm to pixels
    imagesc(filt2(mImg, nSmoothFact, 'lm'));
else
    imagesc(mImg);
end
vAxPos = get(hImAx, 'position');

% Estimate centroid (location of peak activation)
% Simple estimate is find pixels with vals > 2* std, threshold and find median XY of points
mThresh = zeros(size(mImg));
mThresh(mImg >= 2*std(mImg(:))) = 1;
[vY, vX] = find(mThresh);
nX = round(median(vX));
nY = round(median(vY));

% Plot 'centroid'
hold(hImAx, 'on')
plot(hImAx, nX, nY, 'wx', 'markersize', 12)

% Plot ROI
if isfield(FV, 'S1Map_GetROI')
    plot(hImAx, FV.S1Map_ROIx, FV.S1Map_ROIy, 'g.:')
end

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

% Set axis properties
vXInterp = -5:.5:5;
vXInterpVals = interp1(vXlabels, 1:length(vXlabels), vXInterp);
vYInterp = -5:.5:5;
vYInterpVals = interp1(vYlabels, 1:length(vYlabels), vYInterp);

axis(hImAx, 'image')
nMaxImgVal = max(abs(mImg(:)));
set(hImAx, 'tickdir', 'out', 'clim', [-nMaxImgVal nMaxImgVal], ...
    'xtick', vXInterpVals, 'xtickLabel', vXInterp, ...
    'ytick', vYInterpVals, 'ytickLabel', vYInterp, ...
    'ydir', 'normal')
Spiky.main.ThemeObject(hImAx)
grid(hImAx, 'on')

% Set title
hTit = title(sprintf('Response map'));
set(hTit, 'color', 'w', 'interpreter', 'none', 'fontsize', 8);

% Plot barrels on top of map (if available)
if isfield(FV.tData, 'BarrelOutlines')
    if isfield(FV.tData.BarrelOutlines, 'vXmm')
        vInternalInt = [];
        vBoundaryInt = [];
        for b = 1:length(FV.tData.BarrelOutlines)
            %axes(hImAx);
            hold(hImAx, 'on')
            vX = FV.tData.BarrelOutlines(b).vXmm;
            vY = FV.tData.BarrelOutlines(b).vYmm;
            % Contour
            vXi = interp1(vXlabels, 1:length(vXlabels), vY, 'linear', 'extrap');
            vYi = interp1(vYlabels, 1:length(vYlabels), vX, 'linear', 'extrap');
            
            % Emphasize contour if barrel is in cell's receptive field
            if ~g_bPlotS1MapQuiet
                plot(hImAx, vXi, vYi, 'w-', 'linewidth', 1)
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

% Check that figure is current and has not been closed
if ~ishandle(hFig), return; end
if gcf ~= hFig, figure(hFig); end

% Plot average evoked response at all locations
hAx = subplot(2, 1, 2);
hold(hAx, 'on')

vAvgSigAllPos = mean(mAvg');

% Draw thick bar to denote evoked activity window
vYY = [max(abs(vAvgSigAllPos))*-1 max(abs(vAvgSigAllPos))].*5;
hPatch = patch([p_nWinStart p_nWinEnd p_nWinEnd p_nWinStart p_nWinStart], [vYY(1) vYY(1) vYY(2) vYY(2) vYY(1)], [.7 0 0]);
set(hPatch, 'edgeColor', [.7 0 0])

% Plot error
vSigErr = std(mAvg') ./ sqrt(size(mAvg, 2));

% Plot all locations
axes(hAx);
Spiky.main.PlotMeanError(vAvgSigAllPos, vSigErr, [.6 .6 .6], vSigTime)

% Plot average
hold(hAx, 'on')
plot(hAx, [0 0], vYY, 'r--'); % vertical line at time zero

% Set axes properties
vYLim = [max(abs(vAvgSigAllPos))*-1 max(abs(vAvgSigAllPos))].*1.1;
set(hAx, 'xlim', [-p_nPreStimPeriod p_nPostStimPeriod], 'ylim', vYLim)
xlabel('Post-stimulus time (s)')

% Get unit of continuous signal
if isfield(FV.tData, [sCh '_Unit'])
    sUnit = FV.tData.([sCh '_Unit']);
else
    sUnit = 'V'; % default unit is Volts
end
ylabel(sUnit)

% Create colorbar
hCol = colorbar('peer', hImAx);
vPos = get(hImAx, 'position');
set(hCol, 'Position', [.85 vPos(2) .04 vPos(4)], 'tickdir', 'in');
Spiky.main.ThemeObject(hCol)
ylabel(hCol, sUnit)

Spiky.main.ThemeObject(hAx)
colormap(hImAx, fireice(1024))

axis(hAx, 'square')

hTit = title('Average evoked response');
Spiky.main.ThemeObject(hTit)

hHeader = header(['Channel ' sCh ' - GalvoScanner Map'], 10);
hFooter = footer(fullfile(FV.sDirectory, FV.sLoadedTrial), 10);
set([hHeader hFooter], 'color', 'w', 'interpreter', 'none')


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
