function FV = spiky_EEG_Spatial_Map(DATA, varargin)
% Generate a spatial map of EEG/ECoG measurements from an array of point
% sources.
% 
% Usage:
%   spiky_EEG_Spatial_Map(FV)
%   spiky_EEG_Spatial_Map(ERP, FS, PRET, POSTT, COORDS, SKULL, CH, RES, COL)
%    where:  
%       ERP      matrix of ERPs, [time x channel]
%       FS       sample rate (Hz)
%       PRET     pre-stimulus delay (s)
%       POSTT    post-stimulus delay (s)
%       COORDS   matrix of [x y] electrode coordinates (mm; OPTIONAL)
%       SKULL    matrix of [x y] coordinates of skull boundary (mm; OPTIONAL)
%       CH       names of electrodes (cells with strings; OPTIONAL)
%       RES      interpolation resolution (mm)
%       COL      is the colormap to use (e.g. 'fireice', 'jet')
%       METHOD   interpolation method; see inpaint_nans (default 2)
%
%       Note: OPTIONAL parameters cannot be omitted! Provide the empty
%       vector [] instead if they are omitted.
%
% Sub-functions:
%   InterpolateERPFrame()
%   GetParametersUI()
%   GetSpatialCoords()
%   GetElectrodes()
%   GetTimeRange()
%   GetReferenceSignal()
%   GetERPs()
%   DisplayERPs()
%   GetTimeBegin()
%   GetElectrodeCoords()
%   GetData()
%   UpdateFrame()
%
global Spiky
persistent p_tParms

if ~isstruct(DATA)
    % Pre-process data if ERPs are passed directly to function
    % Required parameters if calling function from outside Spiky
    mERPs = DATA;
    nFs = varargin{1};
    p_tParms.nERPPreT = varargin{2};
    p_tParms.nERPPostT = varargin{3};
    mCoords = varargin{4};
    mSkullMM = varargin{5};
    csChDescr = varargin{6};
    p_tParms.nRes = varargin{7}; % interpolation pixel resolution (mm)
    sColMap = varargin{8}; % colormap
    p_tParms.bShowNames = 1;
    if length(varargin) == 9
        p_tParms.nInterpMethod = varargin{9};
    else
        p_tParms.nInterpMethod = 2;
    end
    
    % Check that sample rate makes sense
    nExpLen = (p_tParms.nERPPreT * nFs) + (p_tParms.nERPPostT * nFs);
    if size(mERPs, 1) ~= nExpLen
        error('Sample rate and/or pre/post delay is inconsistent with data length')
    end
    
    % Load spatial coordinates of electrodes and skull
    if isempty(mCoords)
        [mCoords, ~, mSkullMM] = GetSpatialCoords();
        if size(mCoords, 1) ~= size(mERPs, 2)
            error('Number of marked electrodes not the same as number of ERP channels')
        end
    end
    
else
    FV = DATA;
    sColMap = 'fireice';
    
    %% Default parameters
    if isempty(p_tParms)
        p_tParms.cNormMethod = 'z'; % normalization method (z, percent or mean)
        p_tParms.bComputeERP = 1; % compute ERP, 0/1 = no/yes
        p_tParms.nERPPreT = 0.05; % s
        p_tParms.nERPPostT = 0.4; % s
        p_tParms.bShowNames = 1; % show electrode names, 0/1 = no/yes
        p_tParms.sRefCh = 'all'; % reference channel
        p_tParms.nSVDModes = 0; % SVD modes to denoise with (< 1 = no denoising)
        p_tParms.nRes = 0.1; % interpolation pixel resolution (mm)
        p_tParms.nInterpMethod = 2; % interpolation method
    end

    % Get parameters from UI
    p_tParms = GetParametersUI(p_tParms);
    
    % Load spatial coordinates of electrodes and skull
    [mElecMM, csElecCoords, mSkullMM] = GetSpatialCoords();
    
    % Get electrode names
    [csCh, csChDescr] = GetElectrodes(FV);
    
    % Get time-range to analyze
    vTimeRange = GetTimeRange(FV);
    
    % Get electrode coordinates
    mCoords = GetElectrodeCoords(csCh, csChDescr, csElecCoords, mElecMM);
    
    % Get data
    [mChData, vTime, nFs] = GetData(FV, csCh, p_tParms, vTimeRange);

    % Subtract reference signal
    vRef = GetReferenceSignal(mChData, p_tParms);
    if ~isempty(vRef)
        mChData = mChData - repmat(vRef, 1, size(mChData, 2));
    end

    % Compute Event Related Potentials
    if p_tParms.bComputeERP
        mERPs = GetERPs(FV, p_tParms, mChData, vTime, nFs);
    end

    %% Denoise frames by SVD
    % TODO Move into display function
    %if p_nSVDModes > 0
    %    waitbar(0.5, hWaitbar, 'SVD denoising...');
    %    [mMap, ~] = svdmoviefilt(mMap, p_nSVDModes);
    %end
end

% Plot ERPs
vERPTime = DisplayERPs(mERPs, p_tParms, nFs, csChDescr);

%% Save data (used for debugging direct analysis option)
if 1
    nERPPreT = p_tParms.nERPPreT;
    nERPPostT = p_tParms.nERPPostT;
    nRes = p_tParms.nRes;
    save('eeg_spatial_map_vars.mat', 'mERPs', 'nFs', 'nERPPreT', 'nERPPostT', 'mCoords', 'mSkullMM', 'csChDescr', 'nRes', 'sColMap');
    % Run with;
    %load('eeg_spatial_map_vars.mat')
    %spiky_EEG_Spatial_Map(mERPs, nFs, nERPPreT, nERPPostT, mCoords, mSkullMM, csChDescr, nRes, sColMap)
end

%% Play spatial map
hFig = figure;
hAx = axes;

[mImg, vX_lookup, vY_lookup] = InterpolateERPFrame(mERPs(1, :), mCoords, p_tParms.nRes, p_tParms.nInterpMethod);
hIm = imagesc(vY_lookup, vX_lookup, mImg, 'parent', hAx, 'tag', 'EEGSpatialMapImg');
hold(hAx, 'on')

if ~isempty(mSkullMM)
    hSkull = plot(hAx, mSkullMM(:, 2), mSkullMM(:, 1), 'w--');
end
hLabels = plot(hAx, mCoords(:, 2), mCoords(:, 1), 'go');

% Text labels
if p_tParms.bShowNames && ~isempty(csChDescr)
    hTxt = text(mCoords(:, 2) - 1.2, mCoords(:, 1), csChDescr, 'parent', hAx);
    set(hTxt, 'FontSize', 8, 'interpreter', 'none', 'HorizontalAlignment', 'center', 'color', 'w');
end

% Compute appropriate clim's from data
mMax = max(mERPs, [], 1);
mMin = min(mERPs, [], 1);
nStd = std(mERPs(:));
mMin(abs(mMin) < nStd*2) = [];
mMax(abs(mMax) < nStd*2) = [];
cCLim = [mean(mMin(:)) mean(mMax(:))] .* 2.5;
nMax = min(abs(cCLim));
hAx.CLim = [-nMax nMax];

hCol = colorbar(hAx);
axes(hAx)
axis image
hAx.View = [-90 -90];
if ~isempty(mSkullMM)
    hAx.XLim = [min(mSkullMM(:, 2))-.5 max(mSkullMM(:, 2))+.5];
    hAx.YLim = [min(mSkullMM(:, 1))-.5 max(mSkullMM(:, 1))+.5];
end
hTit = title(hAx, 'Anterior');

if ~isempty(Spiky)
    Spiky.main.ThemeObject([hFig hAx hCol hTit])
end
colormap(hFig, eval(sprintf('%s(1024)', sColMap)))

% Navigation buttons
set(hFig, 'units', 'normalized')
uicontrol(hFig, 'Style', 'togglebutton', 'units', 'normalized', 'Position', [0 0 .05 .05], ...
    'Callback', @UpdateFrame, 'String', '<');
uicontrol(hFig, 'Style', 'togglebutton', 'units', 'normalized', 'Position', [0.05 0 .05 .05], ...
    'Callback', @UpdateFrame, 'String', '>');
uicontrol(hFig, 'Style', 'pushbutton', 'units', 'normalized', 'Position', [.1 0 .05 .05], ...
    'Callback', @UpdateFrame, 'String', '<<');
uicontrol(hFig, 'Style', 'pushbutton', 'units', 'normalized', 'Position', [.15 0 .05 .05], ...
    'Callback', @UpdateFrame, 'String', '>>');
uicontrol(hFig, 'Style', 'edit', 'units', 'normalized', 'Position', [.2 0 .8 .05], 'tag', 'EEGSpatialMapStatus');

% Build userdata structure
tD = struct('mERP', mERPs, 'i', 1, 'mSkullMM', mSkullMM, ...
    'nRes', p_tParms.nRes, 'vERPTime', vERPTime, 'mCoords', mCoords, ...
    'nInterpMethod', p_tParms.nInterpMethod );

% Display initial frame
set(hFig, 'userdata', tD, 'KeyPressFcn', @UpdateFrame, 'tag', 'EEGSpatialMap')
UpdateFrame(1);

%%

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateFrame(varargin)
% AUX function for changing the displayed frame
%
% Input:
%   UpdateFrame(F) where F is a valid frame number, or handle of control
%   button.
%
% TODO: Split function into subfunctions; e.g GetFrameNumber()
%
%

global Spiky

hFig = findobj('tag', 'EEGSpatialMap');
tD = get(hFig, 'userdata');
hIm = findobj(hFig, 'tag', 'EEGSpatialMapImg');

% Get frame number to display
if isfield(tD, 'nCurrentFrame')
    iFrame = tD.nCurrentFrame;
else
    iFrame = 1; % default frame to display is first
end

if isnumeric(varargin{1})
    % Frame number is provided
    iFrame = varargin{1};
elseif isobject(varargin{1})
    if strcmp(get(varargin{1}, 'type'), 'axes')
        % ERP axis was clicked
        mCurrPnt = get(varargin{1}, 'CurrentPoint');
        vX = mCurrPnt(1, 1);
        hLines = findobj(varargin{1}, 'type', 'line');
        vXX = get(hLines(end), 'xdata');
        [~, iFrame] = min(abs(vXX - vX));
    elseif strcmp(get(varargin{1}, 'type'), 'uicontrol')
        % Change framenumber depending on button pressed
        hBut = varargin{1};
        if ~ishandle(hBut), return; end
        switch hBut.String
            case '>'
                % Play continuously forward, until last frame
                while 1
                    UpdateFrame(iFrame);
                    iFrame = iFrame + 1;
                    if hBut.Value == 0, break; end
                end
            case '<'
                % Play continuously bachwards, until first frame
                while 1
                    UpdateFrame(iFrame);
                    iFrame = iFrame - 1;
                    if hBut.Value == 0, break; end
                end
            case '>>'
                % Increment frame
                iFrame = iFrame + 1;
            case '<<'
                % Decrement frame
                iFrame = iFrame - 1;
        end
    end
end

% Validate frame number
iFrame = max([1 iFrame]);
iFrame = min([iFrame size(tD.mERP, 1)]);

% Update time markers in ERP window
set(findobj('tag', 'EEGSpatialMapERPMarker'), ...
    'xdata', [tD.vERPTime(iFrame) tD.vERPTime(iFrame)]);

%% Compute frame
[mImg, vX_lookup, vY_lookup] = InterpolateERPFrame(tD.mERP(iFrame, :), tD.mCoords, tD.nRes, tD.nInterpMethod);

%% Compute mask
if ~isfield(tD, 'mMask')
    % Convert mm to pixels
    vX = tD.mSkullMM(:, 1);
    vY = tD.mSkullMM(:, 2);
    vYpx = interp1(vX_lookup, 1:length(vX_lookup), vX, 'linear', 'extrap');
    vXpx = interp1(vY_lookup, 1:length(vY_lookup), vY, 'linear', 'extrap');
    mMask = poly2mask(vXpx, vYpx, size(mImg, 1), size(mImg, 2));
    tD.mMask = mMask;
end

% Mask image
mImg = mImg .* tD.mMask;
mImg(isnan(mImg)) = 0;
mImg(mImg == 0) = nan;

hStat = findobj('tag', 'EEGSpatialMapStatus');
hStat(1).String = sprintf('%.3f s', tD.vERPTime(iFrame));

set(hIm, 'cdata', mImg)
drawnow

% Update figure userdata
tD.nCurrentFrame = iFrame;
if ishandle(hFig)
    set(hFig, 'userdata', tD);
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mERP_i, vX_lookup, vY_lookup] = InterpolateERPFrame(mERP, mCoords, nRes, nInterpMethod)
% Interpolate spatial ERP signal
%
% Usage: InterpolateERPFrame(D, C, R), where
%   D is a time X channel data matrix
%   C is a matrix of the spatial coordinates of electrodes (mm)
%   R is the pixel resolution of the interpolated signal (mm)
%
nSlack = 5;
mRange = [min(mCoords) - nSlack; max(mCoords) + nSlack];
vX_lookup = mRange(1, 1):nRes:mRange(2, 1);
vY_lookup = mRange(1, 2):nRes:mRange(2, 2);
mERP_i = nan(length(vX_lookup), length(vY_lookup));
for i = 1:length(mERP)
    vCoords = mCoords(i, :);
    [~, iMinX] = min(abs(vCoords(1) - vX_lookup));
    [~, iMinY] = min(abs(vCoords(2) - vY_lookup));
    mERP_i(iMinX, iMinY) = mERP(i);
end
mERP_i = inpaint_nans(mERP_i, nInterpMethod);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tParms = GetParametersUI(tParms)
% Open UI to enter analysis parameters
%
cPrompt = {'Normalization method (z, percent, mean)', ...
    'Compute ERP (1=yes)', ...
    'ERP baseline period (s)', ...
    'ERP post-stimulus period (s)', ...
    'Show electrode names (1=yes)', ...
    'Reference channel (string; ''all'')', ... };
    'SVD denoising modes', ...
    'Interpolation resolution (mm)', ...
    'Interpolation method (1-5)' };
cAns = inputdlg(cPrompt, 'EEG Spatial Map', 1, { ...
    tParms.cNormMethod, num2str(tParms.bComputeERP), ...
    num2str(tParms.nERPPreT), num2str(tParms.nERPPostT), ...
    num2str(tParms.bShowNames), tParms.sRefCh, ...
    num2str(tParms.nSVDModes), num2str(tParms.nRes), ...
    num2str(tParms.nInterMethod) } );
if isempty(cAns), return, end
tParms.cNormMethod = cAns{1};
tParms.bComputeERP = max([0 str2num(cAns{2})]);
tParms.nERPPreT = max([0 str2num(cAns{3})]);
tParms.nERPPostT = max([0 str2num(cAns{4})]);
tParms.bShowNames = max([0 str2num(cAns{5})]);
tParms.sRefCh = cAns{6};
tParms.nSVDModes = str2num(cAns{7});
tParms.nRes = str2num(cAns{8});
tParms.nInterMethod = str2num(cAns{9});
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mElecMM, csElecCoords, mSkullMM] = GetSpatialCoords(varargin)
% Load spatial coordinates of electrodes from disk
%
% Usage:
%   GetSpatialCoords(FILE)  Load known file
%   GetSpatialCoords()      Select a file manually
%
if nargin == 1
    sFile = varargin{1};
    load(sFile, 'tData')
else
    persistent p_sFidPath
    [sFile, p_sFidPath] = uigetfile({'*.mat';'*.*'}, 'Select electrode coordinate file', p_sFidPath);
    if sFile == 0, return; end
    load(fullfile(p_sFidPath, sFile), 'tData')
end
csElecCoords = fieldnames(tData);
iCh = regexpi(csElecCoords, '^elec_CH\d*$');
if isempty(cell2mat(iCh))
    iCh = regexpi(csElecCoords, '^elec_\d*$');
end
if isempty(cell2mat(iCh))
    warndlg('No matching electrodes: Check labels created in "Get_Fiducials".')
    return;
end
vRem = [];
for cs = 1:length(iCh)
    if isempty(iCh{cs})
        vRem(end + 1) = cs;
    end
end
csElecCoords(vRem) = [];
csElecCoords = sort(csElecCoords);
mElecCoords = [];
for cs = 1:length(csElecCoords)
    mElecCoords(end+1, :) = tData.(csElecCoords{cs});
end
%% Normalize coordinates to mm relative to bregma
mElecCoordsRel = mElecCoords - repmat(tData.refs_bregma, size(mElecCoords, 1), 1);
mElecMM = mElecCoordsRel .* tData.calibration_mmpix;

% Get skull contours (in mm and pixel coordinates)
mSkull = tData.cont_skull;
mSkullRel = mSkull - repmat(tData.refs_bregma, size(mSkull, 1), 1);
mSkullMM = mSkullRel .* tData.calibration_mmpix;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cmap = fireice(m)
% FIREICE colormap
if ~nargin
    m = 64;
end
clrs = [0.75 1 1; 0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0; 1 1 0.75];
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [csCh, csChDescr] = GetElectrodes(FV)
% Get list of electrodes
%
csCh = fieldnames(FV.tData);
iCh = regexpi(csCh, '^CH\d*$');

for i = fliplr(1:length(iCh))
    if isempty(iCh{i}), csCh(i) = []; end
end
csCh = sort(csCh);

% Get channel descriptions
csChDescr = csCh;
for c = 1:length(csCh)
    if isfield(FV, 'tChannelDescriptions')
        nIndx = strcmp(csCh{c}, {FV.tChannelDescriptions.sChannel});
        if ~isempty(FV.tChannelDescriptions(nIndx))
            csChDescr{c} = FV.tChannelDescriptions(nIndx).sDescription;
        end
    end
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vTimeRange = GetTimeRange(FV)
% Get the start and end points of the timeseries that will be analyzed
% This defaults to the current displayed range in the Spiky UI. If the
% UI is closed, not available etc, then the range will default to the
% entire timeseries range.
%
global Spiky
try
    [nTimeBegin, sCh] = GetTimeBegin(FV);
    vSig = FV.tData.(sCh);
    vTime = (1:length(vSig)) ./ (FV.tData.([sCh '_KHz']) * 1000);
    vXLim = Spiky.main.GetXLim() - nTimeBegin;
    [~, nStart] = min(abs(vTime - vXLim(1))); % samples
    [~, nEnd] = min(abs(vTime - vXLim(2))); % samples
catch
    nStart = 1;
    nEnd = length(FV.tData.(sCh));
end
vTimeRange = nStart:nEnd;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nTimeBegin, sCh] = GetTimeBegin(FV)
% Get begin time (first timepoint)
%
csFields = fieldnames(FV.tData);
iTB = find(~cellfun(@isempty, strfind(csFields', 'TimeBegin')));
nTimeBegin = FV.tData.([csFields{iTB(1)}]);
sCh = csFields{iTB(1)}(1:end-10);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vRef = GetReferenceSignal(mChData, p_tParms)
% Compute reference signal
%
switch lower(p_tParms.sRefCh)
    case 'all'
        % Subtract average of all channels
        vRef = nanmean(mChData, 2);
    case 'none'
        vRef = [];
    otherwise
        % Subtract a reference channel from all other channels
        vRef = [];

        iMatch = find(strcmpi(['elec_' p_sRefCh], csElecCoords));
        if isempty(iMatch)
            Spiky.main.sp_disp(['Unknown name for reference channel elec_' p_sRefCh]);
        else
            % Subtract a reference channel from all other channels
            vRef = mChData(:, iMatch); % TODO Needs to be tested
        end
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vERPTime = DisplayERPs(mERPs, p_tParms, nFs, csCh)
% Display event related potentials as traces
%
global Spiky
hFig = figure;

nERPPreSamp = round(p_tParms.nERPPreT * nFs);
nERPPostSamp = round(p_tParms.nERPPostT * nFs);

%vERPTime = (1:(nERPPreSamp+nERPPostSamp) + 1)' - (nERPPreSamp + 1);
%vERPTime = repmat(vERPTime, 1, size(mERPs, 2)) ./ (nFs); % s
vERPTime = linspace(-p_tParms.nERPPreT, p_tParms.nERPPostT, size(mERPs, 1))';
vERPTime = repmat(vERPTime, 1, size(mERPs, 2));

hAx = axes;
hLines = plot(hAx, vERPTime, mERPs);
xlabel(hAx, 'Time (s)')
ylabel(hAx, ['Event Related Potential'])
axis(hAx, 'tight')
set(hAx, 'XGrid', 'on', 'YGrid', 'on')
hLeg = legend(csCh);
hTit = title('ERP');

if ~isempty(Spiky)
    Spiky.main.ThemeObject([hFig hAx hLeg hTit])
end

% Time marker
hold(hAx, 'on')
hMarker = plot(hAx, [nan nan], [floor(min(mERPs(:))) ceil(max(mERPs(:)))]);
set(hMarker, 'linewidth', 2, 'tag', 'EEGSpatialMapERPMarker', 'color', 'r')

% Enable axis click trigger
set(hAx, 'ButtonDownFcn', @UpdateFrame)

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mERPs = GetERPs(FV, p_tParms, mChData, vTime, nFs)
% Compute event related potentials
%
%
global Spiky
if isempty(Spiky)
    error('Spiky not available')
end
persistent p_csStimChannel
p_csStimChannel = Spiky.main.SelectChannelNumber(FV.csDigitalChannels, 'Select ERP trigger channel', p_csStimChannel);

% Get trigger events in the displayed time range
[vUp, ~] = Spiky.main.GetEventPairs(p_csStimChannel, 'range');

[nTimeBegin, ~] = GetTimeBegin(FV);

vStimT = vUp - nTimeBegin;
nERPPreSamp = round(p_tParms.nERPPreT * nFs);
nERPPostSamp = round(p_tParms.nERPPostT * nFs);

mERPs = [];
nN = 0;
for s = 1:length(vUp)
    % Check that pre/post pulse window is within displayed time range
    [~, iMin] = min(abs(vTime - vStimT(s)));
    nStart = iMin - nERPPreSamp; % samples
    nEnd = iMin + nERPPostSamp; % samples
    if nStart < 1 || nEnd > length(vTime)
        continue;
    end
    
    if isempty(mERPs)
        mERPs = mChData(nStart:nEnd, :);
    else
        mERPs = mERPs + mChData(nStart:nEnd, :);
    end
    nN = nN + 1;
end
if nN == 0
    warndlg('No trigger pulses detected: Try expanding time-range in Spiky UI.')
    return;
elseif nN <= 5
    Spiky.main.sp_disp('Fewer than 5 trigger pulses found: Try expanding time-range in Spiky UI.')
end
nStart = 1;
nEnd = size(mERPs, 1);
mERPs = mERPs ./ nN;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mCoords = GetElectrodeCoords(csCh, csChDescr, csElecCoords, mElecMM)
% Get coordinates of electrodes (returns a subset of input, in millimeters)
%
for ch = 1:length(csCh)
    % Get electrode index (based on channel name first, description second)
    iCh = strcmp(['elec_' csCh{ch}], csElecCoords);
    if ~any(iCh)
        iCh = strcmp(['elec_' csChDescr{ch}], csElecCoords);
    end
    if ~any(iCh)
        disp(['elec_' csCh{ch} ' is missing!'])
        continue
    end
    mCoords(ch, :) = mElecMM(iCh, :);
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mChData, vTime, nFs] = GetData(FV, csCh, p_tParms, vTimeRange)
% Get data
%
global Spiky
hWaitbar = waitbar(0, '');
mChData = [];
for ch = 1:length(csCh)
    waitbar(ch/length(csCh), hWaitbar, 'Computing signal matrix...');

    vSig = FV.tData.(csCh{ch});
    vSig = vSig(vTimeRange);

    % Apply filters
    nFs = FV.tData.([csCh{ch} '_KHz']) * 1000;
    vTime = (1:length(vSig)) ./ nFs;
    if ~isempty(Spiky)
        [vSig, vTime, nFs] = Spiky.main.GetFilteredChannel(csCh{ch}, vSig, vTime, 'decimate');
    end
    
    % Time course normalization
    switch p_tParms.cNormMethod
        case 'z'
            % Z standardization
            vSig = (vSig - mean(vSig)) ./ std(vSig);
            sUnit = 'Z';
        case 'percent'
            % Percent signal change
            vSig = (vSig ./ mean(vSig)) .* 100; % vSig now has mean of 100
            sUnit = '%';
        case 'mean'
            % Mean subtraction
            vSig = vSig - mean(vSig);
    end
    mChData(:, ch) = vSig;
end
close(hWaitbar)
return