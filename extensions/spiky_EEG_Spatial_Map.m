function FV = spiky_EEG_Spatial_Map(FV)
% Generate spatial map of EEG measurements from point sources
%
% TODO Show anatomical coordinates on axes in millimeters 
%      Display when STIM is on in ERP mode
%

global Spiky;

% Hard-coded variables
nDecimateFactor = 1.5; % higher means more or original signal preserved
bInterp = 1; % interpolate
bMask = 1;

persistent p_nLowPass p_nHighPass p_bComputeERP p_nERPPreT p_nERPPostT p_csStimChannel p_bShowNames p_cNormMethod p_sRefCh p_sFidPath p_nSVDModes p_nResMult

% Default parameters
if isempty(p_nLowPass), p_nLowPass = 50; end
if isempty(p_nHighPass), p_nHighPass = 0.1; end
if isempty(p_cNormMethod), p_cNormMethod = 'z'; end % normalization method (z, percent or mean)
if isempty(p_bComputeERP), p_bComputeERP = 1; end % compute ERP, 0/1 = no/yes
if isempty(p_nERPPreT), p_nERPPreT = 0.2; end % s
if isempty(p_nERPPostT), p_nERPPostT = 1; end % s
if isempty(p_bShowNames), p_bShowNames = 1; end % show electrode names, 0/1 = no/yes
if isempty(p_sRefCh), p_sRefCh = 'all'; end % reference channel
if isempty(p_nSVDModes), p_nSVDModes = 0; end % SVD modes to denoise with (< 1 = no denoising)
if isempty(p_nResMult), p_nResMult = 1; end % resolution multiplier

cPrompt = {'Low pass (Hz)', ...
    'High pass (Hz)', ...
    'Normalization method (z, percent, mean)', ...
    'Compute ERP (1=yes)', ...
    'ERP baseline period (s)', ...
    'ERP post-stimulus period (s)', ...
    'Show electrode names (1=yes)', ...
    'Reference channel (string; ''all'')', ... };
    'SVD denoising modes', ...
    'Resolution multiplier' };
cAns = inputdlg(cPrompt, 'EEG Spatial Map', 1, { ...
    num2str(p_nLowPass), num2str(p_nHighPass), ...
    p_cNormMethod, num2str(p_bComputeERP), ...
    num2str(p_nERPPreT), num2str(p_nERPPostT), ...
    num2str(p_bShowNames), p_sRefCh, num2str(p_nSVDModes), num2str(p_nResMult) } );
if isempty(cAns), return, end

p_nLowPass = max([0 str2num(cAns{1})]);
p_nHighPass = max([0 str2num(cAns{2})]);
p_cNormMethod = cAns{3};
p_bComputeERP = max([0 str2num(cAns{4})]);
p_nERPPreT = max([0 str2num(cAns{5})]);
p_nERPPostT = max([0 str2num(cAns{6})]);
p_bShowNames = max([0 str2num(cAns{7})]);
p_sRefCh = cAns{8};
p_nSVDModes = str2num(cAns{9});
p_nResMult = str2num(cAns{10});

%% Load spatial coordinates of electrodes
% Get fiducials image
[sFile, p_sFidPath] = uigetfile({'*.mat';'*.*'}, 'Select electrode coordinate file', p_sFidPath);
if sFile == 0, return; end
load(fullfile(p_sFidPath, sFile), 'tData')

hWaitbar = waitbar(0, '');
centerfig(hWaitbar, Spiky.main.GetGUIHandle())

%% Initialize coordinates matrix
waitbar(0.1, hWaitbar, 'Initialize coordinate system');
csElecCoords = fieldnames(tData);
%iCh = regexpi(csElecCoords, '^elec_EGF_\d*$'); % Axona
iCh = regexpi(csElecCoords, '^elec_CH\d*$'); % Open Ephys
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
waitbar(0.2, hWaitbar);
mElecCoordsRel = mElecCoords - repmat(tData.refs_bregma, size(mElecCoords, 1), 1);
mElecCoordsRelMM = mElecCoordsRel .* tData.calibration_mmpix;

% Get skull contours (in mm and pixel coordinates)
mSkull = tData.cont_skull;
mSkullRel = mSkull - repmat(tData.refs_bregma, size(mSkull, 1), 1);
mSkullRelMM = mSkullRel .* tData.calibration_mmpix;

% Get matrix size to nearest 100 um
nWidthMM = max(abs(mSkullRelMM(:, 1))) * 2;
nHeightMM = max(abs(mSkullRelMM(:, 2)))  * 2;
nWidthMM = round(nWidthMM * p_nResMult) / p_nResMult;
nHeightMM = round(nHeightMM * p_nResMult) / p_nResMult;

nWidth = nWidthMM * p_nResMult + 1;
nHeight = nHeightMM * p_nResMult + 1;

% Bregma coordinates (pixels; in the middle)
vBregmaPos = round([nWidth / 2 nHeight / 2]);
mSkullRelPix = ceil(mSkullRelMM * p_nResMult) + repmat(vBregmaPos, size(mSkullRelMM, 1), 1);

% Electrode positions in pixel coordinates
mElecCoordsPix = ceil(mElecCoordsRelMM * p_nResMult) + repmat(vBregmaPos, size(mElecCoordsRelMM, 1), 1);

%% Get list of electrodes
csCh = fieldnames(FV.tData);
%iCh = regexpi(csCh, '^EGF_\d*$'); % Axona
iCh = regexpi(csCh, '^CH\d*$'); % Open Ephys
for i = fliplr(1:length(iCh))
    if isempty(iCh{i}), csCh(i) = []; end
end
csCh = sort(csCh);

%% Get the start and end points of the timeseries that will be analyzes
% This defaults to the current displayed range in the Spiky UI. If the
% UI is closed, not available etc, then the range will default to the
% entire timeseries range.
try
    nTimeBegin = FV.tData.([csCh{1} '_TimeBegin']);
    vSig = FV.tData.(csCh{1});
    vTime = (1:length(vSig)) ./ (FV.tData.([csCh{1} '_KHz']) * 1000);
    vXLim = Spiky.main.GetXLim() - nTimeBegin;
    [~, nStart] = min(abs(vTime - vXLim(1))); % samples
    [~, nEnd] = min(abs(vTime - vXLim(2))); % samples
catch
    nStart = 1;
    nEnd = length(FV.tData.(csCh{1}));
end
vRange = nStart:nEnd;

%% Assign electrode signals to correct location in mMap
mMap = [];
for ch = 1:length(csCh)
    waitbar(ch/length(csCh), hWaitbar, 'Computing signal matrix...');
    vPos = mElecCoordsPix(strcmp(['elec_' csCh{ch}], csElecCoords), :);
    if isempty(vPos)
        disp(['elec_' csCh{ch} ' is missing!'])
        continue
    end
    vSig = FV.tData.(csCh{ch});
    vSig = vSig(vRange);

    % Normalize signal to microvolts
    nGain = 1;
    if isfield(FV.tGain, csCh{ch})
        nGain = FV.tGain.(csCh{ch});
    end
    vSig = vSig ./ nGain .* 1000000;

    % Bandpass and decimate
    nFs = FV.tData.([csCh{ch} '_KHz']) * 1000;
    vTime = (1:length(vSig)) ./ nFs;
    [vSig, vTime, nNewFs] = Spiky.main.FilterChannel(vSig, vTime, nFs, p_nLowPass, p_nHighPass, 0, 'decimate', nDecimateFactor);
    nTimeBegin = FV.tData.([csCh{ch} '_TimeBegin']);

    % Time course normalization
    switch p_cNormMethod
        case 'z'
            % Z standardization
            vSig = (vSig - mean(vSig)) ./ std(vSig);
        case 'percent'
            % Percent signal change
            vSig = (vSig ./ mean(vSig)) .* 100; % vSig now has mean of 100
        case 'mean'
            % Mean subtraction
            vSig = vSig - mean(vSig);
    end
    
    if isempty(mMap)
        mMap = nan(nWidth, nHeight, length(vSig));
    end
    mMap(vPos(1), vPos(2), :) = vSig;
end
waitbar(0.5, hWaitbar, 'Substituting zeros with NaNs...');
%mMap(mMap(:) == 0) = NaN; % huge mem drain: mMap(:) == 0 duplicates matrix

%% Compute reference signal
waitbar(0, hWaitbar, 'Computing reference signal...');
switch lower(p_sRefCh)
    case 'all'
        % Subtract average of all channels
        mSigs = [];
        for e = 1:size(mElecCoordsPix, 1)
            waitbar(e / size(mElecCoordsPix, 1), hWaitbar, 'Subtracting reference signal...');
            mSigs(:, end+1) = squeeze(mMap(mElecCoordsPix(e, 1), mElecCoordsPix(e, 2), :));
        end
        vRef = mean(mSigs, 2);
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
            vPos = mElecCoordsPix(strcmp(['elec_' csCh{ch}], csElecCoords), :);
            
            vRef = squeeze(mMap(mElecCoordsPix(iMatch, 1), mElecCoordsPix(iMatch, 2), :));
        end
end

%% Subtract reference signal from all channels
if ~isempty(vRef)
    waitbar(0, hWaitbar, 'Subtracting reference signal...');
    if size(mMap, 3) < 50000
        mRef = zeros(size(mMap));
        for j = 1:size(mElecCoordsPix, 1)
            waitbar(j/size(mElecCoordsPix, 1), hWaitbar, 'Subtracting reference signal...');
            mRef(mElecCoordsPix(j, 1), mElecCoordsPix(j, 2), :) = vRef;
        end
        mMap = mMap - mRef;
        clear mRef;
    else
        % Iterate over samples (slow) rather than risking running out of memory
        nLen = size(mMap, 3);
        for j = 1:nLen
            waitbar(j/nLen, hWaitbar, 'Subtracting reference signal...');
            mMap(:, :, j) = mMap(:, :, j) - vRef(j);
        end
    end
end
%%

%% Compute Event Related Potentials
vERPTime = [];
if p_bComputeERP
    p_csStimChannel = Spiky.main.SelectChannelNumber(FV.csDigitalChannels, 'Select ERP trigger channel', p_csStimChannel);
    vUp = FV.tData.([p_csStimChannel '_Up']);
    nFs = FV.tData.([p_csStimChannel '_KHz']);
    vStimT = vUp - nTimeBegin;
    
    nERPPreSamp = round(p_nERPPreT * nNewFs);
    nERPPostSamp = round(p_nERPPostT * nNewFs);

    waitbar(0, hWaitbar, 'Computing event related potentials...');
    
    mMapAvg = [];
    nN = 0;
    for s = 1:length(vUp)
        waitbar(s / length(vUp), hWaitbar);
        [~, iMin] = min(abs(vTime - vStimT(s)));
        nStart = iMin - nERPPreSamp; % samples
        nEnd = iMin + nERPPostSamp; % samples
        if nStart < 1 || nEnd > length(vTime)
            continue;
        end
        if isempty(mMapAvg)
            mMapAvg = mMap(:, :, nStart:nEnd);
        else
            mMapAvg = (mMapAvg + mMap(:, :, nStart:nEnd));
        end
        nN = nN + 1;
    end
    nStart = 1;
    nEnd = size(mMapAvg, 3);
    mMapAvg = mMapAvg ./ nN;
    
    % Plot all ERPs on all channels
    mERPSigs = zeros(size(mMapAvg, 3), length(csCh));
    for ch = 1:length(csCh)
        vPos = mElecCoordsPix(strcmp(['elec_' csCh{ch}], csElecCoords), :);
        mERPSigs(:, ch) = squeeze(mMapAvg(vPos(1), vPos(2), :));
    end
    hFig = figure;
    
    Spiky.main.ThemeObject(hFig)
    vERPTime = (1:(nERPPreSamp+nERPPostSamp) + 1)' - (nERPPreSamp + 1);
    vERPTime = repmat(vERPTime, 1, size(mERPSigs, 2)) ./ (nNewFs); % s
    hAx = axes;
    hLines = plot(hAx, vERPTime, mERPSigs);
    for h = 1:length(hLines)
        set(hLines(h), 'color', FV.mColors(h, :))
    end
    xlabel(hAx, 'Time (s)')
    ylabel(hAx, 'ERP (uV)')
    axis(hAx, 'tight')
    set(hAx, 'XGrid', 'on', 'YGrid', 'on')
    Spiky.main.ThemeObject(hAx)

    hLeg = legend(csCh);
    Spiky.main.ThemeObject(hLeg)
    
    hTit = title(sprintf('ERP, n = %d', nN));
    Spiky.main.ThemeObject(hTit)
    
    % Time marker
    hold(hAx, 'on')
    hMarker = plot(hAx, [nan nan], [floor(min(mERPSigs(:))) ceil(max(mERPSigs(:)))]);
    set(hMarker, 'linewidth', 2, 'tag', 'EEGSpatialMapERPMarker', 'color', 'r')

    mMap = mMapAvg;
    clear mMapAvg
end

%% Interpolate matrices
waitbar(0, hWaitbar, 'Interpolating matrices...');
if bInterp
    nFrames = size(mMap, 3);
    for i = 1:nFrames
        waitbar(i / nFrames, hWaitbar);
        mImg = mMap(:, :, i);
        mMap(:, :, i) = inpaint_nans(mImg, 2);
    end
end

%% Denoise frames by SVD
if p_nSVDModes > 0
    waitbar(0.5, hWaitbar, 'SVD denoising...');
    [mMap, ~] = svdmoviefilt(mMap, p_nSVDModes);
end

close(hWaitbar)

%% Play spatial map
hFig = figure;
Spiky.main.ThemeObject(hFig)
colormap(fireice)
hAx = axes;
hIm = imagesc(mMap(:, :, 1), 'parent', hAx, 'tag', 'EEGSpatialMapImg');
hold(hAx, 'on')
Spiky.main.ThemeObject(hAx)

mSkull = tData.cont_skull;
mSkullRel = mSkull - repmat(tData.refs_bregma, size(mSkull, 1), 1);

%mSkullRelPix
mSkullRelMM = mSkullRel .* tData.calibration_mmpix;
%mElecCoordsRelMM
%mElecCoordsPix

if bInterp
    hLabels = plot(mElecCoordsPix(:, 2), mElecCoordsPix(:, 1), 'wo');
end

% Text labels
if p_bShowNames
    csElec = char();
    for c = 1:length(csElecCoords)
        csElec(c, 1:length(csElecCoords{c})) = char(csElecCoords{c});
    end
    hTxt = text(mElecCoordsPix(:, 2) - 1.2, mElecCoordsPix(:, 1), csElec(:, 6:end));
    set(hTxt, 'FontSize', 8, 'interpreter', 'none', 'HorizontalAlignment', 'center', 'color', 'w');
end

hSkull = plot(mSkullRelPix(:, 2), mSkullRelPix(:, 1), 'w--');

%% Compute appropriate clim's from data
mMax = max(mMap, [], 3);
mMin = min(mMap, [], 3);
cCLim = [mean(mMin(:)) mean(mMax(:))] .* 2;
nMax = max(abs(cCLim));
hAx.CLim = [-nMax nMax];

% Compute clim from the average min/max

hCol = colorbar;
hCol.Label.String = 'uV';
axis(hAx, 'image')
hAx.View = [90 90];
hAx.XLim = [min(mSkullRelPix(:, 2))-.5 max(mSkullRelPix(:, 2))+.5];
hAx.YLim = [min(mSkullRelPix(:, 1))-.5 max(mSkullRelPix(:, 1))+.5];
hTit = title(hAx, 'Anterior');
if bMask
    mMask = poly2mask(mSkullRelPix(:, 2), mSkullRelPix(:, 1), nWidth, nHeight);
else
    mMask = [];
end
Spiky.main.ThemeObject(hTit, 'fontsize', 12)

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
tD = struct('mMap', mMap, 'i', 1, 'vTime', vTime, 'vRange', vRange, 'mMask', mMask, ...
    'bComputeERP', p_bComputeERP, 'vERPTime', vERPTime);

% Display initial frame
set(hFig, 'userdata', tD, 'KeyPressFcn', @UpdateFrame, 'tag', 'EEGSpatialMap')
UpdateFrame(1);


%%
return


function cmap = fireice(m)
% FIREICE
if ~nargin
    m = 64;
end
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


% AUX function for changing the displayed frame
function UpdateFrame(varargin)
% Input:
%   UpdateFrame(F) where F is a valid frame number, or handle of control
%   button.
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
    % Change framenumber depending on button pressed
    hBut = varargin{1};
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

% Validate frame number
iFrame = max([1 iFrame]);
iFrame = min([iFrame length(tD.vRange)]);

% Update time marker
if tD.bComputeERP
    % Update time markers in ERP window
    set(findobj('tag', 'EEGSpatialMapERPMarker'), ...
        'xdata', [tD.vERPTime(tD.vRange(iFrame)) tD.vERPTime(tD.vRange(iFrame))]);
else
    % Update marker in Spiky window
    Spiky.main.SetTimeMarker(tD.vTime(tD.vRange(iFrame)));
end

% Update image
mImg = tD.mMap(:, :, tD.vRange(iFrame));
if ~isempty(tD.mMask)
    mImg = mImg .* tD.mMask;
end
mImg(isnan(mImg)) = 0;

if tD.bComputeERP
    hStat = findobj('tag', 'EEGSpatialMapStatus');
    hStat.String = sprintf('%.3f s', tD.vERPTime(tD.vRange(iFrame)));
end

set(hIm, 'cdata', mImg)
drawnow

% Update figure userdata
tD.nCurrentFrame = iFrame;
set(hFig, 'userdata', tD);

return;

sChar = get(hFig, 'CurrentCharacter');
if ~isempty(sChar)
    switch double(sChar)
        case 113 % 'q'
            disp('quit...')
        case 31
            figure(hFig)
            set(hFig, 'CurrentCharacter', char(0))
        case 30
            iFrame = iFrame + 1;
        case 28
            iFrame = iFrame - 1;
            set(hFig, 'CurrentCharacter', char(0))
        case 29
            iFrame = iFrame + 1;
            set(hFig, 'CurrentCharacter', char(0))
    end
end

return
