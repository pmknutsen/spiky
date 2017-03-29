function FV = spiky_EEG_Source_Analysis(FV)
% Set fidicual coordinates of EEG grid in an interactive UI.
%
% TODO:
%   Auto-rotate image according to bregma/lambda axis
%
%

global Spiky;

% Get fiducials image
[sFile, sPath] = uigetfile({'*.png';'*.jpg';'*.jpeg';'*.*'}, 'Select fiducials image');
if sFile == 0, return; end

% Initialize figure and UI elements
hFig = figure;
set(hFig, 'position', [100 100 800 600]);
centerfig(hFig, Spiky.main.GetGUIHandle())

hImgPanel = uipanel('parent', hFig, ...
    'backgroundcolor', 'k', ...
    'position', [0.2 0 0.8 1] );
hAx = axes('parent', hImgPanel, 'position', [0 0 1 1], 'tag', 'EEGGetFiducials');

%% Display buttons
hButPanel = uipanel('parent', hFig, 'position', [0 0 0.2 1] );
DrawButton(hButPanel, 'Mid-line', [0 1-.1 1 .1], @SetMidline)
DrawButton(hButPanel, 'References', [0 1-.2 1 .1], @SetReferences)
DrawButton(hButPanel, 'Electrodes', [0 1-.3 1 .1], @SetElectrodes)
DrawButton(hButPanel, 'Calibrate', [0 1-.4 1 .1], @SetCalibration)
DrawButton(hButPanel, 'Contours', [0 1-.5 1 .1], @SetContour)
DrawButton(hButPanel, 'Show/Hide Image', [0 1-.6 1 .1], @ShowHideImage)
DrawButton(hButPanel, 'Save', [0 1-.7 1 .1], @Export)
DrawButton(hButPanel, 'Exit', [0 1-.8 1 .1], 'close(gcf)')

% Message box
uicontrol('Parent', hButPanel, ...
    'string', 'Step-through buttons above, starting from the top. You can go back and make changes later.', 'unit', 'normalized', ...
    'style', 'text', ...
    'position', [0 0 1 .15], ...
    'tag', 'EEGGetFiducialsMsg' );

% Display image and data (if exists)
SetData('imagepath', sPath);
SetData('imagefile', sFile);
UpdateImage(1);
tData = get(GetAxis, 'userdata');
LoadData();
UpdateImage(1); % update in case we need to rotate

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateImage(bAllowRotation)
% Update image (rotate if bregma/lambda line has been set)
tData = get(GetAxis, 'userdata');
mImg = imread(fullfile(tData.imagepath, tData.imagefile));

% Rotate image
if bAllowRotation
    if isfield(tData, 'refs_midline')
        % Get angle
        vX = tData.refs_midline(:, 1);
        vY = tData.refs_midline(:, 2);
        nRad = atan2(diff(vX), diff(vY));
        mImg = imrotate(mImg, rad2deg(-nRad), 'bicubic', 'crop');
        SetData('refs_rotation', nRad)
    end
else
    SetData('refs_rotation', [])
end

% Update image
if ~isempty(GetImage)
    set(GetImage, 'cdata', mImg);
else
    hImg = image(mImg, 'parent', GetAxis, 'tag', 'EEGGetFiducialsImg');
end

% For some reason, image() above removes the old axis. Hence, we put back
% tData into UserData here.
set(GetAxis, 'userdata', tData);
axis(GetAxis, 'image')
axis(GetAxis, 'off')
hold(GetAxis, 'on')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawButton(hParent, sStr, vPos, sCallback)
% Draw a button in UI
uicontrol('Parent', hParent, ...
    'string', sStr, ...
    'unit', 'normalized', ...
    'position', vPos, ...
    'callback', sCallback);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetReferences(varargin)
% Set bregma/lambda marks

SetHelpStr('Click on bregma')
[nXbr, nYbr] = ginput(1);
DrawLabel(GetLineHandle('refs_bregma'), nXbr, nYbr, 'refs')
SetData('refs_bregma', [nXbr nYbr])

SetHelpStr('Click on lambda')
[nXlam, nYlam] = ginput(1);
DrawLabel(GetLineHandle('refs_lambda'), nXlam, nYlam, 'refs')
SetData('refs_lambda', [nXlam nYlam])
SetHelpStr('')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetMidline(varargin)
% Set mid-line axis

% Display image un-rotated
UpdateImage(0)

SetHelpStr('Click on two locations along the midline')
[vX, vY] = ginput(2);
SetData('refs_midline', [vX vY])
SetHelpStr('')

UpdateImage(1);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetElectrodes(varargin)
% Set electrode positions
SetHelpStr('Click on an electrode location, or outside image to finish.')
mImg = get(GetImage, 'cdata');

while 1
    [nX, nY] = ginput(1);
    % Abort if click is outside axis
    if nX < 0 || nX > size(mImg, 2) || nY < 0 || nY > size(mImg, 1)
        break
    end
    
    % Get electrode name
    csAns = inputdlg('Electrode name:', 'Name', 1);
    sTag = ['elec_' csAns{1}];
    DrawLabel(GetLineHandle(sTag), nX, nY, 'elec')
    SetData(sTag, [nX nY])
end
SetHelpStr('');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetContour(varargin)
% Set electrode positions
SetHelpStr('Draw a contour, and click outside image to finish.')
mImg = get(GetImage, 'cdata');

% Set contour name
csAns = inputdlg('Contour name (no spaces):', 'Name', 1);
if isempty(csAns), return; end
csAns{1} = strrep(csAns{1}, ' ', '_');

vX = []; vY = [];
while 1
    [nX, nY] = ginput(1);
    % Abort if click is outside axis
    if nX < 0 || nX > size(mImg, 2) || nY < 0 || nY > size(mImg, 1)
        break
    end
    vX(end + 1) = nX;
    vY(end + 1) = nY;

    sTag = ['cont_' csAns{1}];
    DrawLabel(GetLineHandle(sTag), vX, vY, 'cont')
    SetData(sTag, [vX' vY'])
end
SetHelpStr('');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawLabel(hDot, nXbr, nYbr, sType)
switch sType
    case 'refs'
        sCol = 'b';
        sMarker = 'x';
    case 'elec'
        sCol = 'g';
        sMarker = 'o';
    case 'cont'
        sCol = 'r';
        sMarker = '.';
    case 'cali'
        sCol = 'w';
        sMarker = '.';
    otherwise
        sCol = 'k';
        sMarker = '.';
end
set(hDot, 'xdata', nXbr, 'ydata',  nYbr, ...
    'marker', sMarker, 'LineWidth', 1.5, ...
    'color', sCol, 'markersize', 14);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowHideImage(varargin)
set(GetImage, 'AlphaData', ~get(GetImage, 'AlphaData'))
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetHelpStr(sStr)
% Set help string
hTxt = findobj(0, 'tag', 'EEGGetFiducialsMsg');
if ~isempty(hTxt)
    set(hTxt, 'string', sStr)
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hLine = GetLineHandle(sStr)
% Return existing line handle, or create it if not exist.
hLine = findobj(GetAxis, 'tag', sStr);
if isempty(hLine)
    hLine = plot(nan, nan, 'tag', sStr);

    % Attach menu to object
    hMenu = uicontextmenu;
    uimenu(hMenu, 'Label', sStr);
    uimenu(hMenu, 'Label', 'Rename', 'Callback', @RenameObject, 'userdata', hLine);
    uimenu(hMenu, 'Label', 'Delete', 'Callback', @RemoveObject, 'userdata', hLine);
    set(hLine, 'uicontextmenu', hMenu)
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RenameObject(hObj, varargin)
% Rename an object

% Remove entry in data structure
sTag = get(get(hObj, 'userdata'), 'tag');

tData = get(GetAxis, 'userdata');
if ~isfield(tData, sTag), return; end

% Get new label
csLabel = inputdlg('New label:', 'Rename', 1, {sTag});
if isempty(csLabel), return; end

% Abort if label was not changed
if strcmp(sTag, csLabel{1}), return; end

tData.(csLabel{1}) = tData.(sTag);
tData = rmfield(tData, sTag);
vXY = tData.(csLabel{1});

set(GetAxis, 'userdata', tData);

% Remove old and draw new graphical object
delete(get(hObj, 'userdata'))
iSep = strfind(csLabel, '_');
DrawLabel(GetLineHandle(csLabel{1}), vXY(1), vXY(2), csLabel{1}(1:(iSep{1} - 1)))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RemoveObject(hObj, varargin)
% Remove an object and its stored data

% Remove entry in data structure
sTag = get(get(hObj, 'userdata'), 'tag');
SetData(sTag, '')

% Remove graphical object
delete(get(hObj, 'userdata'))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetData(sField, nValue)
% Set data in user-data of image. Fields are overwritten if they already
% exist. Setting a field with no data (e.g. empty matrix) will remove the
% field.
tData = get(GetAxis, 'userdata');

% Append 'num_' to numeric fieldnames
if isstrprop(sField, 'digit')
    sField = sprintf('num_%s', sField);
end

if isempty(tData)
    tData = struct([]);
end
if isempty(nValue)
    if isfield(tData, sField)
        tData = rmfield(tData, sField);
        CmdMsg(sprintf('Removed "%s".', sField))
    end
else
    tData(1).(sField) = nValue;
    CmdMsg(sprintf('Updated "%s".', sField))
end
set(GetAxis, 'userdata', tData);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetCalibration(varargin)
SetHelpStr('Select pairs of features with same distance. Click outside image to finish.');
mImg = get(GetImage, 'cdata');

% Delete old calibration
nI = 1;
while 1
    sTag = sprintf('calib_%d', nI);
    hMarker = GetLineHandle(sTag);
    if isnan(get(hMarker, 'xdata')), break; end
    delete(hMarker)
    SetData(sTag, [])
    nI = nI + 1;
end

nI = 0;
vDist = [];
hCalib = [];
while 1
    [vX, vY] = ginput(2);
    % Abort if click is outside axis
    if any(vX < 0) || any(vX > size(mImg, 2)) || any(vY < 0) || any(vY > size(mImg, 1))
        break
    end
    nI = nI + 1;

    sTag = sprintf('cali_%d', nI);
    DrawLabel(GetLineHandle(sTag), vX, vY, 'cali')
    SetData(sTag, [vX vY])
    hCalib(end + 1) = GetLineHandle(sTag);
    vDist(end + 1) = sqrt(diff(vX)^2 + diff(vY)^2);
end

% Get calibration distance
if nI > 0
    while 1
        csMM = inputdlg('Calibration distance (mm):', 'Distance', 1);
        if ~isempty(csMM)
            if ~isempty(csMM{1}), break; end
        end
    end
    nMMPix = str2num(csMM{1}) / mean(vDist);
    SetData('calibration_mmpix', nMMPix)
end

% Calibration bar
DrawCalibrationBar();

% Hide calibration markers
set(hCalib, 'visible', 'off')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawCalibrationBar(varargin)
tData = get(GetAxis, 'userdata');
if ~isfield(tData, 'calibration_mmpix')
    return;
end
nMMPix = tData.calibration_mmpix;

% Draw 10 mm calibration bar
DrawLabel(GetLineHandle('cali_bar'), [10 10 + (10 / nMMPix)], [20 20], 'cali')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Export(varargin)
% Export data to .mat file
tData = get(GetAxis, 'userdata');
save(GetDataPath, 'tData')
CmdMsg(['Exported data to ' GetDataPath])
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hAx = LoadData(varargin)
if ~exist(GetDataPath, 'file')
    return;
end
load(GetDataPath, 'tData')
% Draw reference, electrode and contour objects on image
csFields = fieldnames(tData);
for f = 1:length(csFields)
    if length(csFields{f}) < 5
        continue;
    end
    if ~strcmp(csFields{f}(1:5), {'cont_' 'elec_' 'refs_' 'cali_'})
        continue
    end
    vX = tData.(csFields{f})(:, 1);
    vY = tData.(csFields{f})(:, 2);
    DrawLabel(GetLineHandle(csFields{f}), vX, vY, csFields{f}(1:4))
end
set(GetAxis, 'userdata', tData);

% Hide calibration and midline markers
set(FindObjects('cali_'), 'visible', 'off')
DrawCalibrationBar();
set(FindObjects('midline'), 'visible', 'off')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hHandles = FindObjects(sStr)
% Wildcard search for objects by tag
hHandles = findobj(GetAxis, '-regexp', 'tag', sStr);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sSavePath = GetDataPath(varargin)
tData = get(GetAxis, 'userdata');
[~, sFile, ~] = fileparts(tData.imagefile);
sSavePath = fullfile(tData.imagepath, sprintf('%s_EEGFiducials.mat', sFile));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CmdMsg(sMsg)
disp(sprintf('EEG_Get_Fiducials: %s', sMsg))
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hAx = GetAxis()
hAx = findobj(0, 'tag', 'EEGGetFiducials');
if isempty(hAx)
    hAx = get(findobj(0, 'tag', 'EEGGetFiducialsImg'), 'parent');
    if isempty(hAx)
        error('Cannot find axis!')
    end
end
hAx = hAx(end);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hImg = GetImage()
hImg = findobj(0, 'tag', 'EEGGetFiducialsImg');
if length(hImg) > 1
    hImg = hImg(end);
end
return
