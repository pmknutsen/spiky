function FV = Convert_GalvoScan_Coordinates_To_Pixels(FV)
% Convert GalvoScanner coordinates from millimeters to pixels
%
% *When/where is the output of this used???*
%
% Pre-requisites:
%   1) Vessel image taken with GalvoScanner that indicates millimeters
%   2) Image processing toolbox
%
%

cd(FV.sDirectory)

% Load the GalvoScanner blue/vessel image WITH stereotactic coordinates
[sFile, sPath] = uigetfile('*.png', 'Select GalvoScanner vessel image w/coordinates');
if ~sFile, return; end
reference = imread(fullfile(sPath, sFile));

% Resize reference image to true size ([480 640])
reference = imresize(reference, [480 640]);

% Display image
hFig = figure;
hAx = axes;
imshow(reference)

% Click on 2 known reference points
title(sprintf('Click on two locations with known anterior/posterior and medial/lateral coordinates.\nIt is assumed that locations right and posterior to bregma are negative.'))

[vX, vY] = ginput(2);
hold(hAx, 'on');
hAnaDots = plot(vX, vY, 'c+');
hAnaTxt = text(vX, vY, {'  1', '  2'}, 'color', 'w');

% Get anatomical coordinates
cAns = inputdlg({'1 - Anterior/Posterior' '1 - Medial/Lateral' ...
    '2 - Anterior/Posterior' '2 - Medial/Lateral'}, 'Coordinates', 1);
delete([hAnaDots; hAnaTxt])
if isempty(cAns), return, end

% Check that supplied coordinates are ok
vAP = [str2num(cAns{1}) str2num(cAns{3})];
vML = [str2num(cAns{2}) str2num(cAns{4})];
if length(vAP) < 2 || length(vML) < 2
    errordlg('You must enter a value in each field.')
    return
end
if diff(vAP) == 0 || diff(vML) == 0
    errordlg('The two clicked points must be unique along each axis.')
    return
end

% Store coordinates
tAnaCoords = struct();
tAnaCoords(1).vAP = vAP;
tAnaCoords(1).vX  = [vX];
tAnaCoords(1).vML = vML;
tAnaCoords(1).vY  = [vY];

% Conversion factor
nXPixMM = diff(vX) / diff(vAP); % pix/mm
nYPixMM = diff(vY) / diff(vML); % pix/mm
nF = mean(abs([nXPixMM nYPixMM]));

% Convert millimeters to pixels
mXY_mm = [-FV.tData.GalvoScanML' -FV.tData.GalvoScanAP'] .* nF;
mXY_pix = TransformCoordsToReference(mXY_mm, tAnaCoords);
plot(mXY_pix(:,1), mXY_pix(:,2), 'k.')

%vBregma_px = TransformCoordsToReference([0 0], tAnaCoords);


% Plot coordinates on reference image
plot(mXY_pix(:,1), mXY_pix(:,2), 'k.')
axis image

title('Black dots are scan positions after re-scaling to pixels')

% Save new coordinates in FV
FV.tData.GalvoScanMLpx = mXY_mm(:, 1)';
FV.tData.GalvoScanMLpx_KHz = FV.tData.GalvoScanML_KHz;
FV.tData.GalvoScanMLpx_TimeBegin = FV.tData.GalvoScanML_TimeBegin;
FV.tData.GalvoScanMLpx_TimeEnd = FV.tData.GalvoScanML_TimeEnd;

FV.tData.GalvoScanAPpx = mXY_mm(:, 2)';
FV.tData.GalvoScanAPpx_KHz = FV.tData.GalvoScanAP_KHz;
FV.tData.GalvoScanAPpx_TimeBegin = FV.tData.GalvoScanAP_TimeBegin;
FV.tData.GalvoScanAPpx_TimeEnd = FV.tData.GalvoScanAP_TimeEnd;

return


%
% Transform [x y] anatomical coordinates to axis coordinates
% Function copied from galvoscanner.m
%
function mXY_out = TransformCoordsToReference(mXY_in, tAnaCoords)
% mXY_in should be in mm

% Clicked coordinates
vX = tAnaCoords.vX; vY = tAnaCoords.vY;

% Corresponding anatomical coordinates
vAP = tAnaCoords.vAP; vML = tAnaCoords.vML;

% Angle to rotate grid by
nA = atan2(diff(vY), diff(vX)); % rotation
nA_r = atan2(diff(vML), diff(vAP)); % rotation

% Find origin in default grid (as we need to rotate around this point)
% i) Rotate user coordinates so that AP/ML axes are at straight angles
nT = -(nA_r - nA);
vX_r2 = ( vX .* cos(nT)) + (vY .* sin(nT));
vY_r2 = (-vX .* sin(nT)) + (vY .* cos(nT));
% ii) Get calibration scale from user points
nXPixMM = diff(vX_r2) / diff(vAP); % pix/mm
nYPixMM = diff(vY_r2) / diff(vML); % pix/mm
nPixMM = (nXPixMM + nYPixMM) / 2;

% iii) Find origin (average of both points, for accuracy)
vO(1,1) = vX_r2(1) - (vAP(1) * nPixMM);
vO(1,2) = vY_r2(1) - (vML(1) * nPixMM);
vO(2,1) = vX_r2(2) - (vAP(2) * nPixMM);
vO(2,2) = vY_r2(2) - (vML(2) * nPixMM);
vO = mean(vO, 1); % origin in default grid

% iv) Rotate origin back to default grid (reverse of above)
vO_x = ( vO(1) .* cos(nT)) - (vO(2) .* sin(nT));
vO_y = ( vO(1) .* sin(nT)) + (vO(2) .* cos(nT));
vO = [vO_x vO_y]; % origin in user-generated grid

% 1) Translate user-generated grid so that origin is at [0 0]
mXY_out(:,1) = mXY_in(:,1) - vO(1);
mXY_out(:,2) = mXY_in(:,2) - vO(2);

% 2) Rotate grid matrix
nTheta = (nA_r - nA); % radians
vXPnts_rot = ( mXY_in(:,1) .* cos(nTheta)) + (mXY_in(:,2) .* sin(nTheta));
vYPnts_rot = (-mXY_in(:,1) .* sin(nTheta)) + (mXY_in(:,2) .* cos(nTheta));
mXY_out(:,1) = vXPnts_rot;
mXY_out(:,2) = vYPnts_rot;

% 3) Translate user-generated grid back
mXY_out(:,1) = mXY_out(:,1) + vO(1);
mXY_out(:,2) = mXY_out(:,2) + vO(2);

return

