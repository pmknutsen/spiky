function FV = Convert_Barrel_Coordinates_To_Metric(FV)
% Convert barrel coordinates to millimeters
%
% Pre-requisites:
%   - Image Processing Toolbox
%   - run ISI_analysisGUI on all IOS files with the Export to MAT option
%   - run the Map_Barrels plugin from ISI_analysisGUI
%   - run the Import_Barrels script in Spiky
% 
% To run this function you also need:
%   - A vessel image taken with GalvoScanner that indicates at least 2 pairs of [x y]
%     coordinates in millimeters
%
% TODO:
%   Barrels convert to GalvoScanner image ok, but rotation of scan grid to
%   anatomical coordinates is buggy: It seems to only work if [x y] pairs
%   with identical numbers is used, e.g. [-2 -2] and [-3 -3].
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
hold on

% Click on 2 known reference points
title(sprintf('Click on two locations with known anterior/posterior and medial/lateral coordinates.\nIt is assumed that locations right and posterior to bregma are negative.'))

[vX, vY] = ginput(2);
hAnaDots = plot(vX, vY, 'c+');
hAnaTxt = text(vX, vY, {'  1', '  2'}, 'color', 'w');

% Get anatomical coordinates interactively
cAns = inputdlg({'1 - Anterior/Posterior' '1 - Medial/Lateral' ...
    '2 - Anterior/Posterior' '2 - Medial/Lateral'}, 'Coordinates', 1);
delete([hAnaDots; hAnaTxt])
if isempty(cAns), return, end

% Check that supplied coordinates are ok
vAP = [str2double(cAns{1}) str2double(cAns{3})];
vML = [str2double(cAns{2}) str2double(cAns{4})];
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
tAnaCoords(1).vX  = vX;
tAnaCoords(1).vML = vML;
tAnaCoords(1).vY  = vY;

% Conversion factor
nXPixMM = diff(vX) / diff(vAP); % pix/mm
nYPixMM = diff(vY) / diff(vML); % pix/mm
nF = mean(abs([nXPixMM nYPixMM]));

% Iterate over barrels
if isfield(FV.tData, 'BarrelOutlines')
    for b = 1:length(FV.tData.BarrelOutlines)
        mXYpx = [FV.tData.BarrelOutlines(b).vX' FV.tData.BarrelOutlines(b).vY']; % px
        figure(hFig); hold on
        plot(mXYpx(:, 1), mXYpx(:, 2), 'r-')
    end
end
axis image

% Plot barrels over all scan positions
hFig = figure;
plot(FV.tData.GalvoScanAP, FV.tData.GalvoScanML, 'o', 'color', [.8 .8 1])
axis image
hold on

% Iterate over barrels
if isfield(FV.tData, 'BarrelOutlines')
    for b = 1:length(FV.tData.BarrelOutlines)
        % Countour
        mXYpx = [FV.tData.BarrelOutlines(b).vX' FV.tData.BarrelOutlines(b).vY']; % px
        mXYmm = TransformCoordsToAnatomy(mXYpx, tAnaCoords); % mm
        mXYmm = mXYmm ./ nF;
        FV.tData.BarrelOutlines(b).vXmm = mXYmm(:, 1);
        FV.tData.BarrelOutlines(b).vYmm = mXYmm(:, 2);
        % Centroid
        nCentroidXYmm = TransformCoordsToAnatomy(FV.tData.BarrelOutlines(b).nCentroidXY, tAnaCoords); % mm
        nCentroidXYmm = nCentroidXYmm ./ nF;
        FV.tData.BarrelOutlines(b).nCentroidXYmm = nCentroidXYmm;
        % Draw barrel and centroid
        figure(hFig)
        plot(mXYmm(:, 2), mXYmm(:, 1), 'k-')
        plot(nCentroidXYmm(:, 2), nCentroidXYmm(:, 1), 'k.')
        hTxt = text(nCentroidXYmm(:, 2), nCentroidXYmm(:, 1), FV.tData.BarrelOutlines(b).sID);
        set(hTxt, 'horizontalAlignment', 'center', 'fontsize', 8)
    end
end
axis image
xlabel('AP (mm)')
ylabel('ML (mm)')

return


% --------------------------------------------------------------------
% Transform [x y] axis coordinates to anatomical coordinates
function mXY_out = TransformCoordsToAnatomy(mXY_in, tAnaCoords)

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
%vO(1,1) = vX_r2(1) - (vAP(1) * nPixMM);
%vO(1,2) = vY_r2(1) - (vML(1) * nPixMM);
%vO(2,1) = vX_r2(2) - (vAP(2) * nPixMM);
%vO(2,2) = vY_r2(2) - (vML(2) * nPixMM);
vO(1,1) = vX_r2(1) - (vAP(1) * nPixMM);
vO(1,2) = vY_r2(1) - (vML(1) * nPixMM);
vO(2,1) = vX_r2(2) - (vAP(2) * nPixMM);
vO(2,2) = vY_r2(2) - (vML(2) * nPixMM);
vO = mean(vO, 1); % origin in default grid

% iv) Rotate origin back to default grid (reverse of above)
vO_x = ( vO(1) .* cos(nT)) - (vO(2) .* sin(nT));
vO_y = ( vO(1) .* sin(nT)) + (vO(2) .* cos(nT));
vO = [vO_x vO_y]; % origin in user-generated grid

% 1) Translate input points in axis coordinates so origin is at [0 0]
mXY_out(:,1) = mXY_in(:,1) - vO(1);
mXY_out(:,2) = mXY_in(:,2) - vO(2);

% 2) Rotate grid matrix
nTheta = -(nA_r - nA); % radians
vXPnts_rot = ( mXY_out(:,1) .* cos(nTheta)) + (mXY_out(:,2) .* sin(nTheta));
vYPnts_rot = (-mXY_out(:,1) .* sin(nTheta)) + (mXY_out(:,2) .* cos(nTheta));
mXY_out(:,1) = vXPnts_rot;
mXY_out(:,2) = vYPnts_rot;
return

