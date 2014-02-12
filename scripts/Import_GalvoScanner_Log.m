function FV = Import_GalvoScanner_Log(FV)
%
%

if isempty(FV.sLoadedTrial)
    sStr = 'No trial loaded. Aborting now.';
    uiwait(warndlg(sStr, 'Spiky!'))
    FV.ScriptError = sStr;
    return
end

% Filename
sTrial = FV.sLoadedTrial;
% Corresponding log-file name (assuming convention is kept...)
sGalvoLog = [sTrial(1:end-4) '_GalvoScannerLog.txt'];
% Change directory to that of currently loaded file
cd(FV.sDirectory)
% Check if log file exists
if ~exist(sGalvoLog, 'file')
    FV.ScriptError = sprintf('Found no corresponding GalvoScanner log for:\n%s', sTrial);
    return
end

% Load log file

% General purpose function for reading GalvoScanner log files
% Input is path to log file. Output is a matrix with data in the format:
%
hFile = fopen(sGalvoLog, 'r');
% Get header
cHeader = {};
while 1
    tLine = fgetl(hFile);
    if strcmpi(tLine, 'header_end') || ~ischar(tLine), break, end
    cHeader{end+1} = tLine;
end
mLog = [];
while 1
    tLine = fgetl(hFile);
    if ~ischar(tLine), break, end
    mLog(end+1,:) = str2num(strrep(tLine, ',', ' '));
end
fclose(hFile);

if isempty(mLog)
    FV.ScriptError = sprintf('Emtpy GalvoScanner log:\n%s', sTrial);
    return;
end

% Get LaserShutter Up/Down times
nLaserIndx = find(strcmpi({FV.tChannelDescriptions.sDescription}, 'LaserShutter'));
sField = FV.tChannelDescriptions(nLaserIndx).sChannel;
vLaserOn = FV.tData.([sField '_Up']); % sec
if isempty(vLaserOn)
    FV.ScriptError = sprintf('Found no data in GalvoScanner log for:\n%s', sTrial);
    return
end
nLaserTimeBegin = vLaserOn(1); % sec
vLaserOn = vLaserOn - nLaserTimeBegin;
vLaserOff = FV.tData.([sField '_Down']); % sec
nLaserTimeEnd = vLaserOff(end); % sec
vLaserOff = vLaserOff - nLaserTimeBegin;
vLaserFs = FV.tData.([sField '_KHz']) * 1000; % Hz

vLaserOn = round(vLaserOn * vLaserFs); % samples
vLaserOff = round(vLaserOff * vLaserFs); % samples

% Make sure vLaserOn and vLaserOff vectors are same length
if length(vLaserOn) > length(vLaserOff)
    vLaserOn = vLaserOn(1:length(vLaserOff));
end
if length(vLaserOn) < length(vLaserOff)
    vLaserOff = vLaserOff(1:length(vLaserOn));
end

% Check that all vLaserOff values are larger than vLaserOn
if any(vLaserOff < vLaserOn)
    FV.ScriptError = sprintf('On/off times do not match in GalvoScanner log for:\n%s', sTrial);
    return
end

% Get all unique positions in preserved order
vAP = mLog(:, 4);
vML = mLog(:, 5);
mPos = [vAP vML];
mPos = fliplr(mPos);
[dum,I] = unique(mPos, 'rows');
mPos = mPos(sort(I), :);
mPos = fliplr(mPos);
vAP = mPos(:,1);
vML = mPos(:,2);

% Check that number of LaserShutter pulses is the same as #pulses in log file
if length(vAP) ~= length(vLaserOn)
    % If the number of laser pulses is smaller than the number of
    % positions, then truncate the length of the latter
    if length(vAP) > length(vLaserOn)
        vAP = vAP(1:length(vLaserOn));
        vML = vML(1:length(vLaserOn));
    end
    % If the number of laser pulses is larger than the number of
    % positions, then truncate the length of the former
    if length(vAP) < length(vLaserOn)
        vLaserOn = vLaserOn(1:length(vAP));
        vLaserOff = vLaserOff(1:length(vAP));
    end
    FV.ScriptError = sprintf('Incorrect number of pulses. Length of dwell positions or laser-pulses has been automatically truncated.');
end

% Generate two continous signals from LaserShutter channel
vAPCont = zeros(vLaserOff(end), 1) .* NaN;
vMLCont = zeros(vLaserOff(end), 1) .* NaN;
for i = 1:length(vLaserOn)
    vRange = [vLaserOn(i):vLaserOff(i)] + 1;
    vAPCont( max([1 (vRange(1)-100)]):end ) = vAP(i);
    vMLCont( max([1 (vRange(1)-100)]):end ) = vML(i);
end

% Decimate vAPCont and vMLCont to 2.5 KHz (i.e. 1 ms resolution)
nNewFs = 2.5; % KHz
nND = round(vLaserFs / (nNewFs * 1000)); % Hz
if nND > 0
    nNewFs = (vLaserFs / nND) / 1000; % actual, new sampling rate (due to rounding)
    vAPCont = downsample(vAPCont, nND);
    vMLCont = downsample(vMLCont, nND);
else
    nNewFs = vLaserFs / 1000; % KHz
end

% Lets insert AP and ML as continuous data vectors
FV.tData.('GalvoScanAP') = vAPCont(:)';
FV.tData.('GalvoScanAP_Imported') = 1;
FV.tData.('GalvoScanAP_KHz') = nNewFs; % KHz
FV.tData.('GalvoScanAP_KHz_Orig') = nNewFs; % KHz
FV.tData.('GalvoScanAP_TimeBegin') = nLaserTimeBegin; % sec
FV.tData.('GalvoScanAP_TimeEnd') = nLaserTimeEnd; % sec
FV.tGain.('GalvoScanAP') = 1;

FV.tData.('GalvoScanML') = vMLCont(:)';
FV.tData.('GalvoScanML_Imported') = 1;
FV.tData.('GalvoScanML_KHz') = nNewFs; % KHz
FV.tData.('GalvoScanML_KHz_Orig') = nNewFs; % KHz
FV.tData.('GalvoScanML_TimeBegin') = nLaserTimeBegin; % sec
FV.tData.('GalvoScanML_TimeEnd') = nLaserTimeEnd; % sec
FV.tGain.('GalvoScanML') = 1;

% Make sure new channels are displayed when Spiky! is redrawn
FV.csDisplayChannels{end+1} = 'GalvoScanAP';
FV.csDisplayChannels{end+1} = 'GalvoScanML';
FV.csDisplayChannels = unique(FV.csDisplayChannels);

return
