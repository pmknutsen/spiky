function FV = import_hd5(sFile, FV)
%Hierarchical Data Format (HDF) files

% Open .hd5 files in Spiky

% Internal Spiky sub-rutines can be called with the syntax:
%  Spiky.SUB(var)
%

global Spiky
[~, sFileOnly] = fileparts(sFile);

% Check that Matlab version 13+ is running (HD5 read otherwise not supported)
tV = ver;
iMat = strcmpi({tV.Name}, 'MATLAB');
if str2num(tV(iMat).Version) < 8.0
    sStr = sprintf('HDF5 import requires at least Matlab version 8.0 (R2012b). You are running version %s %s.', ...
        tV(iMat).Version, tV(iMat).Release);
    FV.sImportError = sStr;
    return
end

try
    % Load NI
    load(sFile, '-mat');
    
    % Load /NI/Data
    mData = h5read(sFile, '/NI/Data')';
catch
    sStr = sprintf('An error occurred when reading the file:\n%s\n\n%s\n\nThis file may be corrupted or truncated.', ...
        sFileOnly, lasterr);
    FV.sImportError = sStr;
    return
end
mData = single(mData); % convert to single precision to conserve memory

if isempty(mData)
    sStr = 'An error occurred during loading of .HD5 file: File appears to be empty.';
    FV.sImportError = sStr;
    return
end
tData = struct([]);

% Import all DAQ channels
cDaqChannels = {NI.Channels.ID};
cDaqChannelNames = {NI.Channels.Name};
for c = 1:length(cDaqChannels)
    % Channel name
    sChName = num2str(cDaqChannels{c}); % use channel's hardware ID
    
    % Channel description (if alternative channel name exists)
    if isempty(str2num(cDaqChannelNames{c})) % use channel name
        sChDescr = num2str(cDaqChannelNames{c}); % use channel's hardware ID
        if ~isfield(FV, 'tChannelDescriptions')
            FV.tChannelDescriptions = struct([]);
        end
    else sChDescr = ''; end
    FV.tChannelDescriptions(end+1).sChannel = sChName;
    FV.tChannelDescriptions(end).sDescription = sChDescr;
    
    % Begin time is counted as the number of seconds elapsed since last midnight
    vTime = NI.AbsStartTime;
    vHourToSec = vTime(4)*60*60;
    vMinToSec = vTime(5)*60;
    vSecSinceMidnight = vHourToSec + vMinToSec + vTime(6);
    
    % Auto-digitize bi-modal signals (putative digital inputs)
    %  1) round
    nMax = max(mData(:,c));
    %  2) normalize range to 0 -> 10
    vNorm = mData(:,c) - min(mData(:,c));
    vNorm = (vNorm / nMax) * 10;
    %  3) round
    vRnd = round( vNorm );
    %  4) get number of unique values. As a heuristic, we assume a
    %     digital signal has less than 5
    nLenU = length(unique(vRnd));
    %  5) get threshold; (max-min) / 2. If signal is a TTL,
    %  assume threshold to be > 1 V
    nThresh = (nMax - min(mData(:,c)) ) / 2;
    if (nLenU < 5) && nThresh > 1
        %  6) digitize through @DigitizeChannel
        nFs = NI.Rate; % Hz
        nTimeBegin = vSecSinceMidnight;
        nTimeEnd = (size(mData, 1) / nFs) + vSecSinceMidnight;
        [vUpTimes, vDownTimes] = Spiky.main.DigitizeChannel(mData(:,c), nThresh, nFs, nTimeBegin, nTimeEnd);
        
        tData(1).([sChName '_Up']) = vUpTimes; % sec
        tData.([sChName '_Down']) = vDownTimes; % sec
        
        % Save threshold value used for digitization
        if ~isfield(FV, 'tEventThresholds')
            FV.tEventThresholds = struct([]);
            nIndx = 1;
        else nIndx = length(FV.tEventThresholds) + 1; end
        FV.tEventThresholds(nIndx).sChannel = sChName;
        FV.tEventThresholds(nIndx).nThreshold = nThresh;
        
        tData.([sChName '_EventThreshold']) = nThresh; % keep this for backwards compatibility...
        
        FV.csDigitalChannels = unique([FV.csDigitalChannels sChName]);
    end
    
    tData(1).(sChName) = mData(:, c)'; % save raw data
    tData.([sChName '_KHz']) = NI.Rate / 1000; % kHz
    tData.([sChName '_KHz_Orig']) = NI.Rate / 1000;
    tData.([sChName '_TimeBegin']) = vSecSinceMidnight;
    tData.([sChName '_TimeEnd']) = (size(mData, 1) / NI.Rate) + vSecSinceMidnight;
end

FV.tData = tData;

return
