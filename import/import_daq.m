function FV = import_daq(sFile, FV)
%Data Acquisition Toolbox files

% Open .daq files in Spiky

% Internal Spiky sub-rutines can be called with the syntax:
%  Spiky.SUB(var)
%

global Spiky
sFile = Spiky.main.CheckFilename(sFile);
[~, sFileOnly] = fileparts(sFile);

try
    [mData, ~, ~, ~, tDAQInfo] = daqread(sFile);
catch
    sStr = sprintf('An error occurred when reading the file:\n%s\n\n%s\n\nThis file may be corrupted or truncated.', ...
        sFileOnly, lasterr);
    FV.sImportError = sStr;
    return
end
mData = single(mData); % convert to single precision to conserve memory

if isempty(mData)
    sStr = 'An error occurred during loading of .DAQ file: File appears to be empty.';
    FV.sImportError = sStr;
    return
end
tData = struct([]);

% Import all DAQ channels
cDaqChannels = {tDAQInfo.ObjInfo.Channel.HwChannel};
cDaqChannelNames = {tDAQInfo.ObjInfo.Channel.ChannelName};
for c = 1:length(cDaqChannels)
    % Channel name
    sChName = ['DAQ_' num2str(cDaqChannels{c})]; % use channel's hardware ID
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
    vTime = tDAQInfo.ObjInfo.InitialTriggerTime;
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
        nFs = tDAQInfo.ObjInfo.SampleRate; % Hz
        nTimeBegin = vSecSinceMidnight;
        nTimeEnd = (size(mData, 1) / tDAQInfo.ObjInfo.SampleRate) + vSecSinceMidnight;
        [vUpTimes vDownTimes] = Spiky.main.DigitizeChannel(mData(:,c), nThresh, nFs, nTimeBegin, nTimeEnd);
        
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
    tData.([sChName '_KHz']) = tDAQInfo.ObjInfo.SampleRate / 1000; % kHz
    tData.([sChName '_KHz_Orig']) = tDAQInfo.ObjInfo.SampleRate / 1000;
    tData.([sChName '_TimeBegin']) = vSecSinceMidnight;
    tData.([sChName '_TimeEnd']) = (size(mData, 1) / tDAQInfo.ObjInfo.SampleRate) + vSecSinceMidnight;
end

% Add DAQ events as digital signals
tEventLog = tDAQInfo.ObjInfo.EventLog;
for e = 1:length(tEventLog)
    sChName = ['DAQ_' tEventLog(e).Type];
    sChDescr = '';
    if isfield(tData, sChName)
        tData.([sChName '_Up']) = [tData.([sChName '_Up']) tEventLog(e).Data];
        tData.([sChName '_Down']) = [tData.([sChName '_Down']) tEventLog(e).Data];
    else
        tData.([sChName '_Up']) = (tEventLog(e).Data.RelSample / tDAQInfo.ObjInfo.SampleRate) + vSecSinceMidnight; % sec
        tData.([sChName '_Down']) = tData.([sChName '_Up']);
        tData.([sChName '_KHz']) = tDAQInfo.ObjInfo.SampleRate / 1000;
        tData.([sChName '_KHz_Orig']) = tDAQInfo.ObjInfo.SampleRate / 1000;
        tData.([sChName '_TimeBegin']) = vSecSinceMidnight;
        tData.([sChName '_TimeEnd']) = (size(mData, 1) / tDAQInfo.ObjInfo.SampleRate) + vSecSinceMidnight;
        FV.csDigitalChannels = unique([FV.csDigitalChannels sChName]);
    end
    FV.tChannelDescriptions(end+1).sChannel = sChName;
    FV.tChannelDescriptions(end).sDescription = sChDescr;
end

FV.tData = tData;

return
