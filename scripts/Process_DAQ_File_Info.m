function FV = Process_DAQ_File_Info(FV)
% This script automatically sets gains, select set channels for
% display, imports GalvoScan logs, digitizes channels. It makes
% specific assumptions about my dataset.
%
% TODO:
%   Enable and configure filter on AccMotion channel
%

hWait = waitbar(0, 'Pre-Processing DAQ file...');

% Default gain on all channels to 1
csFieldnames = fieldnames(FV.tGain);
for fn = 1:length(csFieldnames), FV.tGain.(csFieldnames{fn}) = 1; end
waitbar(1/10, hWait)

% Set gain on electrodes to 8080 if they are currently 1
if isfield(FV.tGain, 'DAQ_0')
    if FV.tGain.DAQ_0 == 1, FV.tGain.DAQ_0 = 8080; end
end
if isfield(FV.tGain, 'DAQ_1')
    if FV.tGain.DAQ_1 == 1, FV.tGain.DAQ_1 = 8080; end
end
if isfield(FV.tGain, 'DAQ_2')
    if FV.tGain.DAQ_2 == 1, FV.tGain.DAQ_2 = 8080; end
end
if isfield(FV.tGain, 'DAQ_3')
    if FV.tGain.DAQ_3 == 1, FV.tGain.DAQ_3 = 8080; end
end
waitbar(2/10, hWait)

% Import GalvoScan log file
if ~isfield(FV.tData, 'GalvoScanML')
    FV = Import_GalvoScanner_Log(FV);
end
waitbar(3/10, hWait)

% Threshold known digital channels
csDigChannels = {'LaserShutter' 'FlashOut' 'ExSyncTrig' 'AirPuff'};
for i = 1:length(csDigChannels)
    csDescr = csDigChannels{i};
    nIndx = find(strcmp({FV.tChannelDescriptions.sDescription}, csDescr));
    sCh = FV.tChannelDescriptions(nIndx).sChannel;
    
    % Select channel
    vContData = FV.tData.(sCh);
    nFs = FV.tData.([sCh '_KHz']) * 1000; % sampling frequency (Hz)
    nBeginTime = FV.tData.([sCh '_TimeBegin']); % start of sampling (sec)
    nEndTime = FV.tData.([sCh '_TimeEnd']); % start of sampling (sec)

    nThresh = 2.5;
    
    % Threshold signal and find UP and DOWN times
    vIndxHIGH = find(vContData >= nThresh);
    vIndxLOW = find(vContData < nThresh);

    % Check UP and DOWN times if any were found
    if ~isempty(vIndxHIGH) && ~isempty(vIndxLOW)
        vI = (find(diff(vIndxHIGH) > 1) + 1);
        vIndxUP = vIndxHIGH(unique([1; vI(:)]));
        vI = (find(diff(vIndxLOW) > 1) + 1);
        vIndxDOWN = vIndxLOW(unique([1; vI(:)]));
        
        % Drop index==1 (because no event could have occurred in first frame)
        vIndxUP(vIndxUP == 1) = [];
        vIndxDOWN(vIndxDOWN == 1) = [];
        
        % Compute corresponding time vector
        vTime = linspace(nBeginTime, nEndTime, length(vContData));
        
        % Convert indices to temporal values
        vUpTimes = vTime(vIndxUP);      % sec
        vDownTimes = vTime(vIndxDOWN);  % sec
        
        % Subtract 1/2 timesample from indices; the mid-point is the best
        % estimate for when the signal actually changed state
        vUpTimes = vUpTimes - ((1/nFs)/2);
        vDownTimes = vDownTimes - ((1/nFs)/2);
    else
        vUpTimes = [];
        vDownTimes = [];
    end

    % Save UP/DOWN times in data structure
    FV.tData.([sCh '_Up']) = vUpTimes;         % UP sample times
    FV.tData.([sCh '_Down']) = vDownTimes;     % DOWN sample times
    FV.csDisplayChannels = {FV.csDisplayChannels{~strcmp(FV.csDisplayChannels, sCh)}};
    FV.csDigitalChannels = unique([FV.csDigitalChannels sCh]);

    % Save threshold value used for digitization
    if ~isfield(FV, 'tEventThresholds')
        FV.tEventThresholds = struct([]);
        nIndx = 1;
    else nIndx = length(FV.tEventThresholds) + 1; end
    FV.tEventThresholds(nIndx).sChannel = sCh;
    FV.tEventThresholds(nIndx).nThreshold = nThresh;        
end

% Display only electrodes
FV.csDisplayChannels = {'DAQ_0', 'DAQ_1', 'DAQ_2', 'DAQ_3'};

% Done and close waitbar
waitbar(1/1, hWait)
close(hWait)

% Command(s) to execute when script terminates
FV.ScriptExitCommand = 'global g_bBatchMode;g_bBatchMode=1;ShowEvents;ViewTrialData;PCACleaning;AutoDetectThresholds;QuickSort;g_bBatchMode=0;';

return