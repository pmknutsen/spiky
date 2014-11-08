% Publish Spiky data into the base workspace of Matlab

% Clear workspace
clear variables

% Get database
FV = get(findobj('tag','Spiky'), 'UserData');

% Variables that will be retained in the workspace
sKeepvars = {};
sDisplayVars = {};

% Unit data
% Format: Double matrix with a row/column of timestamps in sec
if isfield(FV, 'tSpikes')
    csFieldnames = fieldnames(FV.tSpikes);
    hWait = waitbar(0, 'Publishing unit data');
    for nCh = 1:length(csFieldnames) % iterate over channels
        waitbar(nCh/length(csFieldnames), hWait)
        % Channel name
        sCh = csFieldnames{nCh};
        % Check that this channel has been sorted
        if ~isfield(FV.tSpikes.(sCh), 'hierarchy'), continue, end
        tHierarchy = FV.tSpikes.(sCh).hierarchy;
        % Get spiketimes
        vSpiketimes = FV.tSpikes.(sCh).spiketimes; % samples
        % Get sampling rate
        nStrIndx = strfind(sCh, '__');
        if ~isempty(nStrIndx) % stereo/tetrode
            nFs = FV.tData.([sCh(1:nStrIndx-1) '_KHz']) * 1000; % Hz
        else % monotrode
            nFs = FV.tData.([sCh '_KHz']) * 1000; % Hz
        end
        % Iterate over units
        vUnitIDs = unique(tHierarchy.assigns);
        for uid = vUnitIDs'
            vIndx = find(tHierarchy.assigns == uid);
            % Create variable in workspace that contains unit's spiketimes
            eval([sCh '_' num2str(uid) ' = vSpiketimes(vIndx) ./ nFs;']); % sec
            sKeepvars{end+1} = [sCh '_' num2str(uid)];
            sDisplayVars{end+1} = [sCh '_' num2str(uid) ' (' num2str(nFs/1000) ' KHz)'];
        end
    end
end
close(hWait)

% EMG data
% Format: One double vector for each channel
csEMGChannels = FV.csEMGChannels;
hWait = waitbar(0, 'Publishing EMG data');
for nCh = 1:length(csEMGChannels)
    waitbar(nCh/length(csEMGChannels), hWait)
    sCh = csEMGChannels{nCh};
    % get continuous signal
    vContThis = FV.tData.(sCh)'; % continuous trace (mV)    
    nFs = FV.tData.([sCh '_KHz']) * 1000; % sampling frequency (Hz)
    % decimate, filter and rectify UNLESS WE ARE IN MERGE_MODE (in which
    % case EMG has already been filtered)
    global g_bMergeMode
    if ~g_bMergeMode
        nNewFs = FV.nEMGLoPass*2.5;
        nStep = ceil(nFs / nNewFs);
        vContThis = vContThis(1:nStep:length(vContThis)); % decimate
        [b,a] = butter(3, FV.nEMGLoPass/(double(nNewFs)/2) , 'low'); % low-pass
        vContThis  = filtfilt(b, a, vContThis);
        [b,a] = butter(3, FV.nEMGHiPass/(double(nNewFs)/2) , 'high'); % high-pass
        vContThis  = filtfilt(b, a, vContThis);
        if FV.nEMGRectify, vContThis = abs(vContThis); end % rectify
    else, nNewFs = nFs; end

    % create continuous variable that starts at TimeBegin
    nBeginTime = round(FV.tData.([sCh '_TimeBegin']) * nNewFs) + 1; % samples
    nEndTime = FV.tData.([sCh '_TimeEnd']) * nNewFs + 1; % samples
    vCont = (1:nEndTime) .* 0;
    vCont(nBeginTime:(nBeginTime+length(vContThis)-1)) = vContThis;
    
    % create variables that will be accessible from NEX
    sVarName = ['EMG_' sCh];
    eval([sVarName ' = vCont(:);']);
    sKeepvars{end+1} = sVarName;
    sDisplayVars{end+1} = [sVarName ' (resampled at ' num2str(nNewFs/1000) ' KHz; ' sprintf('%.4f sec', 1/nNewFs) ')'];
end
close(hWait)

% Digital timestamps
% TODO


% Show report
% TODO: When showing EMG channels, note also sampling rate AFTER decimating!!!!
helpdlg([{'Published variables' ''} sDisplayVars {'' 'You may now import these variables if this instance of Matlab was initiated as an Engine from within NEX'}], 'Publish to NEX');

% Clear all variables except those that should be available from NEX
cAllVars = who;
for v = 1:length(cAllVars)
    if ~any(strcmp(sKeepvars, cAllVars{v})) & ~any(strcmp({'sKeepvars' 'cAllVars' 'v'}, cAllVars{v}))
        clear(cAllVars{v})
    end
end
clear cAllVars sKeepvars v


return
