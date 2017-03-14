function FV = import_openephys(sFile, FV)
%Open ePhys (openephys) files

% Open .openephys files in Spiky
%
% This filter imports continuous and event data from Open ePhys files.
% 
% 

% TODO
% read spike data
% add support for messages, i.e. save these in FV and display them in the
% data viewer
%

global Spiky

%% Get list of channels names and filenames
oOE = xmlread(sFile);
oChannels = oOE.getElementsByTagName('CHANNEL');
csCh = {};
csFiles = {};
for k = 0:oChannels.getLength-1
   oCh = oChannels.item(k);
   csCh(end+1) = oCh.getAttribute('name');
   csFiles(end+1) = oCh.getAttribute('filename');
end

%% Get experiment number (when multiple files are saved in same directory)
oExp = oOE.getElementsByTagName('EXPERIMENT');
nExpNumber = str2num(oExp.item(0).getAttribute('number'));

% Remove duplicates
[csCh, vI, ~] = unique(csCh);
csFiles = csFiles(vI);

% Rename channels from CHX to CHXX
for c = 1:length(csCh)
    if strcmp(csCh{c}(1:2), 'CH')
        csCh{c} = sprintf('CH%.2d', str2num(csCh{c}(3:end)));
    end
end

% Sort
[csCh, iOrder] = sort(csCh);
csFiles = csFiles(iOrder);
if isempty(csCh), return; end

% Select channels to import
iSelected = Spiky.main.SelectChannelsUI(csCh, csCh);
iRem = setdiff(1:length(csCh), find(iSelected));
csCh(iRem) = [];
csFiles(iRem) = [];

%% Initialize
FV = Spiky.main.SetFVDefaults();
FV.sDirectory = pwd;
FV.tData = struct([]);
FV.sLoadedTrial = fullfile(FV.sDirectory, sFile);

%% Read all .continuous files
hWait = waitbar(0, 'Loading continuous Open ePhys files...');
centerfig(hWait, Spiky.main.GetGUIHandle());
for c = 1:length(csCh)
    if ~ishandle(hWait)
        break;
    end
    waitbar(c/length(csCh), hWait)
    if exist(csFiles{c}, 'file')
        [vData, vTime, tInfo] = load_open_ephys_data_faster(csFiles{c});
        vData = single(vData);
        sCh = csCh{c};
        sTime = tInfo.header.date_created(end-5:end);
        FV.tData(1).(sCh) = vData'; % microvolts
        FV.tData(1).([sCh '_KHz']) = tInfo.header.sampleRate / 1000;
        FV.tData(1).([sCh '_KHz_Orig']) = tInfo.header.sampleRate / 1000;
        FV.tData(1).([sCh '_TimeBegin']) = vTime(1) / tInfo.header.sampleRate;
        FV.tData(1).([sCh '_TimeEnd']) = (length(vData) / tInfo.header.sampleRate) + (vTime(1) / tInfo.header.sampleRate);
        FV.tGain(1).(sCh) = 1;
        
        % Digitize 'ADC' channels (IO board inputs)
        if strcmp(sCh(1:3), 'ADC')
            [vUp, vDown] = Spiky.main.DigitizeChannel(vData', 1, tInfo.header.sampleRate, 0, inf);
            FV.tData(1).([sCh '_Up']) = vUp;
            FV.tData(1).([sCh '_Down']) = vDown;
        end
    end
end
close(hWait)

%% Get name of events file
if nExpNumber > 1
    sEventsFile = sprintf('all_channels_%d.events', nExpNumber);
else
    sEventsFile = 'all_channels.events';
end

%% Read events
tInfo.recNum
if exist(sEventsFile, 'file')
    [vCh, vTime, tInfo] = load_open_ephys_data_faster(sEventsFile);
    for ch = unique(vCh)'
        % Detect what type of event this is (TTL or message)
        vEventType = unique(tInfo.eventType(vCh == ch));
        switch vEventType(1)
            case 3
                sCh = sprintf('TTL%0.2d', ch);
            case 5
                sCh = sprintf('MSG%0.2d', ch);
        end
        
        % Find events
        iUp = (vCh == ch) & (tInfo.eventId == 1);
        iDown = (vCh == ch) & (tInfo.eventId == 0);
        
        % Convert timestamps (s) to samples (tInfo.sampleNum is all zeros)
        vTimeS = vTime * tInfo.header.sampleRate; % sample times
        
        % Save event data
        FV.tData(1).([sCh '_KHz']) = tInfo.header.sampleRate / 1000;
        FV.tData(1).([sCh '_KHz_Orig']) = tInfo.header.sampleRate / 1000;
        FV.tData(1).([sCh '_TimeBegin']) = vTime(1);
        FV.tData(1).([sCh '_TimeEnd']) = vTime(end);
        FV.tData(1).([sCh '_Up']) = vTime(iUp);
        FV.tData(1).([sCh '_Down']) = vTime(iDown);
        FV.csDigitalChannels = unique([FV.csDigitalChannels sCh]);
        FV.tGain(1).(sCh) = 1;
    end
end

return
