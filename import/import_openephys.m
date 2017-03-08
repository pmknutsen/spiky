function FV = import_openephys(sFile, FV)
%Open ePhys (openephys) files

% Open .openephys files in Spiky
%
% This filter imports continuous and event data from Open ePhys files.
% 
% 

% TODO
% remember which channels were imported last time
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

% Remove duplicates
[csCh, vI, ~] = unique(csCh);
csFiles = csFiles(vI);

if isempty(csCh), return; end

%% Ask which channels to load
hFig = figure;
vPos = get(hFig, 'position');
nH = length(csCh) * 20;
set(hFig, 'position', [vPos(1) vPos(2) 200 nH+60], 'menu', 'none', ...
    'Name', 'Select channels', 'NumberTitle', 'off', ...
    'closeRequestFcn', 'global STAT; STAT=[]; delete(gcf)')
csChf = fliplr(csCh);
for nCh = 1:length(csCh)
    nVal = 1;
    sStr = csChf(nCh);
    if strcmp(csChf{nCh}(1:3), 'AUX')
        sStr = sprintf('%s (Accelerometer)', csChf{nCh});
        nVal = 0;
    elseif strcmp(csChf{nCh}(1:3), 'ADC')
        sStr = sprintf('%s (Digital)', csChf{nCh});
        nVal = 0;
    end
    uicontrol(hFig, 'Style', 'checkbox', 'Position', [15 (nCh+2)*20 nH 20], ...
        'String', sStr, 'HorizontalAlignment', 'left', 'value', nVal);
end
persistent p_STAT
if ~isempty(p_STAT)
    hChecks = findobj(gcf, 'style', 'checkbox');
    for c = 1:length(hChecks)
        set(hChecks(c), 'value', p_STAT{c})
    end
end
uicontrol(hFig, 'Style', 'checkbox', 'Position', [15 40 nH 20], ...
    'String', 'Select all', 'HorizontalAlignment', 'left', 'value', 1, ...
    'callback', 'set(findobj(gcf,''style'',''checkbox''),''value'',get(gcbo,''value''))');
uicontrol(hFig, 'Style', 'pushbutton', 'Position', [110 10 80 20], ...
    'Callback', 'global STAT;STAT=get(findobj(gcf,''style'',''checkbox''),''value'');delete(gcf)', ...
    'String', 'OK' );
set(findobj(hFig, 'style', 'checkbox'), 'backgroundcolor', get(hFig, 'color'))
centerfig(hFig)
uiwait(hFig)

% Filter out deselected channels
global STAT
if isempty(STAT)
    FV = false;
    return;
end
p_STAT = STAT;
iRem = setdiff(1:length(csCh), find(cell2mat(STAT(2:end))));
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
    waitbar(c/length(csCh), hWait)
    if exist(csFiles{c}, 'file')
        [vData, vTime, tInfo] = load_open_ephys_data_faster(csFiles{c});
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

%% Read events
if exist('all_channels.events', 'file')
    [vCh, vTime, tInfo] = load_open_ephys_data_faster('all_channels.events');
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
