function FV = import_openephys(sFile, FV)
%Open ePhys (openephys) files

% Open .openephys files in Spiky
%
% This filter imports all continuous data files in an Open ePhys session.
% The list of files to import are given in the .openephys meta file.
% Importing of events is currently not supported.
% 

% Internal Spiky sub-rutines can be called with the syntax:
%  Spiky.SUB(var)
%

% todo
% read .openephys file, xml?
% iterate over and load all referenced files (continuous, events, spikes)
%

global Spiky
[~, sFileOnly] = fileparts(sFile);

%% Initialize
FV = Spiky.main.SetFVDefaults();
FV.sDirectory = pwd;
FV.tData = struct([]);
FV.sLoadedTrial = fullfile(FV.sDirectory, sFile);

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

%% Read all .continuous files
hWait = waitbar(0, 'Loading continuous Open ePhys files...');
centerfig(hWait, Spiky.main.GetGUIHandle());
for c = 1:length(csCh)
    waitbar(c/length(csCh), hWait)
    if exist(csFiles{c}, 'file')
        [vData, ~, tInfo] = load_open_ephys_data_faster(csFiles{c});
        sCh = csCh{c};
        sTime = tInfo.header.date_created(end-5:end);
        vSecSinceMidnight = str2num(sTime(1:2)) * 60 * 60 ...
            + str2num(sTime(3:4)) * 60 ...
            + str2num(sTime(5:6));
        FV.tData(1).(sCh) = vData'; % microvolts
        FV.tData(1).([sCh '_KHz']) = tInfo.header.sampleRate / 1000;
        FV.tData(1).([sCh '_KHz_Orig']) = tInfo.header.sampleRate / 1000;
        FV.tData(1).([sCh '_TimeBegin']) = vSecSinceMidnight;
        FV.tData(1).([sCh '_TimeEnd']) = (length(vData) / tInfo.header.sampleRate) + vSecSinceMidnight;
    end
end
close(hWait)

return
