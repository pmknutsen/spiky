function FV = import_set(sFile, FV)
%Axona dataset

% Open .set files in Spiky

% Internal Spiky sub-rutines can be called with the syntax:
%  Spiky.SUB(var)
%

% Imports the following Axona files:
%   EGF files (raw data)
%   INP files (inputs/outputs)
%
%

global Spiky
sPath = [pwd filesep];
mtint = readAllDACQdata(sPath, sFile);

% Get experiment variables
csFields = {'trial_date', 'trial_time', 'experimenter', 'comments', 'duration', 'sw_version'};
for sField = csFields
    sValue = key_value(sField{1}, mtint.header, 'string', 'exact');
    Spiky.main.NewExperimentVariable(sField{1}, sValue);
end
[FV, ~] = Spiky.main.GetStruct();

%% EGF (raw data)
vSecsSinceMidnight = [];
if isfield(mtint, 'egf')
    for i = 1:length(mtint.egf)
        vEGF = mtint.egf(i).EGF;
        nFs = mtint.egf(i).Fs;
        sHeader = mtint.egf(i).header;

        % Get hardware channel
        nEGFCh = key_value(sprintf('EEG_ch_%d', i), mtint.header, 'num', 'exact');
        sChName = sprintf('EGF_%d', nEGFCh);
        
        % Channel description (if alternative channel name exists)
        if ~isfield(FV, 'tChannelDescriptions')
            FV.tChannelDescriptions = struct([]);
        end
        FV.tChannelDescriptions(end+1).sChannel = sChName;
        FV.tChannelDescriptions(end).sDescription = '';
        
        % Begin time is counted as the number of seconds elapsed since last midnight
        sStartTime = key_value('trial_time', sHeader, 'string');
        csStartTime = strsplit(sStartTime, ':');
        vHourToSec = str2num(csStartTime{1})*60*60;
        vMinToSec = str2num(csStartTime{1})*60;
        vSecSinceMidnight = vHourToSec + vMinToSec + str2num(csStartTime{3});
        vSecsSinceMidnight(end+1) = vSecSinceMidnight;
        
        % Save channel data
        FV.tData(1).(sChName) = vEGF;
        FV.tData.([sChName '_KHz']) = nFs / 1000; % kHz
        FV.tData.([sChName '_KHz_Orig']) = nFs / 1000;
        FV.tData.([sChName '_TimeBegin']) = vSecSinceMidnight;
        FV.tData.([sChName '_TimeEnd']) = (length(vEGF) / nFs) + vSecSinceMidnight;
        
        % Get channel gain
        FV.tGain(1).(sChName) = key_value(sprintf('gain_ch_%d', nEGFCh), mtint.header, 'num', 'exact');
    end
else
    sStr = 'EGF (raw) data was not found in this dataset.';
    FV.sImportError = sStr;
    return
end

vSecSinceMidnight = mode(vSecsSinceMidnight);

%% INP (digital inputs/outputs)
if isfield(mtint, 'inp')
    for i = 1:length(mtint.inp)
        mInputs = mtint.inp(i).Inputs;
        mOutputs = mtint.inp(i).Outputs;
        
        % Parse inputs
        inpch = unique(mInputs(:, 1));
        for inp = inpch'
            sTag = sprintf('INP_%d', inp);
            vUpTimes = mInputs(mInputs(:, 1) == inp & mInputs(:, 3) == 1, 2) + vSecSinceMidnight;
            vDownTimes = mInputs(mInputs(:, 1) == inp & mInputs(:, 3) == 0, 2) + vSecSinceMidnight;
            if isempty(vUpTimes) || isempty(vDownTimes), continue
            else
                vDownTimes(vDownTimes <= vUpTimes(1)) = [];
                FV.tData.([sTag '_Up']) = vUpTimes;         % UP sample times
                FV.tData.([sTag '_Down']) = vDownTimes;     % DOWN sample times
                FV.tData.([sTag '_KHz']) = mtint.inp.Fs;
                FV.csDigitalChannels = unique([FV.csDigitalChannels sTag]);
                FV.tChannelDescriptions(end+1).sChannel = sTag;
                FV.tChannelDescriptions(end).sDescription = sTag;
            end
        end

        % Parse outputs
        outpch = unique(mOutputs(:, 1));
        for outp = outpch'
            sTag = sprintf('OUT_%d', outp);
            vUpTimes = mOutputs(mOutputs(:, 1) == outp & mOutputs(:, 3) == 1, 2) + vSecSinceMidnight;
            vDownTimes = mOutputs(mOutputs(:, 1) == outp & mOutputs(:, 3) == 0, 2) + vSecSinceMidnight;
            if isempty(vUpTimes) || isempty(vDownTimes), continue
            else
                vDownTimes(vDownTimes <= vUpTimes(1)) = [];
                FV.tData.([sTag '_Up']) = vUpTimes;         % UP sample times
                FV.tData.([sTag '_Down']) = vDownTimes;     % DOWN sample times
                FV.tData.([sTag '_KHz']) = mtint.inp.Fs;
                FV.csDigitalChannels = unique([FV.csDigitalChannels sTag]);
                FV.tChannelDescriptions(end+1).sChannel = sTag;
                FV.tChannelDescriptions(end).sDescription = sTag;
            end
        end
        
    end
end

%%


return


