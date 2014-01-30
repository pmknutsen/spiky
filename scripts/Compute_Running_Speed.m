function FV = Compute_Running_Speed(FV)
% Compute running speed from ticks of an incremental rotary encoder.
%
% It is assume that the ticks are in the IncrEncoder field in FV.tData and
% that each each equals 1 angular degree of rotation.
%
% The reported speed is in units of deg/sec.
%
%

global Spiky

% Get channel
csChannels = {FV.tChannelDescriptions.sChannel};
sEncoderCh = csChannels(strcmp({FV.tChannelDescriptions.sDescription}, 'IncrEncoder'));
if isempty(sEncoderCh)
    FV.ScriptError = 'Incremental encoder data could not be found.'; return
end

% Get data
if ~isfield(FV.tData, sEncoderCh{1})
    FV.ScriptError = 'Incremental encoder data could not be found.'; return
end

vData = FV.tData.(sEncoderCh{1});
if isempty(vData)
    FV.ScriptError = 'Incremental encoder data could not be found.'; return
end

% Find rotary ticks (1 per 1 degree of angle)
vData = round(vData./5);
vTicks = [0 diff(vData) == 1];

% Cumulative sum of ticks
vCumTicks = cumsum(vTicks);

% Decimate to 1k Hz sampling rate
nFs = FV.tData.([sEncoderCh{1} '_KHz']); % kHz

[vLpTicks, ~, ~] = Spiky.main.FilterChannel(vCumTicks, 1:length(vCumTicks), nFs*1000, 10, 0, 0, 'none');
nR = (nFs*1000)/1000;
vLpTicks = decimate(double(vLpTicks), nR);
nNewFs = (nFs*1000)/nR;
vLpTicks = vLpTicks * 1000;
vLpTicks = [0 (diff(vLpTicks))];
vLpTicks(vLpTicks < 0) = 0;

FV.tData.RunningSpeed = vLpTicks;
FV.tData.RunningSpeed_KHz = nNewFs / 1000; % kHz
FV.tData.RunningSpeed_TimeBegin = FV.tData.([sEncoderCh{1} '_TimeBegin']);
FV.tData.RunningSpeed_TimeEnd = FV.tData.([sEncoderCh{1} '_TimeEnd']);
FV.tData.RunningSpeed_Unit = 'deg/s';
FV.tChannelDescriptions(end+1) = struct('sChannel', 'RunningSpeed', 'sDescription', 'RunningSpeed');
FV.tGain.RunningSpeed = 1;

% Hide original trace and show new cleaned trace
FV.csDisplayChannels = unique([FV.csDisplayChannels 'RunningSpeed']);

return

