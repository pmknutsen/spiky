function FV = Clean_Whisker_Movements(FV)
% Clean whisker movement traces and configure filters
%
% This function performs several improvements on whisking traces:
%  - interpolate invalid values in WhiskerMotion
%  - configure optimal filters for WhiskerMotion and WhiskerAngle
%
% New channels created:
%   WhiskerMotion_clean
%   WhiskerAngle_clean      TODO
%

global Spiky
csChannels = {FV.tChannelDescriptions.sChannel};

% Get WhiskerMotion channel
sMotionCh = csChannels(strcmp({FV.tChannelDescriptions.sDescription}, 'WhiskerMotion'));
if isempty(sMotionCh), FV.ScriptError = 'WhiskerMotion channel not found.'; return; end

% Get WhiskerAngle channel
sAngleCh = csChannels(strcmp({FV.tChannelDescriptions.sDescription}, 'WhiskerAngle'));
if isempty(sAngleCh), FV.ScriptError = 'WhiskerAngle channel not found.'; return; end

% Get whisking data
if ~isfield(FV.tData, sMotionCh{1}), FV.ScriptError = 'WhiskerMotion data not found.'; return; end
if ~isfield(FV.tData, sAngleCh{1}), FV.ScriptError = 'WhiskerAngle data not found.'; return; end
%%
vMotion = FV.tData.(sMotionCh{1});
vAngle = FV.tData.(sAngleCh{1});
nFs = FV.tData.([sMotionCh{1} '_KHz']);
if isempty(vMotion), FV.ScriptError = 'WhiskerMotion data not found.'; return; end
if isempty(vAngle), FV.ScriptError = 'WhiskerAngle data not found.'; return; end

vRt = sqrt(vAngle.^2 + vMotion.^2);

% Correct motion vectors
% Find indices to replace with NaNs
vRem = vMotion < 0.1;
vAbs = abs(diff(vMotion));
vRem = vRem | ([0 (vAbs > prctile(abs(diff(vMotion)), 99.9))]);
vRem = vRem | (vMotion > prctile(vMotion, 99) );
vRem = vRem | (vMotion < prctile(vMotion, 1) );

vRt = sqrt(vAngle.^2 + vMotion.^2);
vRem = vRem | (vRt > prctile(vRt, 99.9) );
vRem = vRem | (vRt < prctile(vRt, .1) );

vMotion(vRem) = NaN;
vAngle(vRem) = NaN;

% Interpolate NaNs
vX = find(~isnan(vMotion));
vY = vMotion(~isnan(vMotion));
vXi = find(isnan(vMotion));
vYi = interp1(vX, vY, vXi, 'cubic');
vMotion(vXi) = vYi;

vX = find(~isnan(vAngle));
vY = vAngle(~isnan(vAngle));
vXi = find(isnan(vAngle));
vYi = interp1(vX, vY, vXi, 'cubic');
vAngle(vXi) = vYi;

% Low-pass filter at 2.5 KHz
[vMotion, ~, nNewFs] = Spiky.main.FilterChannel(vMotion, 1:length(vMotion), nFs*1000, 1000, 0, 0, 'decimate');
[vAngle, ~, nNewFs] = Spiky.main.FilterChannel(vAngle, 1:length(vAngle), nFs*1000, 1000, 0, 0, 'decimate');

% Invert and normalize to [0 1]
vMotion = vMotion .* -1;
vMotion = vMotion - min(vMotion);
vMotion = vMotion ./ max(vMotion);

% Insert cleaned version of WhiskerMotion into FV
FV.tData.([sMotionCh{1} '_clean']) = vMotion;
FV.tData.([sMotionCh{1} '_clean_KHz']) = nNewFs / 1000; % KHz
FV.tData.([sMotionCh{1} '_clean_TimeBegin']) = FV.tData.([sMotionCh{1} '_TimeBegin']);
FV.tData.([sMotionCh{1} '_clean_TimeEnd']) = FV.tData.([sMotionCh{1} '_TimeEnd']);
FV.tChannelDescriptions(end+1) = struct('sChannel', [sMotionCh{1} '_clean'], 'sDescription', 'WhiskerMotion_clean');

% Insert cleaned version of WhiskerAngle into FV
FV.tData.([sAngleCh{1} '_clean']) = vAngle;
FV.tData.([sAngleCh{1} '_clean_KHz']) = nNewFs / 1000; % KHz
FV.tData.([sAngleCh{1} '_clean_TimeBegin']) = FV.tData.([sAngleCh{1} '_TimeBegin']);
FV.tData.([sAngleCh{1} '_clean_TimeEnd']) = FV.tData.([sAngleCh{1} '_TimeEnd']);
FV.tChannelDescriptions(end+1) = struct('sChannel', [sAngleCh{1} '_clean'], 'sDescription', 'WhiskerAngle_clean');

% Update [0 100] Hz, no rectify
if ~isfield(FV, 'tFilteredChannels')
    FV.tFilteredChannels = struct([]);
else
    FV.tFilteredChannels(strcmp({FV.tFilteredChannels.sChannel}, [sMotionCh{1} '_clean'])) = [];
    FV.tFilteredChannels(strcmp({FV.tFilteredChannels.sChannel}, [sAngleCh{1} '_clean'])) = [];
end
if isempty(FV.tFilteredChannels)
    iNew = 1;
else
    iNew = length(FV.tFilteredChannels) + 1;
end
FV.tFilteredChannels(iNew).sChannel = [sMotionCh{1} '_clean'];
FV.tFilteredChannels(iNew).vBandpass = [0 100];
FV.tFilteredChannels(iNew).bRectify = 0;

FV.tFilteredChannels(iNew+1).sChannel = [sAngleCh{1} '_clean'];
FV.tFilteredChannels(iNew+1).vBandpass = [0 100];
FV.tFilteredChannels(iNew+1).bRectify = 0;

% Hide original trace and show new cleaned trace
FV.csDisplayChannels = setdiff(FV.csDisplayChannels, sMotionCh{1});
FV.csDisplayChannels = unique([FV.csDisplayChannels [sMotionCh{1} '_clean']]);
FV.tGain.([sMotionCh{1} '_clean']) = 1;

FV.csDisplayChannels = setdiff(FV.csDisplayChannels, sAngleCh{1});
FV.csDisplayChannels = unique([FV.csDisplayChannels [sAngleCh{1} '_clean']]);
FV.tGain.([sAngleCh{1} '_clean']) = 1;

return

