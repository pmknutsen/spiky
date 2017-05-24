function FV = spiky_Import_OpenEphys_Dataset(FV)
% Import Open Ephys dataset using the Choi/Korea 40-point thin-film electrode
% 
% Status:
%   Extension works. No major items on todo list.
%
%

global Spiky;

% Get current directory
if isfield(FV, 'sDirectory')
    sDir = [FV.sDirectory filesep];
else
    sDir = pwd;
end

% Choose a dataset
[sFile, sPath] = uigetfile('.openephys', 'Select data file', sDir);
if isempty(sFile), return; end
sFile = fullfile(sPath, sFile);
if ~exist(sFile), return; end

% Load data
Spiky.main.sp_disp('Load data...')
Spiky.main.LoadTrial(sFile);

% Load channel settings
Spiky.main.sp_disp('Import channel settings...')
Spiky.main.ImportChannelSettings('E:\Data\OpenEphys_ChannelSettings.mat');

% Delete hidden channels
Spiky.main.sp_disp('Delete hidden channels...')
Spiky.main.DeleteHiddenChannels();

% Display only first 4 channels
[FV, hWin] = Spiky.main.GetStruct();
FV.csDisplayChannels = FV.csDisplayChannels(1:2);
Spiky.main.SetStruct(FV);

% Display events
set(findobj(hWin, 'Label', 'Show Events'), 'Checked', 'on')

% Update UI
Spiky.main.sp_disp('Update UI...')
Spiky.main.ViewTrialData();

% Save FV to disk
Spiky.main.sp_disp('Save SPB file to disk...')
Spiky.main.SaveResults();

Spiky.main.sp_disp('Done!')

return
