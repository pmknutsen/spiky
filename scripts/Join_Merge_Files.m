function FV = Join_Merge_Files(FV)
%
% ** THIS SCRIPT IS CURRENTLY BROKEN! **
%
% This script joins two existing Merge files. It assumes the Merge files
% have an identical internal structure, where only the data differs across
% the two.
%
% The joined Merge file is loaded and displayed after running the script.
%
% At present, only two Merge files can be joined at a time. To join more
% Merge files, you can run the script multiple times.
%

sPwd = cd;

% Select Merge files
[sFile1, sPath1] = uigetfile( {'*.spb';'*.*'}, 'Pick Merge file #1');
cd(sPath1)
[sFile2, sPath2] = uigetfile( {'*.spb';'*.*'}, 'Pick Merge file #2');

% Select save location and filename
[sOutputFile, sOutputPath] = uiputfile( {'*.spb';'*.*'}, 'Select location to save joined file');

cd(sPwd)

% Load Merge files
tM1 = load([sPath1 sFile1], '-MAT');
tM2 = load([sPath2 sFile2], '-MAT');

% Join Merge files (file #2 is joined to #1)
tDataNew = struct([]);
sFields = fieldnames(tM1.FV.tData);
vTimeBegin = [];
vTimeEnd = [];
for i = 1:length(sFields)
    % Check that field exists in both files
    if ~isfield(tM1.FV.tData, sFields{i}) || ~isfield(tM2.FV.tData, sFields{i})
        waitfor(warndlg(sprintf('The data field %s is missing from one Merge file. Sorry, cannot join these files.', sFields{i})))
        return
    end
    % Get TimeBegin and TimeEnd values (and use these values throughout)
    if isempty(vTimeBegin)
        vInd = findstr('_', sFields{i});
        if isempty(vInd), vInd = length(sFields{i}); end
        if isfield(tM1.FV.tData, [sFields{i}(1:vInd(end)) 'TimeBegin'])
            sBeginField = [sFields{i}(1:vInd(end)) 'TimeBegin'];
            sEndField = [sFields{i}(1:vInd(end)) 'TimeEnd'];
        elseif isfield(tM1.FV.tData, [sFields{i} '_TimeBegin'])
            sBeginField = [sFields{i} '_TimeBegin'];
            sEndField = [sFields{i} '_TimeEnd'];
        end
        vTimeBegin = [tM1.FV.tData.(sBeginField) tM2.FV.tData.(sBeginField)];
        vTimeEnd = [tM1.FV.tData.(sEndField) tM2.FV.tData.(sEndField)];
    end
    % Detect field type
    if ~isempty(strfind(sFields{i}, '_Up'))
        % Sort and concatenate Up events
        vUp1 = tM1.FV.tData.(sFields{i});
        vUp2 = tM2.FV.tData.(sFields{i}) - vTimeBegin(2) + vTimeEnd(1);
        tDataNew(1).(sFields{i}) = [vUp1 vUp2];
    elseif ~isempty(strfind(sFields{i}, '_Down'))
        % Sort and concatenate Down events
        vDown1 = tM1.FV.tData.(sFields{i});
        vDown2 = tM2.FV.tData.(sFields{i}) - vTimeBegin(2) + vTimeEnd(1);
        tDataNew(1).(sFields{i}) = [vDown1 vDown2];
    elseif ~isempty(strfind(sFields{i}, '_KHz'))
        % Check if sampling rates are same. If not, abort script.
        if tM1.FV.tData.(sFields{i}) ~= tM2.FV.tData.(sFields{i})
            waitfor(warndlg(sprintf('Sorry, cannot join these files.\n\nSampling rates %s are different between Merge files. If you expected the sampling-rate to be the same, check if channels were post-filtered with different bandpass filters in Spiky! during merging.', sFields{i})))
            return
        else
            tDataNew(1).(sFields{i}) = tM1.FV.tData.(sFields{i});
        end
    elseif ~isempty(strfind(sFields{i}, '_TimeBegin'))
        % TimeBegin is that of tM1
        tDataNew(1).(sFields{i}) = tM1.FV.tData.(sFields{i});
    elseif ~isempty(strfind(sFields{i}, '_TimeEnd'))
        % TimeEnd is that of vM1 + length of vM2
        tDataNew(1).(sFields{i}) = tM1.FV.tData.(sFields{i}) + (vTimeEnd(2) - vTimeBegin(2));
    elseif any(strcmp(sFields{i}, {'FileStart' 'FileEnd'}))
        % Concatenate sub-structures
        tFiles1 = tM1.FV.tData.(sFields{i});
        tFiles2 = tM2.FV.tData.(sFields{i});
        % Subtract time from tM2 timestamps
        for f = 1:length(tFiles2)
            tFiles2(f).Timestamp = tFiles2(f).Timestamp - vTimeBegin(2) + vTimeEnd(1);
        end
        tDataNew(1).(sFields{i}) = [tFiles1 tFiles2];
    else
        % Merge data vector
        vD1 = tM1.FV.tData.(sFields{i});
        vD2 = tM2.FV.tData.(sFields{i});
        tDataNew(1).(sFields{i}) = [vD1 vD2];
    end
end

% Finalize joined structure
FV = tM1.FV;
FV.tData = tDataNew;
FV.sLoadedTrial = [sOutputPath sOutputFile];
FV.sDirectory = sOutputPath;

% Save new, joined structure
save([sOutputPath sOutputFile], 'FV', '-MAT')
uiwait(msgbox('Marge files joined successfully!', 'Spiky!', 'modal'));

return