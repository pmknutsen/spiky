function FV = import_mat(sFile, FV)
%Matlab files

% Import vectors from .mat files in Spiky

% Internal Spiky sub-rutines can be called with the syntax:
%  Spiky.SUB(var)
%

global Spiky

% Import vectorized data from .mat files
tImportedData = load([sPath sFile], '-MAT');
cFields = fieldnames(tImportedData);

% Select data vector(s)
[cDataFields, nManualValue] = Spiky.main.SelectImportVariable(cFields, 'Set data vector(s)', {'Set data vector(s)'}, 0, NaN, p_cDataFields, 1);
p_cDataFields = cDataFields;
if isempty(cDataFields), return, end % user closed window

% Check data vector is a vector (minimum length 2)
for c = 1:length(cDataFields)
    if any(size(tImportedData.(cDataFields{c})) == 1) % vector data
        if length(tImportedData.(cDataFields{c})) < 2
            uiwait(warndlg(sprintf('Import failed. Input field %s has a length of less than 2 number. Must have a length greater than this.', cDataFields{c})));
            return
        end
    end
end

% Select sampling rate (kHz)
[cFsField, nManualValue] = Spiky.main.SelectImportVariable(cFields, 'Set sampling rate (Hz)', {'Set sampling rate (Hz)'}, 1, 40000, p_sFsField, 0);
if isempty(cFsField), return, end % user closed window
p_sFsField = cFsField{1};
if ~isnan(nManualValue), nFs = double(nManualValue) / 1000;
else nFs = double(tImportedData.(p_sFsField)) / 1000; end

% Select gain
[cGainField, nManualValue] = Spiky.main.SelectImportVariable(cFields, 'Set gain', {'Set gain'}, 1, 10000, p_sGainField, 0);
if isempty(cGainField), return, end % user closed window
p_sGainField = cGainField{1};
if ~isnan(nManualValue), nGain = double(nManualValue);
else nGain = double(tImportedData.(p_sGainField)); end

% Ask whether data should be appended to existing channels.
% Note: Append is not supported with matrix import (i.e. importing multiple channels at once)
bAppend = 1;
for c = 1:length(cDataFields)
    if numel(tImportedData.(cDataFields{c})) > sum(size(tImportedData.(cDataFields{c})))
        bAppend = 0;
    end
end

if ~isempty(FV.sLoadedTrial) && bAppend
    sAns = questdlg('Open in new workspace or append to existing?', 'Import data', 'New Workspace', 'Append', 'Cancel', 'Append');
    switch sAns
        case 'Append'
            % If data is to be appended, get onset time
            bAppend = 1;
            cFields = fieldnames(FV.tData);
            cDescriptions = {};
            vIndx = [];
            for i = 1:length(cFields)
                if ~isempty([strfind(cFields{i}, '_Up') strfind(cFields{i}, '_Down')])
                    if ~isempty(FV.tData.(cFields{i}))
                        vIndx(end+1) = i;
                        % Get channel description
                        nStartIndx = [strfind(cFields{i}, '_Up') strfind(cFields{i}, '_Down')];
                        sStrMatch = cFields{i}(1:nStartIndx-1);
                        nMatchIndx = find(strcmpi({FV.tChannelDescriptions.sChannel}, sStrMatch));
                        cDescriptions{i} = FV.tChannelDescriptions(nMatchIndx).sDescription;
                    end
                end
            end
            
            % Get alignment event
            [cONField, nManualValue, nButton] = SelectImportVariable([cFields(vIndx); 'Append to existing channels'], ...
                'Align to event', {'Align Start', 'Align End'}, 1, 0, p_sONField, 0, ...
                [cDescriptions(vIndx)'; ' ']);
            if isempty(cONField), return, end % user closed window
            p_sONField = cONField{1};
            switch nButton
                case 1 % align start of signal to trigger
                    sAlign = 'start';
                case 2 % align end of signal to trigger
                    sAlign = 'end';
            end
            
        case 'New Workspace', bAppend = 0;
        case 'Cancel', return;
    end
else bAppend = 0; end

% Dont append and reset FV structure
if ~bAppend
    FV = SetFVDefaults();
    FV.sDirectory = sPath;
    FV.tData = struct([]);
    FV.sLoadedTrial = [sPath sFile];
end

% Iterate over all fields to be imported
for c = 1:length(p_cDataFields)
    sDataField = p_cDataFields{c};
    
    % Get data and create inner loop for data matrices
    mData = tImportedData.(sDataField);
    
    % Transpose data vector, if needed, so that column(s) equal channels.
    vSize = size(mData);
    if vSize(2) > vSize(1)
        mData = mData';
    end
    
    for di = 1:size(mData, 2)
        
        % Temporal alignment offset
        if ~isnan(nManualValue)
            nOnsetTime = nManualValue;
        else
            if ~strcmpi(cONField, 'Append to existing channels')
                switch lower(sAlign)
                    case 'start' % align start of signal to trigger
                        nOnsetTime = FV.tData.(p_sONField)(1);
                    case 'end' % align end of signal to trigger
                        nOffsetTime = FV.tData.(p_sONField)(1); % end of signal, sec
                        %nDur = length(tImportedData.(sDataField)(:)) /  (nFs*1000);
                        nDur = length(mData(:,di)) /  (nFs*1000);
                        nOnsetTime = nOffsetTime - nDur; % start of signal, sec
                end
            end
        end
        
        if strcmpi(cONField, 'Append to existing channels')
            % Append new data to an existing field
            if isfield(FV.tData, sDataField)
                %FV.tData.(sDataField) = [FV.tData(1).(sDataField) double(tImportedData.(sDataField)(:)')];
                FV.tData.(sDataField) = [FV.tData(1).(sDataField) double(mData(:, di))];
                %FV.tData(1).([sDataField '_TimeEnd']) = FV.tData(1).([sDataField '_TimeEnd']) + ...
                %    (length(tImportedData.(sDataField)(:)') / (nFs*1000)); % sec
                FV.tData(1).([sDataField '_TimeEnd']) = FV.tData(1).([sDataField '_TimeEnd']) + ...
                    (length(mData(:, di)) / (nFs*1000)); % sec
            else
                uiwait(warndlg(sprintf('The channel %s does not exist and therefore cannot be appended with new data.', sDataField)));
                continue
            end
        else
            % Create field for data
            sChName = sprintf('%s%d', sDataField, di);
            FV.tGain(1).(sChName) = nGain;
            FV.tData(1).([sChName '_KHz']) = nFs;
            FV.tData(1).([sChName '_KHz_Orig']) = nFs;
            %FV.tData(1).([sChName '_TimeEnd']) = nOnsetTime + (length(tImportedData.(sDataField)(:)') / (nFs*1000)); % sec
            FV.tData(1).([sChName '_TimeEnd']) = nOnsetTime + (length(mData(:, di)) / (nFs*1000)); % sec
            %FV.tData(1).(sChName) = double(tImportedData.(sDataField)(:)');
            FV.tData(1).(sChName) = double(mData(:, di));
            FV.tData(1).([sChName '_Imported']) = 1; % Mark field as imported
            FV.csDisplayChannels{end+1} = sChName;
            FV.tData(1).([sChName '_TimeBegin']) = nOnsetTime; % sec
        end
    end % end of channels loop
    
end

return
