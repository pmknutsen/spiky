function [vCont, vTime, nFs] = spikyfilter_common_reference(vCont, vTime, nFs)
% Subtract a common reference signal (average across channels) from
% channel.
% 

global Spiky
persistent p_iSel p_vComRef p_nFs p_vTime

% Select which channels to user for common reference
[FV, ~] = Spiky.main.GetStruct();
if isempty(FV.csChannels), return; end
if isempty(p_iSel) || isempty(vCont)
    csCh = setdiff(FV.csChannels, FV.csDigitalChannels);
    if length(p_iSel) == length(csCh)
        p_iSel = ones(size(csCh));
    end
    p_iSel = Spiky.main.SelectChannelsUI(csCh, {}, p_iSel);
    p_vComRef = []; % reset common average
end
if isempty(vCont), return; end

% Compute common reference - store as persistent variable; only change of
% channels are changed
if isempty(p_vComRef)
    csCh = FV.csChannels(p_iSel);
    nN = 0;
    p_vComRef = [];
    p_nFs = [];
    for c = 1:length(csCh)
        if isfield(FV.tData, csCh{c})
            % Check that reference channels have the same Fs
            nFs = FV.tData.([csCh{c} '_KHz']) * 1000;
            if isempty(p_nFs)
                p_nFs = nFs;
            end
            if p_nFs ~= nFs
                Spiky.main.sp_disp('Common reference channels have different sample rate.');
                return;
            end
            
            if isempty(p_vComRef)
                p_vComRef = FV.tData.(csCh{c});
            else
                p_vComRef = p_vComRef + FV.tData.(csCh{c});
            end
            nN = nN + 1;
        end
    end
    p_vComRef = p_vComRef ./ nN;
    
    % Get same time segment of common reference as vCont
    csFieldnames = fieldnames(FV.tData);
    iTimeBegin = find(~cellfun(@isempty, strfind(csFieldnames, 'TimeBegin')));
    csTimeBegin = cell2mat(csFieldnames(iTimeBegin(1)));
    nTimeBegin = FV.tData.(csTimeBegin);
    p_vTime = (nTimeBegin + 1/nFs):(1/nFs):(nTimeBegin + length(p_vComRef) / p_nFs);
end

% Find corresponding start and end points in common reference signal
[~, iStart] = min(abs(p_vTime - vTime(1)));
[~, iEnd] = min(abs(p_vTime - vTime(end)));

% Subtract from channel
if length(iStart:iEnd) == length(vCont)
    vCont = vCont - p_vComRef(iStart:iEnd);
else
    Spiky.main.sp_disp('Could not subtract common reference.')
end

return