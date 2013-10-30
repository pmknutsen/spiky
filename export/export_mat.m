function FV = export_mat(sFile, FV)
%Matlab file

% Export data to .mat file. Each trace is exported as a vector.
%
% This function exports:
%  Currently viewed continuous channels
%  Analog channels are renamed to their descriptive names
%
% Todo:
%  Export spike data
%

global Spiky

[FV, ~] = Spiky.main.GetStruct();

% Create vectors to export
csFields = FV.csDisplayChannels;
tSpikyExport = struct([]);
tSpikyExport(1).FileStart = FV.tData.FileStart;
tSpikyExport.FileEnd = FV.tData.FileEnd;
tSpikyExport.OriginFile = FV.sLoadedTrial;
tSpikyExport.OriginDir  = FV.sDirectory;

for c = 1:length(csFields)
    csName = Spiky.main.GetChannelDescription(csFields{c});
    tSpikyExport.(csName) = FV.tData.(csFields{c});
    tSpikyExport.([csName 'KHz']) = FV.tData.([csFields{c} '_KHz']);
    tSpikyExport.([csName '_TimeBegin']) = FV.tData.([csFields{c} '_TimeBegin']);
end

save(sFile, 'tSpikyExport')

return
