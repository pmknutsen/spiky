function FV = import_spb(sFile, FV)
%Spiky files

% Internal Spiky sub-rutines can be called with the syntax:
%  Spiky.SUB(var)
%

% TODO
%   Permit this filter to be used to import partial data from one SPB file
%   into another SPB file.
%

global Spiky
FV = Spiky.main.SetFVDefaults();
FV.sDirectory = pwd;
FV.tData = struct([]);
FV.sLoadedTrial = [pwd sFile];
FV.tData = struct([]);

% The rest of the data is loaded with a call to OpenSettings() in LoadTrial()


return