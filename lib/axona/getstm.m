function [ts,type,value,keyval,Fs] = getstm(datafile)
% Read Axona .stm file
%
% Robin Hayman <r.hayman@ucl.ac.uk>
% Per Knutsen <p.m.knutsen@medisin.uio.no>

fid = fopen(datafile,'r');% start file i/o
if (fid == -1)
    error('getinp:fileIO','Could not open file %s',datafile)
end

header = {};
for i = 1:11
    textstring = fgetl(fid);
    sep = findstr(textstring, ' ');
    
    key = textstring(1:(sep-1));
    value = textstring((sep+1):end);
    header{end+1, 1} = strtrim(key);
    header{end, 2} = strtrim(value);
end

% Get timebase
Fs = key_value('timebase', header, 'num', 'exact');

% Get number of stored timestamps
nosamples = key_value('num_stm_samples', header, 'num', 'exact');
nosamples = nosamples + 1;  % bug in dacqUSB means number is one too few

fseek(fid,11,0);
vals = fread(fid, nosamples*4);
fclose(fid); % finish file i/o

vals = reshape(vals, [4, nosamples]);
%%
ts = zeros(nosamples, 1);
tic
if ispc
    for i = 1:numel(ts)
        ts(i,1) = swapbytes(typecast(uint8(vals(1:4, i)), 'uint32'));
    end
elseif isunix
    for i = 1:numel(ts)
        ts(i,1) = typecast(uint8(vals(1:4, i)), 'uint32');
    end
end
toc
%%


ts = ts ./ Fs; % convert ts into seconds

% Post process
% 1. Remove first and last pulse
ts = sort(ts);
ts([1 end]) = [];
% 2. Remove 2nd pulse from pairs with intervals 1/625
dts = diff(ts);
irem = round(dts*1000)/1000 == 0.016;
ts(irem) = [];

%ts = unique(ts);
%ts(end)

%%
figure
plot(ts, ones(size(ts)), 'o')
%%

return
