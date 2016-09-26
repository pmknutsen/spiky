function mtint = readDACQdata ( filepath, flnmroot )
% reads in all files associated with a single DACQ recording session and
% assembles into a structure called mtint.
cd(filepath)
% ---------------- add some metadata to the structure ----------------
mtint.flnmroot = flnmroot;
mtint.filepath = filepath;
mtint.header = getDACQHeader ( [filepath,flnmroot,'.set'], 'set' );

% ---------------- start read pos ---------------------
% find out how many leds were tracked and read pos file appropriately
if exist([filepath,flnmroot,'.pos'], 'file')
    idx = find(strcmpi('colactive_1',mtint.header(:,1)));
    n_leds = numel(find(str2double(char(mtint.header(idx:idx+3,2)))));
    [led_pos,post,led_pix] = rawpos([filepath,flnmroot,'.pos'],n_leds); % raw led positions
    % allocate to structure
    mtint.pos.led_pos = led_pos;
    mtint.pos.ts = post;
    mtint.pos.led_pix = led_pix;
    % get header info
    header = getDACQHeader ( [filepath,flnmroot,'.pos'], 'pos' );
    mtint.pos.header = header;
end
% ---------------- end read pos ---------------------

% list all the files that have the flnmroot
filelist = dir([filepath,flnmroot,'*']);
% extract only the tetrode files from this list
for ifile = 1:numel(filelist)
    type(ifile) = str2double(filelist(ifile).name(strfind(filelist(ifile).name,'.')+1:end));
end
type(isnan(type)) = 0;
filelist = filelist(logical(type));

% ---------------- start read tetrodes ---------------------
% read tetrode files
for ifile = 1:numel(filelist)
    if exist([filelist(ifile).name(1:end-2) '.tet'], 'file')
        current_tet = str2double(filelist(ifile).name(strfind(filelist(ifile).name,'.')+1:end));
        if isfinite(current_tet)
            % get header info
            header = getDACQHeader ( [filepath,flnmroot,'.',num2str(current_tet)], 'tet' );
            mtint.tetrode(ifile).id = current_tet;
            mtint.tetrode(ifile).header = header;
            ts = getspikes([filepath,flnmroot,'.',num2str(current_tet)]);
            mtint.tetrode(ifile).ts = ts;
            %         [ts,ch1,ch2,ch3,ch4] =
            %         getspikes([filepath,flnmroot,num2str(current_tet)]); % uncomment
            %         this line and the 4 below to get the spikes on each channel
            %         mtint.tetrode(current_tet).ch1 = ch1;
            %         mtint.tetrode(current_tet).ch2 = ch2;
            %         mtint.tetrode(current_tet).ch3 = ch3;
            %         mtint.tetrode(current_tet).ch4 = ch4;
            mtint.tetrode(ifile).pos_sample = ceil(ts * 50);
            mtint.tetrode(ifile).cut = [];
        end
    end
end
% ---------------- end read tetrodes ---------------------

% ---------------- start read cuts ---------------------
% read them in and assign to structure
if isfield(mtint, 'tetrode')
    for ifile = 1:numel(mtint.tetrode)
        current_tet = mtint.tetrode(ifile).id;
        if exist([filepath,flnmroot,'_',num2str(current_tet),'.cut'],'file')
            clust = getcut([filepath,flnmroot,'_',num2str(current_tet),'.cut']);
        else
            clust = [];
        end
        mtint.tetrode(ifile).cut = clust;
    end
end
% ---------------- end read cuts ---------------------


% ---------------- start read eeg ---------------------
% list all eeg files
eegfilelist = dir([filepath,flnmroot, '.eeg*']);
% read them in and assign to structure
for ifile = 1:numel(eegfilelist)
    [EEG,Fs] = geteeg([filepath,eegfilelist(ifile).name]);
    mtint.eeg(ifile).EEG = EEG;
    mtint.eeg(ifile).Fs = Fs;
    % get eeg header info
    header = getDACQHeader ( [filepath,eegfilelist(ifile).name], 'eeg' );
    mtint.eeg(ifile).header = header;
end
% ---------------- end read eeg ---------------------


% ---------------- start read egf ---------------------
% list all egf files
egffilelist = dir([filepath,flnmroot, '.egf*']);
% read them in and assign to structure
for ifile = 1:numel(egffilelist)
    [EGF,Fs] = getegf([filepath, egffilelist(ifile).name]);
    mtint.egf(ifile).EGF = EGF;
    mtint.egf(ifile).Fs = Fs;
    % get eeg header info
    header = getDACQHeader([filepath,egffilelist(ifile).name], 'eeg' );

    % Get file suffix (after .efg)
    iS = strfind(egffilelist(ifile).name, '.egf');
    sCh = egffilelist(ifile).name(iS+4:end);
    if isempty(sCh)
        header(end+1, 1:2) = {'hw_channel' '1'}; % ch 1 saved without a suffix
    else
        header(end+1, 1:2) = {'hw_channel' sCh};
    end
    mtint.egf(ifile).header = header;
end
% ---------------- end read egf ---------------------

% ---------------- start read stm ---------------------
% list all stm files
%inpfilelist = dir([filepath,flnmroot, '.stm*']);
% read them in and assign to structure
%for ifile = 1:numel(inpfilelist)
%    [ts, type, value, keyval, Fs] = getstm([filepath, inpfilelist(ifile).name]);
%end


% ---------------- start read inp ---------------------
% list all inp files
inpfilelist = dir([filepath,flnmroot, '.inp*']);
% read them in and assign to structure
for ifile = 1:numel(inpfilelist)
    [ts, type, value, keyval, Fs] = getinp([filepath, inpfilelist(ifile).name]);
    % Parse input/output data as a matrix:
    %   [time channel value]
    % Channel can be input (e.g. I5) or output (O16)

    %% Iterate over all events
    % Create array with [timestamp channel state]
    % (non-changing states will be filtered below)
    mIO = [];
    for t = 1:length(ts)
        % Get status by channels
        vStatus = fliplr(dec2bin(value(t), 16));

        % Determine type (output = 0, input = 1)
        switch type(t)
            case 'O'
                nType = 0;
            case 'I'
                nType = 1;
            otherwise
                nType = -1;
        end
        
        % Iterate over all channels
        for ch = 1:16
            mIO = [mIO; ts(t) ch str2num(vStatus(ch)) nType];
        end
    end
    
    % Remove events that occurred at time zero
    mIO(mIO(:, 1) == 0, :) = [];

    % Assign inputs and outputs
    mOutputs = mIO(mIO(:, 4) == 0, :);
    mInputs = mIO(mIO(:, 4) == 1, :);
    %%
    
    % Remove channels that are never enabled
    for ch = 1:16
        % Outputs
        iCh = mOutputs(:, 2) == ch;
        if ~any(mOutputs(iCh, 3))
            mOutputs(iCh, :) = [];
        end
        % Inputs
        iCh = mInputs(:, 2) == ch;
        if ~any(mInputs(iCh, 3))
            mInputs(iCh, :) = [];
        end
    end
    
    % Remove non-changing states, channel by channel
    for ch = 1:16
        % Outputs
        iIndx = mOutputs(:, 2) == ch;
        if any(iIndx)
            mOutTemp = mOutputs(iIndx, :);
            indxkeep = [1; find(diff(mOutTemp(:, 3))) + 1];
            mOutputs(iIndx, :) = [];
            mOutputs = [mOutputs; mOutTemp(indxkeep, :)];
        end
        
        % Inputs
        iIndx = mInputs(:, 2) == ch;
        if any(iIndx)
            mInpTemp = mInputs(iIndx, :);
            indxkeep = [1; find(diff(mInpTemp(:, 3))) + 1];
            mInputs(iIndx, :) = [];
            mInputs = [mInputs; mInpTemp(indxkeep, :)];
        end
    end
    
    % Sort in temporal order
    [~, indxord] = sort(mOutputs(:, 1));
    mOutputs = mOutputs(indxord, :);
    [~, indxord] = sort(mInputs(:, 1));
    mInputs = mInputs(indxord, :);
    
    %%
    mtint.inp(ifile).Outputs = mOutputs;
    mtint.inp(ifile).Inputs = mInputs;
    mtint.inp(ifile).Fs = Fs;

    % get inp header info
    header = getDACQHeader ( [filepath,egffilelist(ifile).name], 'inp' );
    mtint.inp(ifile).header = header;
end
% ---------------- end read inp ---------------------


[ mtint ] = postprocess_DACQ_data( mtint );

