function mtint = readAllDACQdata ( varargin )

% cd('MY_DATA_DIRECTORY') % Replace this line with the appropriate string

if nargin == 0
    [filename,filepath] = uigetfile('*.set','Select the .set file...',...
        'MultiSelect','on');
elseif nargin == 2
    filepath = varargin{1};
    filename = varargin{2};
end

% check all required files are present
[fileStruct, tetsAvailable] = list_and_check_DACQFiles( filepath,cellstr(filename) );

if numel(fileStruct) > 1
    % process headers
    posHeaders = headerCheck(filepath,cellstr(fileStruct(1).flnmroot),'pos');
    pos = struct('led_pos',{},'ts',{},'led_pix',{},'header',{});
    pos(1).header = posHeaders;
    % for now use only the first trials .set header
    mtint.header = getDACQHeader ( [filepath,fileStruct(1).flnmroot,'.set'], 'set' );
    % load pos
    idx = find(strcmpi('colactive_1',mtint.header(:,1)));
    n_leds = numel(find(str2double(char(mtint.header(idx:idx+3,2)))));
    for ifile = 1:numel(fileStruct)
        if ifile == 1
            duration = 0;
        else
            current_duration = key_value('duration',getDACQHeader([filepath,fileStruct(ifile).flnmroot,'.pos'],'pos'),'num');
            duration = duration + current_duration;
        end
        [led_pos,post,led_pix] = rawpos([filepath,fileStruct(ifile).flnmroot,'.pos'],n_leds); % raw led positions
        pos.led_pos = [pos.led_pos;led_pos];
        pos.led_pix = [pos.led_pix;led_pix];
        ts = post + duration;
        pos.ts = [pos.ts;ts];
    end
    clear duration current_duration
    mtint.pos = pos;

    % load eeg
    % check all eeg files are present first
    for i = 1:numel(fileStruct)
        eegPresent = fileStruct(i).eegFiles;
    end
    if all(eegPresent)
        [EEG,Fs] = geteeg([filepath,fileStruct(1).flnmroot,'.eeg']);
        eeg.eeg = EEG;
        eeg.Fs = Fs;
        clear EEG Fs
        for ifile = 2:numel(fileStruct)
            [EEG,Fs] = geteeg([filepath,fileStruct(ifile).flnmroot,'.eeg']);
            eeg_out.eeg = EEG;
            eeg_out.Fs = Fs;
            clear EEG Fs
            eeg_out.eeg = [eeg.eeg;eeg_out.eeg];
            eeg_out.Fs = [eeg.Fs;eeg_out.Fs];
        end
    end
    mtint.eeg = eeg;
        
    % use the first file as a header
    mtint.flnmroot = filelist;
    mtint.filepath = filepath;
    mtint.header = getDACQHeader([filepath,filelist{1},'.set'],'set');
    mtint.pos.header = headerCheck(filepath,filelist,'pos');
    mtint = postprocess_DACQ_data( mtint );
else
    mtint = readDACQdata ( filepath, fileStruct.flnmroot );
end

% do a final check on the positional data to make sure that the window max
% and min values are enough to contain the x/y coordinate data. if not then
% % scale the data appropriately and inform the user this has happened.
% mtint = checkPosScaling (mtint);

