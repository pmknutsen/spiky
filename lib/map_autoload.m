%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To call this function we need to send the path to the directory of the
% files and the Mapfile.exe (should be in the same directory).
% for example :
% my_path = 'D:\MyDocuments\amir_converter\data\'
% map_autoload(my_path)
%
% written by Amir Monovich, July 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map_autoload(my_path)

old_files = dir (strcat(my_path,'*.map'));
global m
m = 0;
hFig = figure;
my_struct.old_files = old_files;
my_struct.my_path = my_path;
set(hFig, 'Tag', 'load_alphamap_files', 'menu', 'none', ...
    'position', [get(0,'PointerLocation')-50 185 27], ...
    'Name', 'Autoload MAP files', 'NumberTitle', 'off')
uicontrol('Style','pushbutton','string','Start',...
    'Position', [5 5 80 20],'Callback',@OnRun, 'Tag', 'run', 'UserData', my_struct)
uicontrol('Style','pushbutton','string','Stop','Position', [100 5 80 20], 'Callback','global m; m = 1;','tag','stop','enable','on');
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = OnRun(varargin)
set(findobj('Tag', 'run'), 'enable', 'off')
set(findobj('Tag', 'stop'), 'enable', 'on')
my_struct = get(findobj(gcf, 'Tag', 'run'), 'UserData');
my_path = my_struct.my_path;
% Path of Mapfile.exe
sMapfilePath = which('Mapfile.exe');
global m
while m ~= 1
    map_files = dir([my_path '*.map']); % list of all .map files
    mat_files = dir([my_path '*.mat']); % list of all .map files

    for n_map = 1:size(map_files,1) % iterate over .map files
        bConvert = 1;
        for n_mat = 1:size(mat_files,1) % iterate over .mat files
            if strcmp(map_files(n_map).name(1:end-4), mat_files(n_mat).name(1:end-4)) % mat exists
                bConvert = 0;
            end
        end
        if bConvert % convert map file unless corresponding mat exists
            % Convert map file
            [status, error] = dos([sMapfilePath ' =' my_path map_files(n_map).name '=Matlab=' my_path '=']);
            % Load data into spiky
            spiky(['LoadTrial(''' map_files(n_map).name(1:end-3) 'mat' ''')']);
            spiky('ViewTrialData')
            figure(findobj('Tag', 'load_alphamap_files'))
            drawnow
            break
        end
    end
    tic; while toc <= 5
        if (m == 1), break, end
        pause(0.1)
    end
    figure(findobj('Tag', 'load_alphamap_files'))
end
hFig = findobj('Tag', 'load_alphamap_files');
close(hFig)
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
