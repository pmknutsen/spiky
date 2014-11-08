function [Fineal_Wave_Form, Mean_Spike_Times] = concatenate_spikes(Spike_Times, Wave_Forms, Jitter, Pre_Time, Post_Time, fs)

% FUNCTION DESCRIPTION:
%  This function matches Spike Times from N Channels and then Concatenates
%  N continuous waveforms accordingly
% INPUTS DESCRIPTION:
% Spike_Times - a matrix contains n spike time vectors in cols. each vector 
%               contains M Spike Times in samples. if there's less then M spike 
%               times then -1 should be added until M Peak times reached.
% Wave_Forms - a matrix contains n cols. in each column should be the
%              waveform sample by sample.
% Jitter, Pre_Time & Post_Time - should be in mSec.
% fs  - sampling rate, samples/sec (Hz).
% 
% OUTPUT DESCRIPTION:
%   Fineal_Wave_Form    - all concatenated spikes
%   Mean_Spike_Times    - average spiketimes of concatenated spikes
%
% GENERAL NOTES ON USAGE:
%  1. The function assumes that within the jitter borders given, not more then one 
%     matching Spike time per channel should be found. Otherwise only the first spike 
%     is to be considered.
%  2. When concatenating data from the continuous Waveforms the function
%     assumes that the continuous data starts at least pre time before the first spike 
%     and ends at least post time after the last spike. No checking is made
%     to watch for stepping out of the continuous data.
%

% Spike_Times = [200 , 225 , 400 ,80]

% -= FUNCTION STARTS HERE =-

% initializing variables
Global_Spike_Times = [];
Channels = [];
Mean_Spike_Times = [];
Match_Spikes_Matrix = [];
Match_Channels_Matrix = [];
Match_Matrix = [];
Combined_Wave_Form = [];
Fineal_Wave_Form = [];

% updates Jitter, Pre_Time and Post_Time form mSecs to samples
Jitter = round(10^-3*Jitter*fs);
Pre_Time = round(10^-3*Pre_Time*fs);
Post_Time = round(10^-3*Post_Time*fs);

% creates the Global_Spike_Times vector & the corresponding channels vector
[Max_Num_Spikes,Num_Channels] = size(Spike_Times);
Match_Spikes_Matrix=zeros(Num_Channels,1);
Match_Channels_Matrix=zeros(Num_Channels,1);
for n = 1:Num_Channels
    Global_Spike_Times = [Global_Spike_Times; Spike_Times(:,n)];
    Channels = [Channels, repmat(n, 1, Max_Num_Spikes)];
end

% combine the two vector into a matrix and sort it by ascending
% Global_Spike_Times values
Global_Times_Channels = [Global_Spike_Times, Channels'];
Global_Times_Channels = sortrows(Global_Times_Channels,1);

Num_Spikes = size(Global_Times_Channels,1);

% spike times matching & finding average spike time
m = 1;
hWait = waitbar(0, 'Concatenating spikes ...');
while m <= Num_Spikes
    waitbar((m/Num_Spikes)/4, hWait)
    n = 1;
    if Global_Times_Channels(m,1) > -1 % there's no need to deal with junk data
        
        Sum_Times = Global_Times_Channels(m,1); % sums the first spike time to be averaged
        
        Current = size(Mean_Spike_Times,2)+1; % the Current spike
        
        % Filling the first rows in the match matrices
        Match_Spikes_Matrix(n,Current) = Global_Times_Channels(m,1);
        Match_Channels_Matrix(n,Current) = Global_Times_Channels(m,2);
        
        if m < Num_Spikes
            while (Global_Times_Channels(m+n,1) - Global_Times_Channels(m,1)) <= Jitter && (m+n) <= Num_Spikes % if we have a match
                Sum_Times = Sum_Times + Global_Times_Channels(m+n,1);  % sums the later spike times to be averaged
                
                % Filling the later rows in the match matrix
                Match_Spikes_Matrix(n+1,Current) = Global_Times_Channels(m+n,1);
                Match_Channels_Matrix(n+1,Current) = Global_Times_Channels(m+n,2);
                
                n = n+1;
                if (m+n) > Num_Spikes 
                    break; % to avoid stepping out of the matrix 
                end
            end
        end
        Mean_Spike_Times = [Mean_Spike_Times, round(Sum_Times/n)]; % adds the new average to the Mean_Spike_Times vector
    end
    m = m+n;
end

% arranging Match_Spikes_Matrix ,Match_Channels_Matrix by channels
Num_Spikes = size(Mean_Spike_Times,2);
for m = 1:Num_Spikes
    MC_vec=flipud(Match_Channels_Matrix(:,m));
    n=1;
    while MC_vec(n)== 0
        MC_vec(n)=Num_Channels+1-n;
        n=n+1;
    end
    Match=sortrows([Match_Spikes_Matrix(:,m),flipud(MC_vec)],2);
    Match_Spikes_Matrix(:,m)=Match(:,1);
    Match_Channels_Matrix(:,m)=Match(:,2);        
end
waitbar(2/4, hWait)
for n = Num_Channels:-1:1
    for m = 1:Num_Spikes
        Match = Match_Spikes_Matrix(n,m);
        Match_Spikes_Matrix(n,m) = 0;
        if Match_Channels_Matrix(n,m) < 1, continue, end
        Match_Spikes_Matrix(Match_Channels_Matrix(n,m),m) = Match; 
    end
end
waitbar(3/4, hWait)

% concatenating
for m = 1:Num_Spikes
    for n = 1:Num_Channels  
        Match = Match_Spikes_Matrix(n,m);
        if  Match == 0
            Match = Mean_Spike_Times(m);
        end;
        Combined_Wave_Form = [Combined_Wave_Form, (Wave_Forms(Match-Pre_Time:Match+Post_Time,n))'];
    end;
    Fineal_Wave_Form = [Fineal_Wave_Form; Combined_Wave_Form];
    Combined_Wave_Form = [];
end
waitbar(4/4, hWait)
close(hWait)

end