function [mData_dejit, vSNDX] = dejitter_spikes(mData, nContact, nRange, nDimMean)
% DEJITTER_SPIKES Dejitter spikes using algorithm of Aldworth et al (2005).
%
% [mData_dejit, vSNDX] = DEJITTER_SPIKES(mData, nContact, nRange, nDimMean)
%
% Inputs:
%   mData       Matrix containg spike wavewforms (time vs volt)
%   nContact    Index in mData where trigger occurred (i.e. time zero)
%   nRange      Maximal allowed shift of spikes (datapoints)
%   nDimMean    Range after trig to include in dejittering included (datapoints)
%
% Outputs:
%   mData_dejit The dejittered version of mData. Contains NaNs due to shifts
%   vSNDX       Shifts of dejiterred traces (ms)
%
% Dependencies:
%   dejitterEM.m by Alex Dimitrov
%   NETLAB toolbox
%
% Written by Per M Knutsen, July 2006
%

[vI, vJ] = find(isnan(mData)); % remove trials with NaNs
mData(:,vJ) = [];

% RESAMPLE TO 2000 FPS
nResampFact = 4;
mData = ResampleMData(mData, size(mData,1)*nResampFact);
mData = fliplr(mData');

% DEJITTER
nContact = size(mData,2)-(nContact-1)*nResampFact; % frame -> interp. frame
nDimMean = nDimMean * nResampFact; % ms -> interp. frame
tGMM = gmm(nDimMean, 1, 'spherical'); % gaussian model on which the stimulus is aligned
tGMM = gmminit(tGMM, reshape(mData(:,(nContact-nDimMean+1):nContact), size(mData,1), tGMM.nin), []);
nRange = nRange * nResampFact; % # ms -> interp. frames
nS = 10; % initial guess for jitter variance (interp. frames)
nN = 1000; % max number of iterations
nRelDiff = 10; % stopping criterion. stop if |s_n - s_n-1|/s_n-1 < rel_dif
[nm, vSNDX, newstim, Ss, means] = dejitterEM(mData, tGMM, nContact, nRange, nS, nN, nRelDiff);

mData_dejit = [];
nLen = size(mData,2);
for d = 1:length(vSNDX)
    nShift = round(vSNDX(d));
    if nShift > 0
        vInsert = mData(d, abs(nShift):end);
        vInsert = [vInsert zeros(1,abs(nShift)-1)*NaN];
    elseif nShift < 0
        vInsert = mData(d, 1:(end-abs(nShift)));
        vInsert = [zeros(1,abs(nShift))*NaN vInsert];
    else
        vInsert = mData(d, :);
    end
    mData_dejit(d,:) = vInsert;
end
mData_dejit = fliplr(mData_dejit);

% RESAMPLE TO 500 FPS
%mData_dejit = ResampleMData(mData_dejit', size(mData_dejit,2)/4);

%vSNDX = vSNDX / 2; % shifts (ms)

return

% Resample mData traces to length nLen
% Dimensions of mData must be [FRAME TRACE]
function mData = ResampleMData(mData, nLen)
mData_new = [];
for j = 1:size(mData, 2)
    mData_new(:,j) = interp1(1:size(mData,1), mData(:,j), linspace(1,size(mData,1),nLen), 'spline')';
end
mData = mData_new;
return
