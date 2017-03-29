function [vCont, vTime, nFs] = spikyfilter_subtract_mean(vCont, vTime, nFs)
% Subtract mean from the signal
% 

vCont = vCont - nanmean(vCont);

return