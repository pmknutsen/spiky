function [vCont, vTime, nFs] = spiky_filter_invert(vCont, vTime, nFs)
% Invert signal
% 

vCont = vCont .* -1;

return