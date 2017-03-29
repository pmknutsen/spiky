function [vCont, vTime, nFs] = spiky_filter_rectify(vCont, vTime, nFs)
% Invert signal
% 

vCont = abs(vCont);

return