function [vCont, vTime, nFs] = spikyfilter_notch(vCont, vTime, nFs)
% Notch filter
% 

tic
% Interpolate NaN indices
vNaNIndx = isnan(vCont);
if any(vNaNIndx)
    vCont(vNaNIndx) = interp1(find(~vNaNIndx), vCont(~vNaNIndx), find(vNaNIndx), 'linear', 'extrap');
end

% Design butterworth notch filter
% This filter provides up to 45 dB of attenuation.
tBut = designfilt('bandstopiir', 'FilterOrder',2, ...
    'HalfPowerFrequency1', 49, 'HalfPowerFrequency2', 51, ...
    'DesignMethod', 'butter', 'SampleRate', nFs);

vCont = filtfilt(tBut, double(vCont));

return