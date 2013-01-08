function [hLine, hFill] = mean_error_plot(vMean, vError, vColor, vX)
% MEAN_AND_ERROR_PLOT Plot mean and error as a filled area
%
% mean_error_plot(M, E, C, X), where
%       M is the mean vector
%       E is the error vector
%       C is the [R G B] color of the error fill
%       X is the x-vector (optional)
%
% E contains both the upper and lower error ranges for each point in M.
%
% Written by Per M Knutsen, May 2004

vMean = reshape(vMean, length(vMean), 1);
vError = reshape(vError, length(vError), 1);

if ~exist('vX')
    vXt = (1:length(vError))';
else
    vXt = vX';
end

vXb = flipud(vXt);
vYt = vMean + vError;
vYb = flipud(vMean - vError);
hFill = fill([vXt(:);vXb(:)], [vYt(:);vYb(:)], vColor, 'EdgeColor', vColor); hold on;
hLine = plot(vXt, vMean, 'k', 'LineWidth', 1); hold off
set(hFill,'FaceAlpha',.5) % transparency
return;
