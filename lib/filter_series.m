function vSeriesOut = wta_filter_series(vSeriesIn, nFs, nLp)
% FILTER_SERIES Low-pass filter time series data
%
% Syntax: filter_series(A, FS, LP)
%       A   time series data (matrix with series in columns)
%       FS  sampling frequency (Hz)
%       LP  low-pass frequency (Hz)
%
% Written by Per M Knutsen, January 2004

[a,b] = butter(3, double(nLp)/double(nFs/2), 'low');

vSeriesIn = double(vSeriesIn);

vSeriesOut = zeros(size(vSeriesIn))*NaN;
if ~isempty(nLp) & ~(nLp <= 0)
    for i = 1:size(vSeriesIn, 2)
        vFiltSrsIndx = find(~isnan(vSeriesIn(:,i)));
        vFiltSeries = vSeriesIn(vFiltSrsIndx,i);
        if length(vFiltSeries) >= length(a)*3
            vSeriesOut(vFiltSrsIndx,i) = filtfilt(a, b, vFiltSeries);
        else
            vSeriesOut(vFiltSrsIndx,i) = vFiltSeries;
        end
    end
else, vSeriesOut = vSeriesIn; end

return;
