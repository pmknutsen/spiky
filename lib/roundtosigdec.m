function nNum = roundtosigdec(nNum)
% ROUNDTOSIGDEC Round to first significant decimal place
%
% Examples
%   roundtosigdec(0.083)  is rounded to 0.08
%   roundtosigdec(0.16)   is rounded to 0.2
%   roundtosigdec(0.0018) is rounded to 0.002
%
% Per M Knutsen <pmknutsen@gmail.com>, 2014

sNum = num2str(nNum, '%f');
sNum = sNum((findstr(sNum, '.') + 1):end);
iDec = find(sNum ~= '0');

nFactor = 10^iDec(1);
nNum = round(nNum * nFactor) / nFactor;

return
