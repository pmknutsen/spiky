function vCont = intnans(vCont)
% Fill gaps in a vector that contains NaNs by local mirroring and interpolation
%
%

iNaN = isnan(vCont(:));
iTurnsNaN = find([0;diff(iNaN)] == 1) - 1; % start indices of nan segments

if ~isempty(iTurnsNaN)
    % Step 1 - Complete NaN segments mirror images of adjacent data
    
    iNaNDone = find([0;diff(iNaN)] == -1); % end indices of nan segments
    iNaNDone(iNaNDone <= iTurnsNaN(1)) = [];
    iTurnsNaN(iTurnsNaN > iNaNDone(end)) = [];
    iNaNDone = iNaNDone(1:length(iTurnsNaN));

    for i = 1:length(iTurnsNaN)
        nL = ceil((iNaNDone(i)-iTurnsNaN(i)+1)/2);
        if (iTurnsNaN(i)+nL <= length(vCont)) && (iTurnsNaN(i)-nL > 0)
            vCont(iTurnsNaN(i):(iTurnsNaN(i)+nL)) = fliplr(vCont((iTurnsNaN(i)-nL):iTurnsNaN(i))');
        end
        if ((iNaNDone(i)-nL) > 0) && (iNaNDone(i)+nL <= length(vCont))
            vCont((iNaNDone(i)-nL):iNaNDone(i)) = fliplr(vCont(iNaNDone(i):(iNaNDone(i)+nL))');
        end
    end
    
    % Step 2 - Remove remaining NaNs we missed above by interpolation
    iNaNb = isnan(vCont);
    if any(iNaNb)
        vCont(iNaNb) = interp1(find(~iNaNb), vCont(~iNaNb), find(iNaNb), 'linear', 'extrap');
    end
    
end

return
