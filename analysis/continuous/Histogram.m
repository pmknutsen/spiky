function Histogram(FV)
% Histogram of data points
%
global Spiky

[sCh, ~] = Spiky.main.SelectChannelNumber(FV.csChannels);

vData = FV.tData.(sCh);
if isempty(vData), return; end

hFig = figure('name', 'Spiky Histogram', 'NumberTitle', 'off');
Spiky.main.ThemeObject(hFig)

% Remove extreme outliers (0.025% off either edge)
vData = sort(vData);
vData = vData((length(vData)/400):(length(vData)-(length(vData)/400)));
[vN, vX] = hist(vData, 100);
hAx = axes;
nCh = strcmpi({FV.tChannelDescriptions.sChannel}, sCh);
vCol = FV.mColors(nCh,:);

if isempty(vCol)
    vCol = [.8 .8 .8]; % default color, gray
end

hBar = bar(vX, vN);
set(hBar, 'faceColor', vCol, 'edgeColor', vCol)
Spiky.main.ThemeObject(hAx)


% Get unit
if isfield(FV.tData, [sCh '_Unit'])
    sUnit = FV.tData.([sCh '_Unit']);
else
    sUnit = FV.tAmplitudeUnit.sUnit;
end
xlabel(sprintf('Sample value (%s)', sUnit))
ylabel('Number of samples')

hTit = title(sprintf('%s  %s', sCh, FV.tChannelDescriptions(nCh).sDescription));
Spiky.main.ThemeObject(hTit, 'interpreter', 'none')
axis(hAx, 'tight')

return
