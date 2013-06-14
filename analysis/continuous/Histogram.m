function Histogram(FV)
%
%
global Spiky;

[sCh, ~] = Spiky.SelectChannelNumber(FV.csChannels);

hFig = figure('color', [.2 .2 .2]);
set(hFig, 'name', 'Spiky Histogram', 'NumberTitle', 'off');
vData = FV.tData.(sCh);

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
set(hAx, 'Color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'fontsize', 7)
xlabel('V')
ylabel('Number of samples')
hTit = title(sprintf('%s (%s)', FV.tChannelDescriptions(nCh).sDescription, sCh));
set(hTit, 'FontSize', 8, 'FontWeight', 'bold', 'color', vCol, 'backgroundcolor', [.1 .1 .1], 'interpreter', 'none')

return
