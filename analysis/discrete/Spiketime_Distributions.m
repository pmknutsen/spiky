function Spiketime_Distributions(FV)
%
%
%

global Spiky;

% Plot rasters of sorted units
hFig = figure; % create figure
set(hFig, 'color', [.2 .2 .2], 'Name', 'Spiky Spiketime Distributions', 'NumberTitle', 'off', 'menubar', 'none')
hAx = axes('position', [0.1 .075 .88 .9]);
set(hAx, 'color', [.1 .1 .1], 'xcolor', [.6 .6 .6], 'ycolor', [.6 .6 .6], 'fontsize', 8)
hold on
% Iterate over channels
csChannels = fieldnames(FV.tSpikes);
nRow = 0;
sUnits = {}; nU_max = 0;
for nCh = 1:length(csChannels)
    % Iterate over units
    tSpikes = FV.tSpikes.(csChannels{nCh});
    if ~isfield(tSpikes, 'hierarchy'), continue, end
    vUnits = unique(tSpikes.hierarchy.assigns);
    nFs = tSpikes.Fs(1);
    mBoxes = []; mCols = [];
    for nU = 1:length(vUnits)
        nRow = nRow + 1;
        vIndx = find(tSpikes.hierarchy.assigns == vUnits(nU));
        vSpiketimes = tSpikes.spiketimes(vIndx); % samples
        vSpiketimes = vSpiketimes ./ nFs; % sec
        mBoxes = [mBoxes; vSpiketimes repmat(nU, length(vSpiketimes), 1)];
        if vUnits(nU) == 0, mCols(end+1, :) = [.4 .4 .4];
        else mCols(end+1, :) = FV.mColors(nU,:); end
        %plot(vSpiketimes, repmat(nRow,length(vSpiketimes),1), '.', 'color', mCols(end,:), 'linewidth', .5);
        sUnits{end+1} = [csChannels{nCh} '_' num2str(vUnits(nU))];
    end
    boxplot(mBoxes(:,1), mBoxes(:,2), 'orientation', 'horizontal', 'colors', mCols, 'positions', nU_max+(1:nU))
    nU_max = nU_max + mBoxes(end,2);
end
if nRow == 0 % if no units have been sorted
    close(hFig)
    uiwait(warndlg('No units have been sorted.'))
    return
end
xlabel('Time (sec)')
set(hAx, 'ylim', [0 nRow+1], 'ytick', 1:nRow, 'yticklabel', sUnits, 'ygrid', 'off', 'xgrid', 'on')

return
