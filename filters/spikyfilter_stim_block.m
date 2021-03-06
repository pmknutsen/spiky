function [vCont, vTime, nFs] = spikyfilter_stim_block(vCont, vTime, nFs)
% Remove stimulus artefact through blocking. Requires a digital trigger
% event to be defined.
% 

global Spiky
persistent p_sTrigCh p_nStartDelay p_nStopDelay p_nSmooth

%% Select trigger channel
[FV, ~] = Spiky.main.GetStruct();
if isempty(p_sTrigCh) || isempty(vCont)
    [p_sTrigCh, ~] = Spiky.main.SelectChannelNumber(FV.csDigitalChannels, 'Select stimulus channel', p_sTrigCh);
    if isempty(p_nStartDelay)
        p_nStartDelay = 0; % ms
        p_nStopDelay = 0; % ms
        p_nSmooth = 10; % ms
    end
    csAns = inputdlg({'Onset pad (ms)', 'Offset pad (ms)', 'Smoothing (% of window)'}, ...
        'Stimulus block', 1, ...
        {num2str(p_nStartDelay), num2str(p_nStopDelay), num2str(p_nSmooth)} );
    if isempty(csAns), return, end
    p_nStartDelay = str2num(csAns{1});
    p_nStopDelay = str2num(csAns{2});
    p_nSmooth = max(1, min(100, str2num(csAns{3})));
end
if isempty(vCont), return; end

% Get trigger up/down times
if ~isfield(FV.tData, [p_sTrigCh '_Up']), return; end

%% Event timestamps
[vUp, vDown] = Spiky.main.GetEventPairs(p_sTrigCh, 'all');
vUp = vUp - (p_nStartDelay / 1000);
vDown = vDown + (p_nStopDelay / 1000);
vUp = interp1(vTime, 1:length(vTime), vUp);
vDown = interp1(vTime, 1:length(vTime), vDown);
vUp = floor(vUp(~isnan(vUp))); % samples
vDown = floor(vDown(~isnan(vDown))); % samples

%% Block stimulus artefact
for e = 1:length(vUp)
    vRange = vUp(e):vDown(e);
    vRangePre = (vUp(e) - length(vRange)):(vUp(e) - 1);
    if vRangePre(1) > 1
        vCont(vRange) = fliplr(vCont(vRangePre));
    end
    % Smooth inserted trace through a moving average
    vSmoothRange = vRangePre(1):(vRange(end)+length(vRangePre));
    vCont(vSmoothRange) = smooth(vCont(vSmoothRange), length(vRange)*(p_nSmooth/100), 'moving');
end

return
