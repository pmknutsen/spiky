function Welchs_Power_Spectrum(FV)
% Compute power spectral density estimates via Welch's method (pwelch)
%
% Spectra of all analog channels are computed and displayed.
%
% Requirements:
%   Signal Processing Toolbox
%

global Spiky g_bBatchMode

% Get channel
persistent p_sContCh
if isempty(p_sContCh) || ~g_bBatchMode
    [p_sContCh, bResult] = Spiky.main.SelectChannelNumber(FV.csDisplayChannels, 'Select continuous signal', p_sContCh);
    if ~bResult; return, end
end
drawnow

%% Fetch data
vCont = double(FV.tData.(p_sContCh)');
if all(size(vCont) > 1); return; end
nFs = FV.tData.([p_sContCh '_KHz']) * 1000;
nDur = length(vCont) / nFs; % signal duration, s

% Get parameters interactively
persistent p_nMinF p_nMaxF p_nD
if isempty(p_nMinF) || ~g_bBatchMode
    if isempty(p_nMinF), p_nMinF = 0; end
    if isempty(p_nMaxF), p_nMaxF = nFs/2; end
    if isempty(p_nD), p_nD = 1; end
    
    cPrompt = {'Min frequency (Hz):', 'Max frequency (Hz):', 'Signal derivative:'};
    cAnswer = inputdlg(cPrompt, 'Options', 1, {num2str(p_nMinF), num2str(p_nMaxF), num2str(p_nD), });
    if isempty(cAnswer), return, end
    p_nMinF = str2num(cAnswer{1}); % hz
    p_nMaxF = str2num(cAnswer{2}); % hz
    p_nD = str2num(cAnswer{3});
end

p_nMinF = max(0, p_nMinF);
p_nMaxF = min(nFs, p_nMaxF);

% Get channel descriptive string
sDescr = Spiky.main.GetChannelDescription(p_sContCh);
if isempty(sDescr); sDescr = p_sContCh; end

% Compute derivatives
if p_nD > 0
    vCont = diff(vCont, p_nD);
end

% Fill in NaN's with nearest non-NaN neighbour
vCont(isnan(vCont)) = interp1(find(~isnan(vCont)), vCont(~isnan(vCont)), find(isnan(vCont)), 'linear');

%[vPxx, vF] = pwelch(vCont, kaiser(length(vCont),4), 0, 2^16, nFs);
[vPxx, vF] = pwelch(vCont, [], [], [], nFs);

% Plot
hFig = figure;
set(hFig, 'name', 'Welch Power Spectrum')
Spiky.main.ThemeObject(hFig);
hAx = axes('position', [.1 .125 .85 .8], 'xscale', 'li', 'yscale', 'log');
Spiky.main.ThemeObject(hAx);
hold(hAx, 'on')
axes(hAx);
box(hAx, 'on')
grid(hAx, 'on')
set(hAx, 'xlim', [p_nMinF p_nMaxF])

hPw = plot(hAx, vF, 10*(vPxx));

xlabel('Frequency (Hz)')
ylabel('Spectral Density (Hz^{-1})')
hTit = title(sprintf('Welch Power Spectrum:  %s', p_sContCh), 'interpreter', 'none');
Spiky.main.ThemeObject(hTit);
Spiky.main.AttachAxisCrossHairs(hAx)

%%


return

