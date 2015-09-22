% Default theme for Spiky

% The syntax of a Spiky theme file is simple. Each object (figure, axes
% etc) is defined by a sub-structure named by its Type string (e.g. 'figure' for
% figure, 'axes' for axes). Each field in the structure defines the value for
% a property. You can set any number of properties, or omit any (in which
% case the Matlab default property value is used). See the default theme
% file for examples. The parent structure should be named p.

% Axes properties
p.axes.color = [.15 .15 .15];
p.axes.xcolor = [.7 .7 .7];
p.axes.ycolor = [.7 .7 .7];
p.axes.tickdir = 'out';
p.axes.fontsize = 8;
p.axes.ticklength = [.005 .005];

mGradient = uint8([ 20 20 50 50; ... % [B B T T] red
    20 20 50 50; ... % green
    20 20 50 50; ... % blue
    255 255 255 255]);
p.backdrop.face.ColorData = mGradient;

% Figure properties
p.figure.color = [.2 .2 .2];
p.figure.NumberTitle = 'off';
if ispc
    p.figure.colormap = jet(2^7);
else
    p.figure.colormap = jet(2^12);
end
    
% Text properties (e.g. figure titles)
p.text.color = [.7 .7 .7];

% Lines (typically requested only for high-contrast lines)
p.line.color = [.8 .8 .8];

% uicontrol (checkboxes etc)
p.uicontrol.backgroundcolor = p.figure.color;
p.uicontrol.foregroundColor = [1 1 1];

% uipanel
p.uipanel.backgroundcolor = p.figure.color;
p.uipanel.foregroundcolor = p.figure.color;
p.uipanel.highlightcolor  = p.figure.color;
p.uipanel.shadowcolor     = p.figure.color;

% hggroup (e.g. bar graphs)
p.hggroup.facecolor = [.85 .85 .85];
p.hggroup.edgecolor = [.85 .85 .85];

% colors on lines
p.colors = distinguishable_colors(100, p.axes.color);
