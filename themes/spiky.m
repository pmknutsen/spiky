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

% Figure properties
p.figure.color = [.2 .2 .2];
p.figure.NumberTitle = 'off';
p.figure.colormap = jet(2^12);

% Text properties (e.g. figure titles)
p.text.color = [.7 .7 .7];

% Lines (typically requested only for high-contrast lines)
p.line.color = [.8 .8 .8];

% uicontrol (checkboxes etc)
p.uicontrol.backgroundcolor = [.2 .2 .2];
p.uicontrol.foregroundColor = [1 1 1];

% uipanel
p.uipanel.backgroundcolor = [.2 .2 .2];
p.uipanel.foregroundcolor = [.2 .2 .2];
p.uipanel.highlightcolor  = [.2 .2 .2];
p.uipanel.shadowcolor     = [.2 .2 .2];

% hggroup (e.g. bar graphs)
p.hggroup.facecolor = [.85 .85 .85];
p.hggroup.edgecolor = [.85 .85 .85];
