% Matlab theme for Spiky

% The syntax of a Spiky theme file is simple. Each object (figure, axes
% etc) is defined by a sub-structure named by its Type string (e.g. 'figure' for
% figure, 'axes' for axes). Each field in the structure defines the value for
% a property. You can set any number of properties, or omit any (in which
% case the Matlab default property value is used). See the default theme
% file for examples. The parent structure should be named p.

% Axes properties
p.axes.color = [1 1 1];
p.axes.xcolor = [0 0 0];
p.axes.ycolor = [0 0 0];
p.axes.tickdir = 'out';

% Figure properties
p.figure.color = [.8 .8 .8];
p.figure.name  = 'Spiky';
p.figure.NumberTitle = 'off';
p.figure.colormap = jet(2^12);

% Text properties (e.g. figure titles)
p.text.color = [0 0 0];

% Lines (typically requested only for high-contrast lines)
p.line.color = [.1 .1 .1];

% uicontrol (checkboxes etc)
p.uicontrol.backgroundcolor = [.8 .8 .8];
p.uicontrol.foregroundColor = [0 0 0];

% uipanel
p.uipanel.backgroundcolor = [.8 .8 .8];
p.uipanel.foregroundcolor = [.8 .8 .8];
p.uipanel.highlightcolor  = [.8 .8 .8];
p.uipanel.shadowcolor     = [.8 .8 .8];

% hggroup (e.g. bar graphs)
p.hggroup.facecolor = [.15 .15 .15];
p.hggroup.edgecolor = [.15 .15 .15];

