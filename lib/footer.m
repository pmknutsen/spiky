function y = footer(x, fontsize)
%FOOTER	Puts a text string on bottom of the current figure
%	Y = FOOTER(X, FONTSIZE)
%
%	X: String
%	FONTSIZE: optional, default = 10
%
%	Output: handle to the text object
%	
%	See also: HEADER
%
%	Doron, Sep 5 1995.

if nargin < 2; fontsize = 10; end;

t = findobj(gcf, 'Tag', 'footer');
if length(t) == 0;
    ax = gca;
    axes('Position', [0 0 1 1], 'Visible', 'off');
    t = text(0.5, 0.02, x, ...
        'Units', 'normalized', ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center', ...
        'Tag', 'footer', ...
        'FontSize', fontsize);

    % Bugfix that now allows zoom to work, Per M Knutsen (jan 2005)
    vChildren=get(gcf,'children');
    vChildren = [vChildren(2:end); vChildren(1)];
    set(gcf, 'children', vChildren);
    
    %axes(ax);
else;
    set(t, 'String', x, 'FontSize', fontsize);
end;

if nargout > 0; y = t; end;
