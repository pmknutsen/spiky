function y = header(x, fontsize)
%HEADER Puts a text string on top of the current figure
%	Y = HEADER(X, FONTSIZE)
%
%	X: String
%	FONTSIZE: optional, default = 20
%
%	Output: handle to the text object
%
%	See also: FOOTER
%
%	Doron, Sep 5 1995.

if nargin < 2; fontsize = 20; end;

t = findobj(gcf, 'Tag', 'header');
if length(t) == 0;
    ax = gca;
    axes('Position', [0 0 1 1], 'Visible', 'off');
    t = text(0.5, 0.98, x, ...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'center', ...
        'Tag', 'header', ...
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

