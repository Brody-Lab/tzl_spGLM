function [] = stylizeaxes(varargin)
if nargin < 1
    a = gca;
else
    a = varargin{:};
end
set(a, 'FontSize',18)
set(a, 'ActivePositionProperty', 'outerposition')
set(a, 'Color', 'None')
set(a, 'TickDir', 'out')
set(a, 'NextPlot', 'Add')
set(a, 'LineWidth', 1)
set(a, 'XColor', [0,0,0])
set(a, 'YColor', [0,0,0])
set(a, 'TitleFontWeight', 'normal')
set(a, 'TitleFontSizeMultiplier', 1)
set(a, 'LabelFontSizeMultiplier', 1)
set(a, 'fontname', 'Arial')
set(get(a, 'Parent'), 'Renderer', 'Painters')
