function [ax] = axis_font_sizes(ax,tickLabelFontSize,labelFontSize)
ax.XAxis.FontSize = tickLabelFontSize;
ax.YAxis.FontSize = tickLabelFontSize;
ax.XLabel.FontSize = labelFontSize;
ax.YLabel.FontSize = labelFontSize;