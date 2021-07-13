function h = xyheatmap(xval,yval,V,xticks,yticks)
%XYHEATMAP Heatmap of scalar field on x-y plane
%   Detailed explanation goes here

h = heatmap(xval,yval,squeeze(V)');
h.Colormap = parula;
grid off



% Labels - always display 10
xlabs = string(xval);
%xlabs(~xtickl) = " ";
ylabs = string(yval);
%Ylabs(~ytickl) = " ";

h.XDisplayLabels = strings(size(xval));
h.YDisplayLabels = strings(size(yval));

h.YDisplayData = flipud(h.YDisplayData);
h.XLabel = "$x$";
h.YLabel = "$y$";

end

