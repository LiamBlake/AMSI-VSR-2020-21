function c = colourplot(x,y,C,stitle,xlab,ylab,clab,cmap)
%COLORPLOT Visualises 2D scalar field with color plot
%   Detailed explanation goes here

% Default arguments
if nargin < 8
    cmap = parula;
end
if nargin < 7
    clab = "";
end
if nargin < 5
    xlab = "";
    ylab = "";
end
if nargin < 4
    stitle = "";
end

h = imagesc(x([1 end]), y([1 end]),squeeze(C)');
set(gca,'YDir','normal')
c = colorbar;
colormap(cmap);
set(h, 'AlphaData', ~isnan(squeeze(C)'))
%colormap([0 0 0; parula(256)])

title(stitle,'Interpreter', 'LaTeX');
xlabel(xlab,'Interpreter', 'LaTeX');
ylabel(ylab,'Interpreter', 'LaTeX');
c.Label.String = clab;


end

