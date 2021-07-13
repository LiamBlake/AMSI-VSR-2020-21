function plot_memberships(mu, suptitle,xval, yval, xlab, ylab)
%PLOT_MEMBERSHIPS 

% Default arguments
if nargin < 5
    xlab = "$x$";
    ylab = "$y$";
end

smu = squeeze(mu);

K = size(smu,3);
figure;
t = tiledlayout(ceil(K/4),4);

for k = 1:K
    nexttile; colourplot(xval, yval, smu(:,:,k), sprintf('Cluster %i', k), xlab, ylab);
    colorbar('off');
    xlabel([]);
    ylabel([]);
    
end

% Shared colorbar
v = regexp(version,'(?<=\().*(?=\))','match');      % Check MATLAB version
v = v{1};
if strcmp(v(1:end-1), '2020')
    c = colorbar;
    c.Layout.Tile = "east";
else
    c = colorbar("EastOutside");
end
c.Label.String = "Membership Value";

% Global labels
title(t, suptitle, 'Interpreter', 'LaTeX')

end