% Load and visualised ocean data taken from Copernicus EU

% Include additional directories
addpath ocean;

% View info on data
dat = ncinfo("../data/eastaus.nc");
temp = squeeze(ncread("../data/eastaus.nc", "to"));
u = squeeze(ncread("../data/eastaus.nc", "ugo"));
v = squeeze(ncread("../data/eastaus.nc", "vgo"));

% Calculate absolute speed
abs = (u.^2 + v.^2).^(0.5);

%% Heatmap of Temperature
for i = 1:size(temp,3)
    h = heatmap(temp(:,:,i)');
    h.Colormap = parula;
    h.YDisplayData = flipud(h.YDisplayData);
    h.XDisplayData = fliplr(h.XDisplayData);
    grid off
    pause(.1)
end


%% Heatmap of Speed
for i = 1:size(temp,3)
    h = heatmap(abs(:,:,i)');
    h.Colormap = parula;
    h.YDisplayData = flipud(h.YDisplayData);
    h.XDisplayData = fliplr(h.XDisplayData);
    grid off
    pause(.1)
end



%% Interpolate to improve resolution
% Set NaN velocities to zero
u(isnan(u)) = 0;
v(isnan(v)) = 0;

