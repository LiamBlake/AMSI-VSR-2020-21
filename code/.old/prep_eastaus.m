% Load and visualised ocean data taken from Copernicus EU

% Include additional directories
addpath ocean;

% View info on data
dat = ncinfo("../data/eastaus.nc");
temp = squeeze(ncread("../data/eastaus.nc", "to"));
u = squeeze(ncread("../data/eastaus.nc", "ugo"));
v = squeeze(ncread("../data/eastaus.nc", "vgo"));

lat = double(ncread("../data/eastaus.nc", "latitude"));
lon = double(ncread("../data/eastaus.nc", "longitude"));
tspan = double(ncread("../data/eastaus.nc", "time"));

% Calculate absolute speed
abs = (u.^2 + v.^2).^(0.5);


% Region of interest
xilim = [15 65];
yilim = [15 100];


%% Interpolate to improve resolution
% Set NaN velocities to zero
u(isnan(u)) = 0;
v(isnan(v)) = 0;

% Desired region resolution improvement x
Mx = 4;
My = 3;

xmp = linspace(lon(xilim(1)), lon(xilim(2)), Mx*(diff(xilim)) + 1);
ymp = linspace(lat(yilim(1)), lat(yilim(2)), My*(diff(yilim)) + 1);

% Interpolate for meshgrid
[x,y,t] = meshgrid(lon, lat, tspan);
p = @(A) permute(A, [2 1 3]);
x = p(x); y = p(y); t = p(t);

uq = interp3(lat, lon, tspan, u, ymp, xmp, tspan, 'cubic');
vq = interp3(lat, lon, tspan, v, ymp, xmp, tspan, 'cubic');
tempq = interp3(lat, lon, tspan, temp, ymp, xmp, tspan, 'cubic');

% Insert to original data
%insdat = @(v,n,lims) [v(1:(lims(1)-1)); n; v(1:(lims(2)+1));];
%imlon = insdat(lon, xmp', xilim);
%imlat = insdat(lat, ymp', yilim);
imlon = xmp;
imlat = ymp;

% Save as .mat
save('../data/eastaus.mat', 'temp', 'u', 'v', 'uq', 'vq', 'imlon', 'imlat', 'tspan', 'xilim', 'yilim')



%% Heatmap of Temperature
for i = 1:size(temp,3)
    h = heatmap(tempq(:,:,i)');
    h.Colormap = parula;
    h.YDisplayData = flipud(h.YDisplayData);
    h.XDisplayData = fliplr(h.XDisplayData);
    grid off
    pause(.1)
end


%% Heatmap of Speed
for i = 1:size(temp,3)
    h = heatmap(abs(xilim(1):xilim(2),yilim(1):yilim(2),i)');
    h.Colormap = parula;
    h.YDisplayData = flipud(h.YDisplayData);
    h.XDisplayData = fliplr(h.XDisplayData);
    grid off
    pause(.1)
end





