%% Load data
temp = squeeze(ncread("../data/eastaus.nc", "to"));
u = squeeze(ncread("../data/eastaus.nc", "ugo"));
v = squeeze(ncread("../data/eastaus.nc", "vgo"));

lat = double(ncread("../data/eastaus.nc", "latitude"));
lon = double(ncread("../data/eastaus.nc", "longitude"));
tspan = double(ncread("../data/eastaus.nc", "time"));



%% Plots
% Temperature
if 1
    figure;
    for i = 1:size(temp,3)
        h = heatmap(temp(:,:,i)');
        h.Colormap = parula;
        h.YDisplayData = flipud(h.YDisplayData);
        h.XDisplayData = fliplr(h.XDisplayData);
        grid off
        pause(.1)
    end
end


%% Calculate fields





%% Extract LAVD regions


%% Cluster trajectories and temperature



% Store temperature as concatenated vector
tempc = reshape(temp, size(temp,1)*size(temp,2), []);
tempc(isnan(tempc)) = -100;
[~,U_temp] = fcm(tempc, 3);

