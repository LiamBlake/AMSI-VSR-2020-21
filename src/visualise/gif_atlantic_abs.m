% Required directories
addpath visualise

conFigure(11, 4/3);

% Load data
ud = squeeze(ncread("../data/atlantic_daily.nc", "uo"));
vd = squeeze(ncread("../data/atlantic_daily.nc", "vo"));
lat = double(ncread("../data/atlantic_daily.nc", "latitude"));
lon = double(ncread("../data/atlantic_daily.nc", "longitude"));
time = double(ncread("../data/atlantic_daily.nc", "time"));


% Absolute speed
s = (ud.^2 + vd.^2).^(1/2);
xilim = [241 541]; yilim = [121 361];

h = figure;
axis tight manual
filename = '../figures/atlantic_speed.gif';
for t = 1:length(time)
    % Draw plot for y = x.^n
    cbar = colourplot(lon(xilim(1):xilim(2)),lat(yilim(1):yilim(2)),s(xilim(1):xilim(2),yilim(1):yilim(2),t),"Atlantic Absolute Speed", "Longitude", "Latitude", "Speed (m/s)");
    cbar.Limits = [min(s,[],'all') max(s,[],'all')];
    drawnow
    break
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if t == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.01);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.01);
    end
end
