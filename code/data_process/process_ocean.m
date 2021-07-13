% Preprocess Atlantic ocean velocities and temperatures. Produces
% trajectories, LCS diagnostic fields and interpolated temperatures at each
% initial point.

addpath fields visualise

% Globals
NO_DAYS = 20;

% Load velocity and spatiotemperal data
fprintf("Loading data...\n");
ud = squeeze(ncread("../data/atlantic_daily.nc", "uo"));
vd = squeeze(ncread("../data/atlantic_daily.nc", "vo"));
lat = double(ncread("../data/atlantic_daily.nc", "latitude"));
lon = double(ncread("../data/atlantic_daily.nc", "longitude"));
time = double(ncread("../data/atlantic_daily.nc", "time"));
fprintf("Success\n");

% Adjust latitude, longitude for rounding errors to create uniform grid
fprintf("Preprocessing data...\n");
lstep = 1/12;
lat = (lat(1):lstep:lat(end))';
lon = (lon(1):lstep:lon(end))';

% Convert m/s velocities to deg/s
% WGS 84 model of Earth
a = 6378137.0;
f = 1/298.257223563;
e2 = 2*f - f^2;
ulat = @(p) (180*(1 - e2*(sind(p)).^2).^(3/2))./(pi*a*(1-e2));
ulon = @(p) (180*(1 - e2*(sind(p)).^2).^(1/2))./(pi*a*cosd(p));

ud = ud.*repmat(ulon(lat)',length(lon),1,length(time));
vd = vd.*repmat(ulat(lat)',length(lon),1,length(time));

% Set NaN to zero velocity
landu = isnan(ud);
landv = isnan(vd);
ud(landu) = 0;
vd(landv) = 0;

% Time span - shifted to start at 0
time_sc = (time - time(1));
time_sc = time_sc(1:NO_DAYS);
ud = ud(:,:,1:NO_DAYS)*3600;
vd = vd(:,:,1:NO_DAYS)*3600;

% Velocity functions
u_interp = griddedInterpolant({lon,lat,time_sc}, ud,'cubic','none');
v_interp = griddedInterpolant({lon,lat,time_sc}, vd,'cubic','none');

u = @(t,x,y) u_interp(x,y,t);
v = @(t,x,y) v_interp(x,y,t);


% Region of interest
xilim = [241 541]; yilim = [121 361];
tlim = [0 time_sc(NO_DAYS)];
fprintf("Success\n");


%% Integrate for trajectories
Nx = 300;
Ny = 150;
Nt = 100;
tspan = linspace(0,time_sc(end),Nt);

% Initial conditions - arrange in matrix
[xp,yp] = meshgrid(linspace(lon(xilim(1)), lon(xilim(2)), Nx), linspace(lat(yilim(1)), lat(yilim(2)), Ny));
xp = xp'; yp = yp';
init = [xp(:); yp(:)];
n = Nx*Ny;

% Integrate trajectories
fprintf("Integrating for trajectories...\n");
[~,solved_trajs] = ode45(@(t,x) atlantic_vel(t,x,n,u,v), tspan, init);

traj = NaN(Nx,Ny,Nt,2);
traj(:,:,:,1) = reshape(solved_trajs(:,1:end/2)',Nx,Ny,Nt);
traj(:,:,:,2) = reshape(solved_trajs(:,end/2+1:end)',Nx,Ny,Nt);
fprintf("Success\n");

%% Visualise trajectories
if 0
    animate_trajectories(traj,'r.',.1,1);
end


%% Temperature interpolation
% SEEMS TO BE FIXED WHEN ANIMATING, NOT SURE WHAT TO DO
fprintf("Interpolating temperature field...\n");
tempd = squeeze(ncread("../data/atlantic_daily.nc", "thetao"));

% Set NaN (land) to -inf
tempd(isnan(tempd)) = -inf;

temp_interp = griddedInterpolant({lon,lat,time_sc}, tempd(:,:,1:NO_DAYS),'cubic','none');
temp_f = @(t,x,y) temp_interp(x,y,t);

fodder = NaN(Nx,Ny,Nt);
for i = 1:Nx
   for j = 1:Ny
       fodder(i,j,:) = tspan;
   end
end

temp_vec = temp_f(fodder,traj(:,:,:,1),traj(:,:,:,2));
fprintf("Success\n");


%% Calculate fields and save
X = squeeze(traj(:,:,:,1));
Y = squeeze(traj(:,:,:,2));

% Calculate Cauchy-Green tensors
fprintf("Calculating Cauchy-Green tensors...\n")
C = NaN(Nx,Ny,2,2);
for i = 2:(Nx-1)
    for j = 2:(Ny-1)
        C(i,j,:,:) = cauchy_green(X, Y, i, j);
    end
end
fprintf("Success\n");

% Calculate FTLE field
fprintf("Calculating FTLE field...\n");
ftle = FTLE(C(2:(end-1),2:(end-1),:,:),tspan(end) - tspan(1));
fprintf("Success\n");

% Calculate vorticity and LAVD
fprintf("Calculating vorticity and LAVD field...\n");
[lavd,w] = LAVD(X,Y,u,v,tspan);
fprintf("Success\n");


%% Output data
% Collate in structs
resol = struct('x_ip', xp, 'y_ip', yp, 'Nx', Nx, 'Ny', Ny, ...
               'Nt', length(tspan), 'unan', landu, 'vnan', landv, ...
               'lon', lon, 'lat', lat, 'tlim', tlim, 'no_days', NO_DAYS);
vel = struct('u', u, 'v', v);
fields = struct('cauchygreen', C, 'ftle', ftle, 'vort', w, 'lavd', lavd);
ocean = struct('vel', vel, 'res', resol);
temp_ex = struct('temp', temp_vec);%,'temp_f', temp_f);

% Save data
fprintf("Saving data...\n");
save('../data/atlantic_res.mat', 'ocean');
save('../data/atlantic_traj_X.mat', 'X');
save('../data/atlantic_traj_Y.mat', 'Y');
save('../data/atlantic_fields.mat', 'fields');
save('../data/atlantic_temp.mat', 'temp_ex');
fprintf("Success\n");
fprintf("Processing complete!\n");


%% ODE Function
function dXdt = atlantic_vel(t,x,n,u,v)
    
    ts = t*ones(n,1);

    dXdt = zeros(2*n,1);
    dXdt(1:n,1) = u(ts,x(1:n,1),x(n+1:2*n,1));
    dXdt(n + 1:2*n,1) = v(ts,x(1:n,1),x(n+1:2*n,1));
    
end
