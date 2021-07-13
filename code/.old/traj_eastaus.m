% Create trajectories from East Australia Current velocities

global u v lat lon time_sc

% Load data
% load('../data/eastaus.mat')
% u = permute(uq, [2 1 3]);
% v = permute(vq, [2 1 3]);
% time = tspan;
% lat = imlat; lon = imlon;
% 
temp = squeeze(ncread("../data/eastaus.nc", "to"));
ud = squeeze(ncread("../data/eastaus.nc", "ugo"));
vd = squeeze(ncread("../data/eastaus.nc", "vgo"));

ud = permute(ud, [2 1 3]);
vd = permute(vd, [2 1 3]);

lat = double(ncread("../data/eastaus.nc", "latitude"));
lon = double(ncread("../data/eastaus.nc", "longitude"));
time = double(ncread("../data/eastaus.nc", "time"));



% Set NaN to zero velocity
ud(isnan(ud)) = 0;
vd(isnan(vd)) = 0;

% Time span - scaled to be within 0 and 1
time_sc = 1*(time - min(time))/(max(time) - min(time));
tspan = linspace(0,1,length(time));

% Velocity functions
u = @(t,x,y) interp3(lon,lat,time_sc,ud,x,y,t,'cubic');
v = @(t,x,y) interp3(lon,lat,time_sc,vd,x,y,t,'cubic');

% Region of interest
xilim = [20 80]; yilim = [10 100];



%% Integrate for trajectories
Nx = 150;
Ny = 150;
traj = NaN(Nx,Ny,length(tspan),2);

% Initial conditions - arrange in matrix
[xp,yp] = meshgrid(linspace(lon(xilim(1)), lon(xilim(2)), Nx), linspace(lat(yilim(1)), lat(yilim(2)), Ny));
init = [xp(:)'; yp(:)'];
n = size(init,2);

% Integrate trajectories
[~,solved_trajs] = ode45(@(t,x) eastaus_vel(t,x,n), tspan, init);
traj(:,:,:,1) = reshape(solved_trajs(:,1:2:end)', Nx,Ny,[]);
traj(:,:,:,2) = reshape(solved_trajs(:,2:2:end)', Nx,Ny,[]);


%% Visualise trajectories
if 1
    animate_trajectories(traj,'r.',.1,1);
    
end



%% Calculate fields and save
X = squeeze(traj(:,:,:,1));
Y = squeeze(traj(:,:,:,2));

% Structure for saving data
resol = struct('x_ip', xp, 'y_ip', yp, 'tlim', [0 1], 'Nx', Nx, 'Ny', Ny, 'Nt', length(tspan));
traj_out = struct('X', X, 'Y', Y);

% Calculate Cauchy-Green tensors
C = NaN(Nx,Ny,2,2);
igrid = NaN(Nx,Ny,2);
igrid(:,:,1) = xp';
igrid(:,:,2) = yp';
for i = 2:(Nx-1)
    for j = 2:(Ny-1)
        C(i,j,:,:) = cauchy_green(X, Y, igrid, i, j);
    end
end

% Calculate FTLE field
ftle = FTLE(C(2:(end-1), 2:(end-1),:,:),tspan(end) - tspan(1));

% Calculate vorticity and LAVD
[lavd,w] = LAVD(X,Y,u,v,tspan);

% Collate in struct
vel = struct('u', u, 'v', v);
fields = struct('cauchygreen', C, 'ftle', ftle, 'vort', w, 'lavd', lavd);
eastaus = struct('vel', vel, 'traj', traj_out, 'res', resol, 'fields', fields);




%% Function
function dXdt = eastaus_vel(t,x,n)
    global lon lat time_sc ud vd
    
    % Resized passed values
    x = reshape(x,[],n);

    % Vectorised equations
    dxdt = diag(interp3(lon,lat,time_sc,ud,x(1,:),x(2,:),t, 'cubic'));
    dydt = diag(interp3(lon,lat,time_sc,vd,x(1,:),x(2,:),t, 'cubic'));
    
    % Linearise output
    dXdt = NaN(2*numel(dxdt),1);
    dXdt(1:2:end) = dxdt(:);
    dXdt(2:2:end) = dydt(:);
end
