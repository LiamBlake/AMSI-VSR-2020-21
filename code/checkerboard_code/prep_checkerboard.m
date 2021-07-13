global ui vi

% Load frames
[u,v,x0,y0,ts] = load_checkerboard();

% Trim 50 pixels from each boundary - account for unreliable measurements
xilim = [50 size(u,1) - 50];
yilim = [50 size(u,2) - 50];

% Velocity functions
ui = @(t,x,y) diag(interp3(x0(1,:),y0(:,1)',ts,u,x,y,t));
vi = @(t,x,y) diag(interp3(x0(1,:),y0(:,1)',ts,v,x,y,t));

%% Integrate trajectories

Nx = 5;
Ny = 5;
Nt = 50;
tspan = linspace(0,ts(end),Nt);
traj = NaN(Nx,Ny,length(tspan),2);

% Initial conditions - arrange in matrix
[xp,yp] = meshgrid(linspace(x0(xilim(1),yilim(1)), x0(xilim(2), yilim(2)), Nx), linspace(y0(xilim(1), yilim(1)), y0(xilim(2), yilim(2)), Ny));
init = [xp(:)'; yp(:)'];
n = size(init,2);

%% Integrate trajectories
[~,solved_trajs] = ode45(@(t,x) cboard_ode(t,x,n), tspan, init);
traj(:,:,:,1) = reshape(solved_trajs(:,1:2:end)', Nx, Ny,[]);
traj(:,:,:,2) = reshape(solved_trajs(:,2:2:end)', Nx, Ny,[]);


%% Visualise trajectories
if 1
    animate_trajectories(traj,'r.',.1,1);
end

%% Fields
X = squeeze(traj(:,:,:,1));
Y = squeeze(traj(:,:,:,2));

% Structure for saving data
resol = struct('x_ip', xp, 'y_ip', yp, 'tlim', [0 1], 'Nx', Nx, 'Ny', Ny, 'Nt', Nt);
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
cboard = struct('vel', vel, 'res', resol, 'fields', fields);

% Save data
save('../data/checkerboard_traj.mat', 'traj_out');
save('../data/checkerboard_fields.mat', 'cboard');



%% Function
function dXdt = cboard_ode(t,x,n)
    global ui vi
    
    % Resized passed values
    x = reshape(x,[],n);

    % Vectorised equations
    dxdt = ui(t,x(1,:),x(2,:));
    dydt = vi(t,x(1,:),x(2,:));
    
    % Linearise output
    dXdt = NaN(2*numel(dxdt),1);
    dXdt(1:2:end) = dxdt(:);
    dXdt(2:2:end) = dydt(:);
    
end


