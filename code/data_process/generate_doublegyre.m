% Generate and saves trajectories for double gyre flow. Also calculates
% FTLE, Cauchy-Green tensors, vorticity and LAVD fields for the flow. Data
% is saved to a .mat file and stored in a struct.

% Resolutions
n = 200;
m = 100;
tn = 100;
tlim = [0 1];
xlim = [0 2];
ylim = [0 1];
xr = linspace(xlim(1),xlim(2),n);
yr = linspace(ylim(1),ylim(2),m);
tr = linspace(tlim(1),tlim(2),tn);
[x, y] = meshgrid(xr, yr);

% Constants
A = 1;
eps = 0.3;
w = 10*pi;

% Functions
h = @(x,t) eps*sin(w*t).*x.^2 + (1 - 2*eps*sin(w*t)).*x;
dhdx = @(x,t) 1 - 2*eps*sin(w*t).*(1 - x);
u = @(t,x,y) -pi*A*sin(pi*h(x,t)).*cos(pi*y);
v = @(t,x,y) pi*A*cos(pi*h(x,t)).*sin(pi*y).*dhdx(x,t);

% Generate trajectories
dXdT = @(t,x) [u(t,x(1),x(2)); v(t,x(1),x(2))];
X = NaN(length(xr), length(yr), length(tr));
Y = NaN(length(xr), length(yr), length(tr));
% Numerically integrate for path trajectories
for i = 1:length(xr)
    for j = 1:length(yr)
        [t,tmp] = ode45(dXdT, tr, [xr(i), yr(j)]);
        X(i,j,:) = tmp(:,1);
        Y(i,j,:) = tmp(:,2);
    end
end

% Structure for saving data
resol = struct('xlim', xlim, 'ylim', ylim, 'tlim', tlim, 'Nx', n, 'Ny', m, 'Nt', tn);
traj = struct('X', X, 'Y', Y);
params = struct('A', A, 'eps', eps, 'w', w);

% Calculate Cauchy-Green tensors
C = NaN(n,m,2,2);
igrid = NaN(n,m,2);
igrid(:,:,1) = x';
igrid(:,:,2) = y';
for i = 2:(n-1)
    for j = 2:(m-1)
        C(i,j,:,:) = cauchy_green(X, Y, i, j);
    end
end

% Calculate FTLE field
ftle = FTLE(C(2:(end-1), 2:(end-1),:,:),tr(end) - tr(1));

% Calculate vorticity and LAVD
[lavd,w] = LAVD(X,Y,u,v,tr);

% Strain fields
[l1, l2, s] = strain_eigs(C(2:(end-1), 2:(end-1),:,:));

% Collate in struct
vel = struct('u', u, 'v', v, 'h', h, 'dhdx', dhdx);
fields = struct('cauchygreen', C, 'ftle', ftle, 'vort', w, 'lavd', lavd, 'stretch', l1, 'shrink', l2, 'shear', s);
doublegyre = struct('vel', vel, 'traj', traj, 'res', resol, 'params', params, 'fields', fields);

% Save data
save('../data/doublegyre.mat', 'doublegyre');
