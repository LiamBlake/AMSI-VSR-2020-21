% Generates trajectories for Bickley Jet flow

global u v

% Constants
U = 62.66;
L = 1770e3;
A = [0.0075 0.15 0.30];
k = pi*(2:2:6);
c = [0.1446*U 0.205*U 0.461*U];
p = c - c(3);

% Functions
usyL = @(y) U*sec(y./L).^2;
u = @(t,x,y) c(1) - usyL(y) - 2*usyL(y).*tanh(y./L).*(A(1)*cos(k(1)*(x - p(1)*t)) + A(2)*cos(k(2)*(x - p(2)*t)) + A(3)*cos(k(3)*(x - p(3)*t)));
v = @(t,x,y) -L*usyL(y).*(A(1)*k(1)*cos(k(1)*(x - p(1)*t)) + A(2)*k(2)*cos(k(2)*(x - p(2)*t)) + A(3)*k(3)*cos(k(3)*(x - p(3)*t)));

% Resolutions
Nx = 5;
Ny = 5;
Nt = 5;
xr = linspace(0,20e6,Nx);
yr = linspace(-3e6,3e6, Ny);
tr = linspace(0,40*24*60*60,Nt);

X = NaN(length(xr), length(yr), length(tr));
Y = NaN(length(xr), length(yr), length(tr));


% Numerically integrate for path trajectories
traj = NaN(Nx,Ny,Nt,2);

% Initial conditions - arrange in matrix
[xp,yp] = meshgrid(xr,yr);
init = [xp(:)'; yp(:)'];
n = size(init,2);

% Integrate trajectories
[~,solved_trajs] = ode45(@(t,x) bickleyjet(t,x,n), tr, init);
traj(:,:,:,1) = reshape(solved_trajs(:,1:2:end)', Nx,Ny,[]);
traj(:,:,:,2) = reshape(solved_trajs(:,2:2:end)', Nx,Ny,[]);


% for i = 1:length(xr)
%     for j = 1:length(yr)
%         [t,tmp] = ode45(dXdT, tr, [xr(i), yr(j)]);
%         X(i,j,:) = tmp(:,1);
%         Y(i,j,:) = tmp(:,2);
%     end
% end

% Save data
%save('../data/bickley_jet.mat', 'X', 'Y', 'xr', 'yr', 'tr');


function dXdt = bickleyjet(t,x,n)
    global u v
    
    % Resized passed values
    x = reshape(x,[],n);
    tp = repmat(t,1,size(x,2));
    
    dxdt = u(t,x(1,:),x(2,:)); 
    dydt = v(t,x(1,:),x(2,:));
    
     % Linearise output
    dXdt = NaN(2*numel(dxdt),1);
    dXdt(1:2:end) = dxdt(:);
    dXdt(2:2:end) = dydt(:);
    
end
