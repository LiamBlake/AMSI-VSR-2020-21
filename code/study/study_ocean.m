% Performs clustering analysis on Atlantic ocean data, using the saved
% output from traj_ocean.m

% Required directories
addpath clustering fields lavd_files visualise

% Figure setup


% Load data
load('../data/atlantic_res.mat')
load('../data/atlantic_traj_X.mat')
load('../data/atlantic_traj_Y.mat')
load('../data/atlantic_fields.mat')
load('../data/atlantic_temp.mat')

ftle = fields.ftle;
lavd = fields.lavd;
temp = temp_ex.temp;
time = ocean.res.tlim;

Nx = ocean.res.Nx;
Ny = ocean.res.Ny;
Nt = ocean.res.Nt;

xp = squeeze(X(:,1,1));
yp = squeeze(Y(1,:,1));
xp2 = xp(2:end-1);
yp2 = yp(2:end-1);

% Scale LAVD to time units
lavd = 1/(time(end) - time(1))*lavd;

% Normalise FTLE and LAVD
ftle = ftle./max(ftle,[],'all');
nlavd = lavd./max(lavd,[],'all');

%% NAN the LANd
% Reshaped trajectories
Z = NaN(Nx*Ny,2*Nt);
Z(:,1:2:end) = reshape(X,Nx*Ny,Nt);
Z(:,2:2:end) = reshape(Y,Nx*Ny,Nt);

% Remove trajectories which remain constant
isconstx = sum(abs(diff(Z(:,1:2:end),1,2)), 2) == 0;
isconsty = sum(abs(diff(Z(:,2:2:end),1,2)), 2) == 0;
isconst = isconstx & isconsty;
Z(isconst,:) = [];

% Interpolate NaNs in temperature field - regions corresponding to land are
% removed based on trajectory analysis above
for t = 1:Nt
    temp(:,:,t) = inpaint_nans(temp(:,:,t), 3);
end

% Reshaped temperature
temp = reshape(temp, Nx*Ny,Nt);
temp(isconst,:) = NaN;
temp = reshape(temp, Nx,Ny,Nt);

T = reshape(temp, Nx*Ny,Nt);
T(isconst,:) = [];

% Reshape FTLE
isconst2 = reshape(isconst, Nx, Ny);
isconst2 = isconst2(2:end-1, 2:end-1);
isconst2 = reshape(isconst2, (Nx-2)*(Ny-2),1);

ftle_l = ftle(:);
ftle_l(isconst2) = [];

lavd_l = nlavd(:);
lavd_l(isconst) = [];
lavd2 = nlavd(2:end-1, 2:end-1);
lavd2 = lavd2(:);
lavd2(isconst2) = [];

% Normalise temperature
temp = temp./max(temp,[],'all');

%% Plot FTLE Field
conFigure(8.5);
figure;
cbar = colourplot_land(X(2:end-1,1,1),Y(1,2:end-1,1),ftle, "Atlantic FTLE Field", "Longitude", "Latitude", "FTLE", isconst2, winter);
cbar.Limits = [0 1];

%% Rudimentary FTLE ridge extraction
% Extract ridge from FTLE by taking values which are at least 80% of the
% maximum
ftleridge = ftle > 0.8;
ftleridge(:,:,2) = ~ftleridge;
%ftleridge(:,:,2) = ftle < 0.2;
%ftleridge(:,:,3) = ~(ftleridge(:,:,1) | ftleridge(:,:,2));
%ftleridge = double(ftleridge);
%ftleridge(repmat(reshape(isconst2,Nx-2,Ny-2), 1,1,3)) = NaN;

plot_memberships_crisp(ftleridge, "Atlantic FTLE Ridges", xp2,yp2, "Longitude", "Latitude");


%% Better FTLE Ridge extraction
[fridge,iridge,lridge] = tfridge(ftle',yp2);
rvals = ftle(lridge);

figure;
cbar = colourplot_land(X(2:end-1,1,1),Y(1,2:end-1,1),ftle, "Atlantic FTLE Field", "Longitude", "Latitude", "FTLE", isconst2, winter);
cbar.Limits = [0 1];
hold on
plot(xp2, fridge, 'ko')%, 'LineWidth', 4)

%% Plot LAVD Field
figure;
colourplot_land(X(:,1,1),Y(1,:,1),lavd, "Atlantic LAVD Field", "Longitude", "Latitude", "LAVD", isconst, parula);


%% Plot Temperature Field
figure;
colourplot_land(xp,yp,temp(:,:,1), "Atlantic Temperature on 25/06/2018", "Longitude", "Latitude", "(Normalised) Temperature", isconst, jet);

figure; 
colourplot_land(xp,yp,temp(:,:,end), "Atlantic Temperature on 25/07/2018", "Longitude", "Latitude", "Temperature (Celsius)", isconst, jet);

%% Global Variables
MAX_K = 15;



%% Individual Clusters
% Set seed for reproducible results
rng(1109200121)

% Membership probabilities
fcm_traj = NaN(MAX_K-1,Nx,Ny,MAX_K);
fcm_ftle = NaN(MAX_K-1,Nx-2,Ny-2,MAX_K);
fcm_lavd = NaN(MAX_K-1,Nx,Ny,MAX_K);
fcm_temp = NaN(MAX_K-1,Nx,Ny,MAX_K);
fcm_ftlelavd = NaN(MAX_K-1,Nx-2,Ny-2,MAX_K);

% Cluster centers
cen_traj = NaN(MAX_K-1,MAX_K,2*Nt);
cen_ftle = NaN(MAX_K-1,MAX_K,1);
cen_lavd = NaN(MAX_K-1,MAX_K,1);
cen_temp = NaN(MAX_K-1,MAX_K,Nt);
cen_ftlelavd = NaN(MAX_K-1,MAX_K,2);

% Entropy of each object
ent_traj = NaN(MAX_K-1,Nx,Ny);
ent_ftle = NaN(MAX_K-1,Nx-2,Ny-2);
ent_temp = NaN(MAX_K-1,Nx,Ny);

for k = 2:MAX_K
    % Trajectories
    [fc,cen_traj(k-1,1:k,:)] = fcm_wrapper(Z,k,2);
    fc_full = NaN(Nx*Ny,k);
    fc_full(~isconst,:) = fc;
    fcm_traj(k-1,:,:,1:k) = reshape(fc_full,Nx,Ny,[]);
    
    % Fuzzy c-means on FTLE (Euclidean distance)
    [fc,cen_ftle(k-1,1:k,:)] = fcm_wrapper(ftle_l, k, 2);
    fc_full = NaN((Nx-2)*(Ny-2),k);
    fc_full(~isconst2,:) = fc;
    fcm_ftle(k-1,:,:,1:k) = reshape(fc_full,Nx-2,Ny-2,[]);
    
    % Fuzzy c-means on temperature (Euclidean distance)
    [fc,cen_temp(k-1,1:k,:)] = fcm_wrapper(T, k, 2);
    fc_full = NaN(Nx*Ny,k);
    fc_full(~isconst,:) = fc;
    fcm_temp(k-1,:,:,1:k) = reshape(fc_full,Nx,Ny,[]);
    
    % Fcm on LAVD field
    [fc,cen_lavd(k-1,1:k,:)] = fcm_wrapper(lavd_l, k, 2);
    fc_full = NaN(Nx*Ny,k);
    fc_full(~isconst,:) = fc;
    fcm_lavd(k-1,:,:,1:k) = reshape(fc_full,Nx,Ny,[]);
    
    % FCM on FTLE and LAVD field
    [fc,cen_ftlelavd(k-1,1:k,:)] = fcm_wrapper([ftle_l lavd2], k, 2);
    fc_full = NaN((Nx-2)*(Ny-2),k);
    fc_full(~isconst2,:) = fc;
    fcm_ftlelavd(k-1,:,:,1:k) = reshape(fc_full,Nx-2,Ny-2,[]);
    
end


%% Metric plots
conFigure(8.5);
% Trajectories
figure;
plot_metrics(fcm_traj, "Trajectory FCM", 2:MAX_K);
xlim([2,MAX_K])
xticks(2:MAX_K)
% FTLE


%nexttile; plot_metrics(fcm_ftle, "FTLE", 2:MAX_K);

% Temperature
figure;
plot_metrics(fcm_temp, "Temperature FCM", 2:MAX_K);
xlim([2,MAX_K])
% LAVD 
%nexttile; plot_metrics(fcm_lavd, "LAVD (Clustering)", 2:MAX_K);

% LAVD & Temperature
%nexttile; plot_metrics(fcm_ftlelavd, "FTLE and LAVD", 2:MAX_K);


%% "Optimal" No. Clusters
N_traj = 4;             % Clear minimum
N_ftle = 4;             % Elbow criterion
N_temp = 3;             % Heuristic, need to think about this one
N_lavd_clust = 6;       % Elbow
N_ftlelavd_clust = 2;   % Elbow


%% Visualise cluster memberships
conFigure(11.33)

% Membership plots
plot_memberships(fcm_traj(N_traj-1,:,:,1:N_traj), "Trajectory FCM", xp, yp, "Longitude", "Latitude")

plot_memberships(fcm_ftle(N_ftle-1,:,:,1:N_ftle), "FTLE Cluster Membership", xp2, yp2, "Longitude", "Latitude")

plot_memberships(fcm_temp(N_temp-1,:,:,1:N_temp), "Temperature FCM", xp, yp, "Longitude", "Latitude")

plot_memberships(fcm_lavd(N_lavd_clust-1,:,:,1:N_lavd_clust), "LAVD Cluster Membership", xp, yp, "Longitude", "Latitude")

plot_memberships(fcm_ftlelavd(N_ftlelavd_clust-1,:,:,1:N_ftlelavd_clust), "FTLE and LAVD Cluster Membership", xp2, yp2, "Longitude", "Latitude")


%% Plot trajectory cluster centres 
plot_trajectory_centres(cen_traj(N_traj-1,1:N_traj,:))
xlim([xp(1), xp(end)]);
ylim([yp(1), yp(end)]);


%% Extract LAVD contours
% Using code from https://github.com/Hadjighasem/Lagrangian-Averaged-Vorticity-Deviation-LAVD
% G. Haller, A. Hadjighasem, M. Farazamand & F. Huhn, Defining coherent vortices objectively from the vorticity. submitted (2015): arXiv:1506.04061 [physics.flu-dyn]
Nct = 100;
min_length = 1.5;
def_thres = 1;
xi = squeeze(X(:,:,1))'; yi = squeeze(Y(:,:,1))';
lavd_bnd = ContourExtraction(lavd',xi,yi,Nct,min_length,def_thres);


% Convert to membership values
[M_lavd, N_lavd] = contour_membership(xi,yi,lavd_bnd);
M_lavd = permute(M_lavd, [2 1 3]);

% Deal with Land
M_lavd = reshape(M_lavd, Nx*Ny, []);
M_lavd(isconst) = NaN;
M_lavd = reshape(M_lavd, Nx, Ny, []);
plot_memberships_crisp(M_lavd, "Atlantic LAVD Partitions", xp,yp, "Longitude", "Latitude");



%% Consensus Clustering





%% Trajectories and Temperature
membership = {reshape(fcm_traj(N_traj-1,:,:,1:N_traj), Nx*Ny, []), reshape(fcm_temp(N_temp-1,:,:,1:N_temp), Nx*Ny, [])};
membership{1}(isconst,:) = []; membership{2}(isconst,:) = []; 
Ki = [N_traj N_temp];
wts = [0.5 0.5];

cons_trajtemp = NaN(MAX_K-1,Nx,Ny,MAX_K);

for k = 2:MAX_K
    [fc,~] = fcc_wrapper(membership,Ki,wts,"Euclidean",k,2);
    fc_full = NaN(Nx*Ny,k);
    fc_full(~isconst,:) = fc;
    cons_trajtemp(k-1,:,:,1:k) = reshape(fc_full,Nx,Ny,[]);
end


%% Visualise results
conFigure(8.5);
figure;
plot_metrics(cons_trajtemp, "Trajectory \& Temperature FCC", 2:MAX_K);
xlim([2 MAX_K])
N_trajtemp = 8;

%%
conFigure(18/1.3, 1.8)
plot_memberships(cons_trajtemp(N_trajtemp-1,:,:,1:N_trajtemp), "Trajectory \& Temperature FCC", xp, yp, "Longitude", "Latitude")




%% Temperature & LAVD contours 
membership = {reshape(M_lavd, Nx*Ny, []), reshape(fcm_temp(N_temp-1,:,:,1:N_temp), Nx*Ny, [])};
membership{1}(isconst,:) = []; membership{2}(isconst,:) = []; 
Ki = [N_lavd N_temp];
wts = [0.5 0.5];
N_lavdtemp = 6;%N_lavd + N_temp;

[fc,~] = fcc_wrapper(membership,Ki,wts,"Euclidean",N_lavdtemp,2);
fc_full = NaN(Nx*Ny,N_lavdtemp);
fc_full(~isconst,:) = fc;
tmp = reshape(fc_full,Nx,Ny,[]);


%% Plot results
plot_memberships(tmp, "Temperature & LAVD Consensus Cluster Membership", xp, yp, "Longitude", "Latitude")






%% Trajectory & Temperature with convex weightings
Na = 20;
consalph_trajtemp = NaN(Na,Nx,Ny,N_trajtemp);
alphas = linspace(0,1,Na);

for i = 1:Na
    a = alphas(i);
    wts = [a 1 - a];
    
    [fc,~] = fc_consensus(membership,Ki,wts,f,N_trajtemp,2);
    fc_full = NaN(Nx*Ny,N_trajtemp);
    fc_full(~isconst,:) = fc;
    consalph_trajtemp(i,:,:,:) = reshape(fc_full,Nx,Ny,[]);
end


%% Plot results
plot_metrics(consalph_trajtemp, "", alphas);
xlabel("$\alpha$", "Interpreter", "LaTeX")

plot_memberships(consalph_trajtemp(floor(3*Na/4),:,:,:), "Trajectory Cluster Membership for $\alpha = " + alphas(floor(3*Na/4)) + "$", xp, yp, "Longitude", "Latitude")

%% Combining FTLE ridges and LAVD regions

fcc_ftlelavd = NaN(MAX_K-1,Nx-2,Ny-2,MAX_K);

wts = [1/2 1/2];
Ki = [2 N_lavd];
membership = {reshape(ftleridge, (Nx-2)*(Ny-2), 2), ...
    reshape(M_lavd(2:end-1,2:end-1,:), (Nx-2)*(Ny-2), N_lavd)};

for i = 1:2
    membership{i}(isconst2,:) = []; 
end

for k = 2:MAX_K
    [fc,~] = fcc_wrapper(membership,Ki,wts,"Euclidean",k,2);
    fc_full = NaN((Nx-2)*(Ny-2),k);
    fc_full(~isconst2,:) = fc;
    fcc_ftlelavd(k-1,:,:,1:k) = reshape(fc_full, Nx-2, Ny-2, []);
end



%% Plot metrics

if 1
    conFigure(8.5);
    figure;
    
    plot_metrics(fcc_ftlelavd(1:MAX_K-1,:,:,:), "FTLE \& LAVD FCC", 2:MAX_K);
    xlim([2 15])
    set(gca, 'YTickLabel', get(gca, 'YTick'));
end

N_ftlelavd = 8;


%% Plot memberships
conFigure(18/1.3, 1.8)
plot_memberships(fcc_ftlelavd(N_ftlelavd-1,:,:,1:N_ftlelavd), "FTLE \& LAVD FCC", xp2,yp2);




%% Combining FTLE ridge, temperature clustering and LAVD regions into one
temp_clust = squeeze(fcm_temp(N_temp-1,2:end-1,2:end-1,1:N_temp));
sensible = NaN(MAX_K-1,Nx-2,Ny-2,MAX_K);

wts = [1/3 1/3 1/3];
Ki = [N_temp 2 N_lavd];
membership = {reshape(temp_clust, (Nx-2)*(Ny-2), N_temp), ...
    reshape(ftleridge, (Nx-2)*(Ny-2), 2), ...
    reshape(M_lavd(2:end-1,2:end-1,:), (Nx-2)*(Ny-2), N_lavd)};

for i = 1:3
    membership{i}(isconst2,:) = []; 
end

for k = 2:MAX_K
    [fc,~] = fcc_wrapper(membership,Ki,wts,"Euclidean",k,2);
    fc_full = NaN((Nx-2)*(Ny-2),k);
    fc_full(~isconst2,:) = fc;
    sensible(k-1,:,:,1:k) = reshape(fc_full, Nx-2, Ny-2, []);
end



%% Plot metrics

if SHOW_PLOTS
    figure;
    plot_metrics(sensible, "", 2:MAX_K);
end

N_sensible = 5;


%% Plot memberships
plot_memberships(sensible(N_sensible-1,:,:,1:N_sensible), "Consensus of all three fields", xp2,yp2);


%% Combining all four
MAX_K = 20;
temp_clust = squeeze(fcm_temp(N_temp-1,2:end-1,2:end-1,1:N_temp));
traj_clust = squeeze(fcm_traj(N_traj-1,2:end-1,2:end-1,1:N_traj));
all4 = NaN(MAX_K-1,Nx-2,Ny-2,MAX_K);

wts = [1/4 1/4 1/4 1/4];
Ki = [N_traj N_temp 2 N_lavd];
membership = {reshape(traj_clust, (Nx-2)*(Ny-2), N_traj), ...
    reshape(temp_clust, (Nx-2)*(Ny-2), N_temp), ...
    reshape(ftleridge, (Nx-2)*(Ny-2), 2), ...
    reshape(M_lavd(2:end-1,2:end-1,:), (Nx-2)*(Ny-2), N_lavd)};

for i = 1:4
    membership{i}(isconst2,:) = []; 
end

for k = 2:MAX_K
    [fc,~] = fcc_wrapper(membership,Ki,wts,"Euclidean",k,2);
    fc_full = NaN((Nx-2)*(Ny-2),k);
    fc_full(~isconst2,:) = fc;
    all4(k-1,:,:,1:k) = reshape(fc_full, Nx-2, Ny-2, []);
end



%% Plot metrics

if 1
    figure;
    conFigure(8.5);
    plot_metrics(all4(1:MAX_K-1,:,:,:), "All Four Basic Partitions", 2:MAX_K);
end

N_all4 = 17;


%% Plot memberships
conFigure(18/1.3, 1.8)
plot_memberships(all4(N_all4-1,:,:,1:N_all4), "All Four Basic Partitons FCC", xp2,yp2);





%% Save Results
save('../data/ocean_indivclusters.mat', 'fcm_traj', 'fcm_temp', 'fcm_ftle', 'fcm_lavd', 'M_lavd', 'fcm_ftlelavd');




