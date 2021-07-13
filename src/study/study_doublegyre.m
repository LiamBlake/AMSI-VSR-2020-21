% Explores different combinations of clusterings on the double gyre flow
% and the corresponding LCS fields.
tic;
addpath clustering lavd_files visualise
%set(groot, 'defaultAxesTickLabelInterpreter','LaTeX');
%set(groot, 'defaultLegendInterpreter','LaTeX');

conFigure(11, 4/3);

% Global variables
SHOW_PLOTS = 1;

% Load data and parameters
load('../data/doublegyre.mat')

X = doublegyre.traj.X;
Y = doublegyre.traj.Y;
T = linspace(doublegyre.res.tlim(1),doublegyre.res.tlim(2),doublegyre.res.Nt);

[Nx,Ny] = size(X, [1 2]);
Nt = doublegyre.res.Nt;

ftle = doublegyre.fields.ftle;
lavd = doublegyre.fields.lavd;

% Maximum number of clusters to consider
MAX_K = 15;

% Scale LAVD to time units
lavd = 1/(T(end) - T(1))*lavd;

% Normalise FTLE and LAVD
ftle = ftle./max(ftle,[],'all');
lavd = lavd./max(lavd,[],'all');

xp = squeeze(X(:,1,1));
yp = squeeze(Y(1,:,1));
xp2 = xp(2:end-1);
yp2 = yp(2:end-1);


% Reshaped trajectories
Z = NaN(Nx*Ny,2*Nt);
Z(:,1:2:end) = reshape(X,Nx*Ny,Nt);
Z(:,2:2:end) = reshape(Y,Nx*Ny,Nt);

lavd2 = lavd(2:end-1,2:end-1);
Z2 = NaN((Nx-2)*(Ny-2),2*Nt);
Z2(:,1:2:end) = reshape(X(2:end-1,2:end-1,:),(Nx-2)*(Ny-2),Nt);
Z2(:,2:2:end) = reshape(Y(2:end-1,2:end-1,:),(Nx-2)*(Ny-2),Nt);
% Scale x-axis to ensure all values of Z are in [0,1]
Z2(:,1:2:end) = 0.5*Z2(:,1:2:end);

%% Plot FTLE Field
if SHOW_PLOTS
    figure;
    cbar = colourplot(X(2:end-1,1,1),Y(1,2:end-1,1),ftle, "Double Gyre FTLE Field", "$x$", "$y$", "FTLE", winter);
    cbar.Limits = [0 1];
end

%% Plot LAVD Field
if SHOW_PLOTS
    figure;
    cbar = colourplot(X(:,1,1),Y(1,:,1),lavd, "Double Gyre LAVD Field", "$x$", "$y$", "(Scaled) LAVD", winter);
    cbar.Limits = [0 1];
end

%% Clustering applied to each field
% Set seed for reproducible results
rng(945210121)

% Membership probabilities
fcm_traj = NaN(2,MAX_K-1,Nx,Ny,MAX_K);
fcm_ftle = NaN(2,MAX_K-1,Nx-2,Ny-2,MAX_K);
fcm_lavd = NaN(2,MAX_K-1,Nx,Ny,MAX_K);
fcm_ftlelavd_concat = NaN(2,MAX_K-1,Nx-2,Ny-2,MAX_K);
fcm_all_concat = NaN(2,MAX_K-1,Nx-2,Ny-2,MAX_K);

% Trajectory Centers
cent_traj = NaN(2,MAX_K-1,MAX_K,2*Nt);


for k = 2:MAX_K
    
    for m = 1:2
        fprintf("Clustering fields with m = %i and %i clusters...\n", m, k)
        
        % Trajectories
        [fc,cent_traj(m,k-1,1:k,:)] = fcm_wrapper(Z,k,m);
        fcm_traj(m,k-1,:,:,1:k) = reshape(fc,Nx,Ny,[]);
        
        % Fuzzy c-means on FTLE (Euclidean distance)
        [fc,~] = fcm_wrapper(ftle(:),k,m);
        fcm_ftle(m,k-1,:,:,1:k) = reshape(fc,Nx-2,Ny-2,[]);
        
        % Fuzzy c-means on LAVD (Euclidean distance)
        [fc,~] = fcm_wrapper(lavd(:), k, m);
        fcm_lavd(m,k-1,:,:,1:k) = reshape(fc,Nx,Ny,[]);
        
        % Fuzzy c-means on FTLE and LAVD
        [fc,~] = fcm_wrapper([ftle(:) lavd2(:)], k, m);
        fcm_ftlelavd_concat(m,k-1,:,:,1:k) = reshape(fc,Nx-2,Ny-2,[]);
        
        
        % Let's just shove all the data together
        [fc,~] = fcm_wrapper([ftle(:) lavd2(:) Z2], k, m);
        fcm_all_concat(m,k-1,:,:,1:k) = reshape(fc,Nx-2,Ny-2,[]);
        
    end
    break;
end


%% Entropy plots
if SHOW_PLOTS
    figure;
    t = tiledlayout(2,2);
    nexttile; [ent_traj, ~] = plot_metrics(fcm_traj(2,:,:,:,:), "Trajectories", 2:MAX_K);
    nexttile; [ent_ftle, ~] = plot_metrics(fcm_ftle(2,:,:,:,:), "FTLE", 2:MAX_K);
    nexttile; [ent_lavd, ~] = plot_metrics(fcm_lavd(2,:,:,:,:), "LAVD", 2:MAX_K);
    nexttile; [ent_ftlelavd_concat,~] = plot_metrics(fcm_ftlelavd_concat(2,:,:,:,:), "FTLE \& LAVD", 2:MAX_K);

    figure;
    plot_metrics(fcm_all_concat(2,:,:,:,:), "All three", 2:MAX_K);
end

%% Optimal no. Clusters
N_traj = 2;
N_ftle = 6;
N_lavd_clust = 6;
N_ftlelavd_clust = 2;
N_all_clust = 4;        % seems to be overpowered by trajectory clusters - makes sense, perhaps should do some scaling

%% Membership Plots - optimal no. clusters
if SHOW_PLOTS
    % Trajectories
    plot_memberships(fcm_traj(2,N_traj-1,:,:,1:N_traj), "Trajectory Membership Values",xp,yp);
    
    
    % FTLE
    plot_memberships(fcm_ftle(2,N_ftle-1,:,:,1:N_ftle), "FTLE Membership Values",xp2,yp2);
    
    
    % LAVD
    plot_memberships(fcm_lavd(2,N_lavd_clust-1,:,:,1:N_lavd_clust), "LAVD (Clustering) Membership Values", xp, yp);
    
    % FTLE AND LAVD
    plot_memberships(fcm_ftlelavd_concat(2,N_ftlelavd_clust-1,:,:,1:N_ftlelavd_clust), "FTLE \& LAVD Membership Values",xp2,yp2);
   
    % All three
    plot_memberships(fcm_all_concat(2,N_all_clust-1,:,:,1:N_all_clust), "FTLE \& LAVD Membership Values",xp2,yp2);

end


%% Trajectory Centres
if SHOW_PLOTS
    plot_trajectory_centres(cent_traj(N_traj-1,1:N_traj,:))
    
end



%% LAVD Vortex Extraction
% Using code from https://github.com/Hadjighasem/Lagrangian-Averaged-Vorticity-Deviation-LAVD
% G. Haller, A. Hadjighasem, M. Farazamand & F. Huhn, Defining coherent vortices objectively from the vorticity. submitted (2015): arXiv:1506.04061 [physics.flu-dyn]
Nct = 50;
min_length = 1;
def_thres = 0.5;
xi = squeeze(X(:,:,1))'; yi =  squeeze(Y(:,:,1))';
lavd_bnd = ContourExtraction(lavd',xi,yi,Nct,min_length,def_thres);

% Convert to membership values
[M_lavd, N_lavd] = contour_membership(xi,yi,lavd_bnd);
M_lavd = permute(M_lavd, [2 1 3]);
if SHOW_PLOTS
    plot_memberships_crisp(M_lavd, "Double Gyre LAVD Contours", xp,yp,"$x$", "$y$");
end


if SHOW_PLOTS
    figure;
    cbar = colourplot(X(:,1,1),Y(1,:,1),lavd, "Double Gyre LAVD Field", "$x$", "$y$", "(Scaled) LAVD", winter);
    for kk=1:numel(lavd_bnd.xc); hold on; plot(lavd_bnd.xc{kk},lavd_bnd.yc{kk},'r','linewidth',3); end
    plot(lavd_bnd.xp,lavd_bnd.yp,'or','MarkerFaceColor','r','MarkerSize',4);
    cbar.Limits = [0 1];
end


%% Combinations of cluster configurations using FCC


%% Trajectories & FTLE
% Combine trajectories and FTLE clusters with 2 to 15 clusters
membership = {reshape(fcm_traj(2,N_traj-1,2:(end-1),2:(end-1),1:N_traj), (Nx-2)*(Ny-2), []), reshape(fcm_ftle(2,N_ftle-1,:,:,1:N_ftle), (Nx-2)*(Ny-2), [])};
Ki = [N_traj N_ftle];
wts = [0.5 0.5];

cons_trajftle = NaN(2,MAX_K-1,Nx-2,Ny-2,MAX_K);

for k = 2:MAX_K
    
    for m = 1:2
        [fc,~] = fcc_wrapper(membership,Ki,wts,"Euclidean",k,2);
        cons_trajftle(m,k-1,:,:,1:k) = reshape(fc, Nx-2, Ny-2, []);
    end
end

%% Metric plots
if SHOW_PLOTS
    % Entropy
    figure;
    plot_metrics(cons_trajftle(2,:,:,:,:), "Consensus of Trajectories \& FTLE", 2:MAX_K, "");
    
end

N_trajftle = 6;


%% Membership plots
if SHOW_PLOTS
    plot_memberships(cons_trajftle(2,N_trajftle-1,:,:,1:N_trajftle), "LAVD partitions", xp,yp);
end

%% Trajectories & LAVD

% Combine trajectories and FTLE clusters with 2 to 15 clusters
membership = {reshape(fcm_traj(2,N_traj-1,:,:,1:N_traj), Nx*Ny, []), reshape(M_lavd, Nx*Ny, [])};
Ki = [N_traj N_lavd];
wts = [0.5 0.5];

cons_trajlavd = NaN(2,MAX_K-1,Nx,Ny,MAX_K);

for k = 2:MAX_K
    % Soft clusters, m = 2
    [fc,~] = fcc_wrapper(membership,Ki,wts,"Euclidean",k,2);
    cons_trajlavd(2,k-1,:,:,1:k) = reshape(fc, Nx, Ny, []);
   
end


%% Metrics plot
if SHOW_PLOTS
    figure;
    plot_metrics(cons_trajlavd(2,:,:,:,:), "Consensus of Trajectory Clusters and LAVD Contours", 2:MAX_K);
    
end

N_trajlavd = 4;

%% Membership values
if SHOW_PLOTS
    plot_memberships(cons_trajlavd(2,N_trajlavd-1,:,:,1:N_trajlavd), "Consensus of Trajectory Clusters and LAVD Contours", xp,yp);
end


%% FTLE & LAVD
% Combine LAVD and FTLE clusters with 2 to 15 clusters
membership = {reshape(fcm_ftle(2,N_ftle-1,:,:,1:N_ftle), (Nx-2)*(Ny-2), []), reshape(M_lavd(2:(end-1),2:(end-1),:), (Nx-2)*(Ny-2), [])};
Ki = [N_ftle N_lavd];
wts = [0.5 0.5];

cons_ftlelavd = NaN(MAX_K-1,Nx-2,Ny-2,MAX_K);

for k = 2:MAX_K
    
    % Soft clusters, m = 2
    [fc,~] = fcc_wrapper(membership,Ki,wts,"Euclidean",k,2);
    cons_ftlelavd(k-1,:,:,1:k) = reshape(fc, Nx-2, Ny-2, []);
    
end

%% Metrics plot
if SHOW_PLOTS
    figure;
    plot_metrics(cons_ftlelavd(:,:,:,:), "Consensus of FTLE Clusters and LAVD Contours", 2:MAX_K);
    
end

N_ftlelavd = 5;

%% Membership values
if SHOW_PLOTS
    plot_memberships(cons_ftlelavd(N_ftlelavd-1,:,:,1:N_ftlelavd), "Consensus of FTLE Clusters and LAVD Contours", xp2,yp2);
end


%% All three

wts = [1/3 1/3 1/3];
Ki = [N_traj N_ftle N_traj];
membership = {reshape(fcm_traj(2,N_traj-1,2:(end-1),2:(end-1),1:N_traj), (Nx-2)*(Ny-2), []), ...
    reshape(fcm_ftle(2,N_ftle-1,:,:,1:N_ftle), (Nx-2)*(Ny-2), []), ...
    reshape(M_lavd(2:(end-1),2:(end-1),:), (Nx-2)*(Ny-2), [])};

% Variable number of clusters
cons_all = NaN(MAX_K-1,Nx-2,Ny-2,MAX_K);
for k = 2:MAX_K
    [fc,~] = fcc_wrapper(membership,Ki,wts,"Euclidean",k,2);
    cons_all(k-1,:,:,1:k) = reshape(fc, Nx-2, Ny-2, []);
    
end

%% Metrics plot
if SHOW_PLOTS
    figure;
    plot_metrics(cons_all(:,:,:,:), "Consensus of all three fields", 2:MAX_K);
    
end

N_all = 8;

%% Membership values
if SHOW_PLOTS
    plot_memberships(cons_all(N_all-1,:,:,1:N_all), "Consensus of all three fields", xp2,yp2);
end


%% FTLE Ridge Extraction
% Rudimentary method - get points with FTLE value which is 80% of maximum
ftleridge = ftle > 0.75;
ftleridge(:,:,2) = ~ftleridge;
%ftleridge(:,:,1) = ~ftleridge(:,:,2);

ftleridge = double(ftleridge);

if SHOW_PLOTS
plot_memberships_crisp(ftleridge, "Double Gyre FTLE Ridges", xp2,yp2, "$x$", "$y$");
end


%% FTLE & LAVD with convex weighting
N_ftlelavd = 5;
membership = {reshape(fcm_ftle(2,N_ftle-1,:,:,1:N_ftle), (Nx-2)*(Ny-2), []), reshape(M_lavd(2:(end-1),2:(end-1),:), (Nx-2)*(Ny-2), [])};
Ki = [N_ftle N_lavd];

Na = 20;
consalph_ftlelavd = NaN(Na,Nx-2,Ny-2,N_ftlelavd);
alphas = linspace(0,1,Na);

for i = 1:Na
    a = alphas(i);
    wts = [a 1 - a];
    
    [fc,~] = fcc_wrapper(membership,Ki,wts,"Euclidean",N_ftlelavd,2);
    consalph_ftlelavd(i,:,:,:) = reshape(fc,Nx-2,Ny-2,[]);
end


%% Plot results
if SHOW_PLOTS
    figure;
    plot_metrics(consalph_ftlelavd, "", alphas);
    xlabel("$\alpha$", "Interpreter", "LaTeX")
end
%% Memberships
if SHOW_PLOTS
    for i = [3 5 10 15 18]
        plot_memberships(consalph_ftlelavd(i,:,:,:), "$\alpha = " + alphas(i) + "$",xp2,yp2, "", "")
    end
end

%% Save results
save('../data/doublegyre_indivclusters.mat', 'fcm_ftle', 'fcm_traj', 'fcm_lavd', 'fcm_ftlelavd_concat')
save('../data/doublegyre_consclusters.mat', 'cons_all', 'cons_ftlelavd', 'cons_trajftle', 'cons_trajlavd')
save('../data/doublegyre_otherclusters.mat', 'consalph_ftlelavd');

toc;



%% Combining FTLE ridge, clustering and LAVD regions into one
traj_clust = squeeze(fcm_traj(2,1,2:end-1,2:end-1,1:2));
sensible = NaN(MAX_K-1,Nx-2,Ny-2,MAX_K);

wts = [1/3 1/3 1/3];
Ki = [2 2 N_lavd];
membership = {reshape(traj_clust, (Nx-2)*(Ny-2), 2), ...
    reshape(ftleridge, (Nx-2)*(Ny-2), 2), ...
    reshape(M_lavd(2:end-1,2:end-1,:), (Nx-2)*(Ny-2), N_lavd)};

for k = 2:MAX_K
    [fc,~] = fcc_wrapper(membership,Ki,wts,"Euclidean",k,2);
    sensible(k-1,:,:,1:k) = reshape(fc, Nx-2, Ny-2, []);
end

%% Plot metrics

if SHOW_PLOTS
    figure;
    plot_metrics(sensible, "", 2:MAX_K);
end

N_sensible = 5;


%% Plot memberships
plot_memberships(sensible(N_sensible-1,:,:,1:N_sensible), "Atlantic Consensus Membership Values", xp2,yp2, "Longitude", "Latitude");
