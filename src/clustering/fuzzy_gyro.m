% Application of fuzzy clustering with a variety of weighted objective
% functions

global ftle

% Include additional directories
addpath fields fuzzy;

% Load trajectories
load('../data/double_gyro.mat');

n = size(X,1);
m = size(X,2);

% Maximum number of clusters to consider
K = 15;

%% No Weighting, Euclidean distance, K = 2,...,8 clusters
Z = NaN(m*n, 2*length(tr));
Z(:,1:2:end) = reshape(X, m*n, length(tr));    % x-coords
Z(:,2:2:end) = reshape(Y, m*n, length(tr));    % y-coords
[~,U_nwe] = fcm(Z,2);

%% No weighting, Euclidean distance, K = 2 membership probabilities
figure;
h = heatmap(xr,yr,reshape(U_nwe(1,:), n,m)');
h.Colormap = parula;
grid off
h.YDisplayData = flipud(h.YDisplayData);



%% No weighting, dynamic distance, K = 2,...,8 clusters
rng(5);
fc_nw = NaN(K-1,n,m);
U_nw = NaN(K-1,n,m,K);
for k = 2:K
    [fc,~] = fuzzy_lcs(X,Y,tr,k,@dyn_dist_fc,@(x,y,c,u) 1,2,0.1);
    U_nw(k-1,:,:,1:k) = fc;
    
    % Get hard clusters
    [~,fc_nw(k-1,:,:)] = max(fc,[],3);
    
    % Set incoherent cluster
    fc_nw(k, fc_nw(k-1,:,:) < 0.5) = 0;
end

%% No weighting, K = 2 membership probabilities
figure;
h = heatmap(xr,yr,squeeze(U_nw(1,:,:,1))');
h.Colormap = parula;
grid off
h.YDisplayData = flipud(h.YDisplayData);

%% Entropy, no weighting
H_nw = NaN(K-1,n,m);
for k = 2:K
    H_nw(k-1,:,:) = clust_entropy(squeeze(U_nw(k-1,:,:,1:k)));
end
entop_nw = mean(H_nw, [2 3]);

% Plot entropy
figure
plot(2:K, entop_nw)
xlabel('No. Clusters')
ylabel('Mean Normalised Entropy')
title('Unweighted Fuzzy Clustering on Double Gyro')


%% Plot, no weighting
% Set K here
vK = 8;
figure;
gscatter(squeeze(reshape(X(:,:,1),[],1)),squeeze(reshape(Y(:,:,1),[],1)), fc_nw(vK-1,:), 'rgbcmyk', '.', 10);


%% 3-Cluster membership plots
figure;
tiledlayout(3,1);
nexttile; fuzzy_heatmap(xr,yr,squeeze(U_nw(2,:,:,1))', parula);
nexttile; fuzzy_heatmap(xr,yr,squeeze(U_nw(2,:,:,2))', parula);
nexttile; fuzzy_heatmap(xr,yr,squeeze(U_nw(2,:,:,3))', parula);


%% Calculate FTLE field
% Cauchy-Green tensor
[x,y] = meshgrid(xr,yr);
C = NaN(n,m,2,2);
igrid = NaN(n,m,2);
igrid(:,:,1) = x';
igrid(:,:,2) = y';
for i = 2:(n-1)
    for j = 2:(m-1)
        C(i,j,:,:) = cauchy_green(X, Y, igrid, i, j);
    end
end

% FTLE calculation
ftle = FTLE(C(2:(end-1), 2:(end-1),:,:),tr(end) - tr(1));


%% Plot FTLE field
figure;
h = heatmap(xr(2:(end-1)), yr(2:(end-1)), ftle');
h.YDisplayData = flipud(h.YDisplayData);
h.Colormap = parula;
grid off


%% Just FTLE, K = 1,2,...,8 clusters
rng(5);
fc_ftle = NaN(K-1,n-2,m-2);
U_ftle = NaN(K-1,n-2,m-2,K);
opts = [NaN, NaN, NaN, false];
for k = 2:K
    [~,fc] = fcm(ftle(:), k, opts);
    fc = reshape(fc',n-2,m-2,[]);
    U_ftle(k-1,:,:,1:k) = fc;
    
    % Get hard clusters
    [~,fc_ftle(k-1,:,:)] = max(fc,[],3);
    
    % Set incoherent cluster
    fc_ftle(k-1, fc_ftle(k-1,:,:) < 0.5) = 0;
end


%% Entropy, FTLE only
H_ftle = NaN(K-1,n-2,m-2);
for k = 2:K
    H_ftle(k-1,:,:) = clust_entropy(squeeze(U_ftle(k-1,:,:,1:k)));
end
entop_ftle = mean(H_ftle, [2 3]);

% Plot entropy
figure
plot(2:K, entop_ftle)
xlabel('No. Clusters')
ylabel('Mean Normalised Entropy')
title('Fuzzy Clustering of FTLE only on Double Gyro')

%% FTLE only, K = 2 membership probabilities
figure;
h = heatmap(xr(2:(end-1)),yr(2:(end-1)),squeeze(U_ftle(1,:,:,1))');
h.Colormap = parula;
grid off
h.YDisplayData = flipud(h.YDisplayData);





%% Weighted by Variance of FTLE in that cluster, assuming hard clustering
rng(5);

[U_wf,~] = fuzzy_lcs(X(2:(end-1),2:(end-1),:),Y(2:(end-1),2:(end-1),:),tr,2,@dyn_dist_fc, @var_ftle, 2,1e-5);

%% 2-cluster membership plot
figure;
fuzzy_heatmap(xr(2:(end-1)),yr(2:(end-1)),U_wf(:,:,1)', parula);


%% 3-Cluster membership plots
figure;
tiledlayout(3,1);
nexttile; fuzzy_heatmap(xr(2:(end-1)),yr(2:(end-1)),U_wf(:,:,1)', parula);
nexttile; fuzzy_heatmap(xr(2:(end-1)),yr(2:(end-1)),U_wf(:,:,2)', parula);
nexttile; fuzzy_heatmap(xr(2:(end-1)),yr(2:(end-1)),U_wf(:,:,3)', parula);

%% Crude intersection of no weighting and ftle only
U_nw2 = squeeze(U_nw(1,2:(end-1),2:(end-1),1:2));
U_fl2 = squeeze(U_ftle(1,:,:,1:2));

U_crude = NaN(n-2,m-2,4);
U_crude(:,:,1) = U_nw2(:,:,1).*U_fl2(:,:,1);
U_crude(:,:,2) = U_nw2(:,:,1).*U_fl2(:,:,2);
U_crude(:,:,3) = U_nw2(:,:,2).*U_fl2(:,:,1);
U_crude(:,:,4) = U_nw2(:,:,2).*U_fl2(:,:,2);

% Calculate entropy
fprintf("Entropy: %f\n", mean(clust_entropy(U_crude), 'all'))

%% Plot
figure;
tiledlayout(2,2);
nexttile; fuzzy_heatmap(xr(2:(end-1)),yr(2:(end-1)),U_crude(:,:,1)', parula);
nexttile; fuzzy_heatmap(xr(2:(end-1)),yr(2:(end-1)),U_crude(:,:,2)', parula);
nexttile; fuzzy_heatmap(xr(2:(end-1)),yr(2:(end-1)),U_crude(:,:,3)', parula);
nexttile; fuzzy_heatmap(xr(2:(end-1)),yr(2:(end-1)),U_crude(:,:,4)', parula);


%% LAVD Field
lavd = LAVD(X,Y,u,v,tr);


%% Functions

function weight = var_ftle(X,~,C,U)
    global ftle

    K = size(C,1);
    
    % Determine hard clustering
    [~,hard] = max(U,[],3);
    
    % FTLE distributions for each cluster
    weight = NaN(size(X,1), size(X,2), K);
    for k = 1:K
        weight(:,:,k) = var(ftle(hard == k));
    end
    
    
end









% This was all a bit dumb - any weighting which is constant w.r.t. the
% clustering of the trajectories is not going to change the results - the
% objective function is always being multiplied by the same thing.
% % Weighted by FTLE, K = 2,...,8 Clusters
% rng(5);
% U_wf = NaN(K-1,n-2,m-2,K);
% for k = 2:K
%     res_ftle = repmat(ftle,1,1,k);
%     weight_ftle = @(X,Y) res_ftle;
% 
%     [fc,~] = fuzzy_lcs(X(2:(end-1),2:(end-1),:),Y(2:(end-1),2:(end-1),:),tr,k, weight_ftle, 2,1);
%     U_wf(k-1,:,:,1:k) = fc;
% end
% 
% 
% % Weighted by FTLE, K = 2 membership probabilities
% figure
% fuzzy_2heatmap(xr(2:(end-1)),yr(2:(end-1)),squeeze(U_wf(1,:,:,1))', parula);
% 
% 
% % Entropy, weighted by FTLE
% H_wf = NaN(K-1,n-2,m-2);
% for k = 2:K
%     H_wf(k-1,:,:) = clust_entropy(squeeze(U_wf(k-1,:,:,1:k)));
% end
% entop_wf = sum(H_wf, [2 3]);
% 
% Plot entropy
% figure;
% plot(2:K, entop_wf)
% xlabel('No. Clusters')
% ylabel('Mean Normalised Entropy')
% title('Fuzzy Clustering Weighted by FTLE on Double Gyro')
