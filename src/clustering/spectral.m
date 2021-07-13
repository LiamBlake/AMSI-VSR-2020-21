function [nc,clusters] = spectral(X, Y, T, eps, K, nc, dist, max_cluster, show_plots, varargin)
%SPECTRAL

% Distance metric
if dist == "dyn_dist" || dist == "dyn_dist_gauss"
    dist_handle = @dyn_dist_met;
elseif dist == "kin_dissim"
    dist_handle = @kin_dissim;
end

% Dimensions
n = size(X,1);
m = size(X,2);
tn = length(T);

% Reshape for use with ExhaustiveSearcher
Z = NaN(m*n, 2*tn);
Z(:,1:2:end) = reshape(X, m*n, tn);    % x-coords
Z(:,2:2:end) = reshape(Y, m*n, tn);    % y-coords

% Calculate similarity matrix
W = inf*ones(n*m, n*m);
[idx,dists] = rangesearch(Z, Z, eps, 'Distance', dist_handle);
for i = 1:length(idx)
    %logic = (idx{i} <= i);
    %W(i,idx{i}(logic)) = dists{i}(logic);
    W(i,idx{i}) = dists{i};
end

if dist == "dyn_dist"
    % Invert and set diagonal elements
    W = 1./real(W);
    W(W == inf) = K;
elseif dist == "dyn_dist_gauss"
    W = exp(-W.^2./(2*varargin{1}));
end

% Create graph
%G = graph(W, 'omitselfloops');

% Graph Laplacian & Degree Matrix
%D = diag(degree(G));
D = diag(sum(W,2));
L = D -W;
%n = size(L,1);
%L(1:(n+1):end) = diag(D);
%L = laplacian(G);

% Solve generalised eigenproblem and order
[V,S] = eig(L,D);
S(1,1) = 0;         % Correct for floating point error

%[d,ind] = sort(diag(D));
%D = D(ind,ind);
%V = V(:,ind);

% Determine the optimal number of clusters if not passed
eigengaps = diff(diag(S(1:max_cluster, 1:max_cluster)));
if isnan(nc)
    if show_plots
        figure;
        plot(eigengaps);
    end

    [~,nc] = max(eigengaps);
end
% K-means on U-matrix
U = V(:,1:nc);
clusters = kmeans(U, nc+1);



function D2 = dyn_dist_met(ZI,ZJ)
% Calculation of dynamical distance for

p = size(ZJ,1);
D2 = NaN(p,1);

% Time scaling
tf = 1/(T(end) - T(1));

% Shifted times and ZI for vectorisation
tk = T(1:end-1); tk1 = T(2:end);
ZIk = ZI(1:end-2); ZIk1 = ZI(3:end);

% TODO: vectorise this :/
for j = 1:p
    ZJk = ZJ(j,1:end-2); ZJk1 = ZJ(j,3:end);
    D2(j) = tf*sum(0.5*(tk1 - tk).*( ((ZIk1(1:2:end) - ZJk1(1:2:end)).^2 - (ZIk1(2:2:end) - ZJk1(2:2:end)).^2).^0.5 ...
        + ((ZIk(1:2:end) - ZJk(1:2:end)).^2 - (ZIk(2:2:end) - ZJk(2:2:end)).^2).^0.5 ));
end


end

function D2 = kin_dissim(ZI,ZJ)
% Calculation of dynamical distance for



end



    
end




