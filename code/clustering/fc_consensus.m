function [mu,v] = fc_consensus(membership,Ki,wi,f,K,m,stop_opt)
%FC_CONSENSUS - Fuzzy consensus clustering
% usage:  [mu,v] = fc_consensus(membership,Ki,f,K,m)
% usage:  [mu,v] = fc_consensus(membership,Ki,f,K,m,stop_opt)
%
% Inputs:
%    membership - cell array of membership matrices for each basic
%                 partition
%    Ki - vector containing the number of clusters (number of columns) in 
%         membership matrix. Must be the same length as membership.
%
%    d - distance function, as function handle with arguments 
%        (l,k,Y,v,inds,wi), where l,k are the indices of the passed 
%        vectors, Y is the multimembership matrix, v is the matrix of 
%        centers, inds is the index vector and wi is the vector of weights.
%
% Outputs:
%    output1 - membership matrix for the consensual clustering
%    output2 - centres 
%
% See also: none
% Author: Liam Blake

% Default arguments
if (nargin < 7)
    thres = 1e-5;
    max_iter = 100;
else
    thres = stop_opt(1);
    max_iter = stop_opt(2);
    
end
% Dimensions
n = size(membership{1}, 1);
r = length(Ki);
sKi = sum(Ki);

% Verify input dimensions


% Construct multimembership data
Y = [];
inds = [1 cumsum(Ki) + 1];
inds(end) = sKi + 1;
for i = 1:r 
    Y = [Y membership{i}];
end

% Initialise mu uniformly
%mu = ones(n,K)*1/K;
% Initialise mu randomly
mu = rand(n,K);
mu = diag(1./sum(mu,2))*mu;

% Initialise piecewise centres v
v = []; update_centers();

% Iterate until stopping criterion is met
in = 0;
obj_old = inf;
while 1
    in = in + 1;
    
    % Update memberships
    update_membership();
    
    % Update centers
    update_centers();
    
    % Calculate objective function
    obj_new = obj();
    
    % Check for convergence
    if (abs(obj_new - obj_old) < thres)
        fprintf("Reached threshold in %i iterations\n", in);
        break;
    end
    obj_old = obj_new;
    
    % Check if reached max iterations
    if (in > max_iter)
        fprintf("Warning: Maximum number of iterations %i exceeded. Solution may not have converged.\n", max_iter)
        break
    end
    
end

% Algorithm functions
    function update_membership()
        % Calculates membership probabilities
        mu = NaN(n,K);
        
        % Distance term
        for k = 1:K
            mu(:,k) = (d(1:n, k).^(-1/(m-1)));
        end
        
        % Normalise rows
        mu = diag(1./sum(mu,2))*mu;
    end

    function update_centers()
        v = NaN(K,sKi);
        for j = 1:r
            call = inds(j):(inds(j+1)-1);
            mukL1 = sum(mu.^m, 1)';
            for k = 1:K
                v(k,call) = sum(mu(:,k).^m.*Y(:,call),1);
            end
            v(:,call) = v(:,call)./mukL1;
        end
    end

    function dist = d(l,k)
        dist = zeros(size(l'));
        for j = 1:r
            call = inds(j):(inds(j+1)-1);
            dist = dist + wi(j)*f(Y(l,call), v(k,call).*ones(length(l),length(call)));
        end
    end

    function Jfm = obj()
        Jfm = 0;
        for k = 1:K
            Jfm = Jfm + sum(mu(:,k).^m.*d(1:n, k),'all');
        end
    end

end