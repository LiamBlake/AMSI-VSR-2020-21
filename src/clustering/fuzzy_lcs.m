function [U,C] = fuzzy_lcs(X,Y,T,Nc,dist_func, weight_func,m,thres,varargin)
%FUZZY_LCS Applies fuzzy

% Failsafe in case of lack of convergence
MAX_ITER = 100;

% Randomly initialise centers from initial points
C = NaN(Nc, 2*size(X,3));
x_idx = randi([1 size(X,1)],1,Nc);
y_idx = randi([1 size(X,2)],1,1);
C(:,1:2:end) = X(x_idx, y_idx,:);
C(:,2:2:end) = Y(x_idx, y_idx,:);


% Ensures at least one update
obj_old = inf;

% Update and assign until change in objective is below threshold or
% MAX_ITER is reached
in = 1;
%wts = ones(Nc,1);
while 1
    % Calculate distances
    D = dist_func(X,Y,T,C,varargin{:});

    % Calculate membership probabilities
    U = membership(D,m);
    
    % Evaluate objective function
    %tmp = weight_func(X,Y,C,U);
    %wts(:,in) = tmp(1,1,:);
    obj_new = sum(weight_func(X,Y,C,U).*U.*D, 1:3);
    
    % Check for insignificant improvement
    if (abs(obj_new - obj_old) < thres)
        fprintf("Reached threshold in %i iterations\n", in);
        break;
    end
    obj_old = obj_new;
    
    in = in + 1;
    if (in > MAX_ITER)
        fprintf("Warning: Maximum number of iterations %i reached. Solution may not have converged.\n", MAX_ITER)
        break
    end
    
    % Calculate centers
    C = centres(U,X,Y);
end


%figure;
%plot(1:in, wts)

end

function C = centres(U,X,Y)
    C = NaN(size(U,3), 2*size(X,3));
    Cnorm = sum(U,1:2);
    U = reshape(U, [], size(U,3));
    
    C(:,1:2:end) = U'*reshape(X,[],size(X,3))./squeeze(Cnorm);
    C(:,2:2:end) = U'*reshape(Y,[],size(Y,3))./squeeze(Cnorm);
end

function U = membership(D,m)
    U = NaN(size(D));
    for i = 1:size(D,1)
        for j = 1:size(D,2)
            U(i,j,:) = (1./(D(i,j,:).^(1/(m-1))))./(sum(1./(D(i,j,:).^(1/(m-1)))));
        end
    end
    U(isnan(U)) = 1;
end