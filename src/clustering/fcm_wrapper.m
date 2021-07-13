function [U,c] = fcm_wrapper(X,k,m,print)
% FCM_WRAPPER Fuzzy c-means wrapper with k-means

% Default arguments
if nargin < 4
    print = 0;
end

if m == 1
    % k-means clustering
    if print
        passstr = "iter";
    else
        passstr = "final";
    end
    
    [idx,c] = kmeans(X,k,'Display',passstr);
    U = zeros(size(X,1), k);
    
    for i = 1:size(X,1)
        U(i,idx(i)) = 1; 
    end
    
else
    % fuzzy c-means clustering
    opts = [2, NaN, NaN, print];
    [c,U] = fcm(X, k, opts);
    U = U';
end

end