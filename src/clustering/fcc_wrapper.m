function [mu,v] = fcc_wrapper(membership,Ki,wi,f,K,m)


if (f == "Euclidean") && (sum(diff(wi)) == 0)
    % Degenerate case - reduces to regular FCM
    fprintf("Degenerate case - reduces to regular fuzzy c-means\n");
    
    % Construct multimembership data
    Y = [];
    for i = 1:length(membership)
        Y = [Y membership{i}];
    end
    
    [mu,v] = fcm_wrapper(Y,K,m,0);
    
else
    if f == "Euclidean"
        % Distance function - norm
        f = @(x,y) vecnorm(x - y,2,2).^2;
    end
    % TODO: implement other distance norms
    
    [mu,v] = fc_consensus(membership,Ki,wi,f,K,m);
end


end