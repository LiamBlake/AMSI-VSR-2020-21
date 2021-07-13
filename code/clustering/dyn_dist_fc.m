function D = dyn_dist_fc(X,Y,T,C)
    D = NaN(size(X,1), size(X,2), size(C,1));
    tk = T(1:(end-1));
    tk1 = T(2:end);
    tf = 1/(T(end) - T(1));
    
    nc = size(C,3);
    
    % Centres
    ck = C(:,1:(end-2));
    ck1 = C(:,3:end);
    
    for i = 1:size(D,1)
        for j = 1:size(D,2)
            xk = repmat(squeeze(X(i,j,1:(end-1)))', nc,1);
            xk1 = repmat(squeeze(X(i,j,2:end))', nc,1);
            yk = repmat(squeeze(Y(i,j,1:(end-1)))', nc,1);
            yk1 = repmat(squeeze(Y(i,j,2:end))', nc,1);
            D(i,j,:) = (tf*sum(0.5*(tk1-tk).*(((xk1 - ck1(:,1:2:end)).^2 + (yk1 - ck1(:,2:2:end)).^2).^(0.5) + ((xk - ck(:,1:2:end)).^2 + (yk - ck(:,2:2:end)).^2).^(0.5)),2));
        end
    end
    D = D.^2;
end