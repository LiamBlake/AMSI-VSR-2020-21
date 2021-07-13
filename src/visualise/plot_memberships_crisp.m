function plot_memberships_crisp(mu, suptitle,xval, yval, xlab, ylab)
%PLOT_MEMBERSHIPS 

smu = squeeze(mu);


K = size(smu,3);

% Determine indices
idx = NaN(size(smu,1), size(smu,2));
vals = 1:K;
for i = 1:size(smu,1)
    for j = 1:size(smu,2)
        if sum(isnan(smu(i,j,:))) > 0
            continue
        end
        idx(i,j) = vals(logical(smu(i,j,:)));
    end
end

% Figure
figure;
colourplot(xval, yval, idx, suptitle, xlab, ylab, "", colorcube)
colorbar('off')


end