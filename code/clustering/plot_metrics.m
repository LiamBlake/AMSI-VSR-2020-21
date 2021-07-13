function [ent, mpc] = plot_metrics(U_allk, part_title, kspan, ylab)
%PLOT_METRICS Plot fuzzy clustering performance metrics for a range of
%cluster numbers.
%   Detailed explanation goes here

% Default arguments
if nargin < 4
    ylab = "No. Clusters";
end
if nargin < 2
    kspan = 2:size(U_allk,1);
end

% Remove redundant dimensions
U_allk = squeeze(U_allk);

% Calculate metrics for each cluster configuration
ent = NaN(size(U_allk,1), size(U_allk,2), size(U_allk,3));
mpc = NaN(length(kspan));
% Add any additional metrics here as required
for i = 1:length(kspan)
    ent(i,:,:) = clust_entropy(U_allk(i,:,:,:));
    %mpc(i) = mp_coeff(U_allk(i,:,:,:));
end

% Plot mean metrics
plot(kspan, mean(ent, [2 3], 'omitnan'), 'r-')
hold on;
plot(kspan, mpc, 'b-')

%legend('Mean entropy');%, 'Modified partition coefficient');
title("Mean Entropy for " + part_title, 'Interpreter', 'LaTeX')
xlabel(ylab);
ylabel("Mean Entropy");


end

