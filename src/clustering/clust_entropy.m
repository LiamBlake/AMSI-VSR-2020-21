function h = clust_entropy(U)
%CLUST_ENTROPY Calculates normalised entropy of the probability vector of
%membership probabilities for a given clustering.
%   Detailed explanation goes here

% Remove any redundant dimensions and recurrent NaN entries
U = squeeze(U);
nNan = max(sum(~isnan(U),3), [], 'all');
U = U(:,:,1:nNan);


% Sizes
dim = size(U);
K = dim(end);

lU = log(U);
% Correct for zero probabilities
lU(U == 0) = 0;

% Entropy
h = -sum(U.*lU,length(dim), 'omitnan')./log(K);

end

