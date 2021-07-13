function ftle = FTLE(C,tdiff)
%FTLE Calculates finite-time Lypanouv exponent given a 4D array of
%Cauchy-Green tensors.
%   Detailed explanation goes here

% Calculate eigenvalues
max_eigs = NaN(size(C,1:2));
for i = 1:size(C,1)
    for j = 1:size(C,2)
        max_eigs(i,j) = max(eig(squeeze(C(i,j,:,:))));
    end
end

ftle = 1/(2*tdiff) * log(max_eigs);

end

