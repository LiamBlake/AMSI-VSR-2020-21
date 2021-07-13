function value = mp_coeff(U)
%MP_COEFF Modified partition coefficient
%   Calculates the modified partition coefficient for a 

% Remove redundant dimensions and recurrent NaN entries 
U = squeeze(U);
nNan = max(sum(~isnan(U),3), [], 'all');
U = U(:,:,1:nNan);

% Calculate standard partition coefficient
pc = 1/(size(U,1)*size(U,2)) * sum(U.^2, 'all', 'omitnan');

% Modified to lie within [0,1]
c = size(U,3);
value = 1 - c/(c-1) * (1 - pc);

end

