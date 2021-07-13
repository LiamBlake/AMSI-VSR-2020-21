function Dkl = relative_entropy(P,Q)
%RELATIVE_ENTROPY Summary of this function goes here
%   Detailed explanation goes here

terms = P.*log(P./Q);
terms(isnan(terms)) = 0;
Dkl = sum(terms);

end

