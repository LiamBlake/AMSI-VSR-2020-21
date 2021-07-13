function U = alpha_ftle_lavd(ftle_U,lavd_U,k,m,f,alpha)


K1 = size(ftle_U,2);
K2 = size(lavd_U,2);

w = [alpha 1 - alpha];

membership = {ftle_U, lavd_U};

U = fc_consensus(membership,[K1 K2],w,f,k,m);


end