function [l1,l2,sp,sm] = strain_eigs(C)
%STRAIN_EIGENVECTORS Eigenvalues of Cauchy-Green tensor
% usage: [l1,l2,sp,sm] = strain_eigenvectors(C)
%   
% 

Nx = size(C,1);
Ny = size(C,2);

l1 = NaN(Nx,Ny);
l2 = NaN(Nx,Ny);
sp = NaN(Nx,Ny,2);
sm = NaN(Nx,Ny,2);

for i = 1:Nx
	for j = 1:Ny
		[v,D] = eigs(squeeze(C(i,j,:,:)));
        	[d,ind] = sort(diag(D));
		v = v(:,ind);
        
		% Shrink and stretch lines
		l1(i,j) = d(1);
		l2(i,j) = d(2);

		% Shear fields
		t1 = sqrt(sqrt(d(2))/(sqrt(d(1)) + sqrt(d(2))))*v(:,1);
		t2 = sqrt(sqrt(d(1))/(sqrt(d(1)) + sqrt(d(2))))*v(:,2);
        	sp(i,j,:) = t1 + t2;
		sm(i,j,:) = t1 - t2;
                   
	end
end


end

