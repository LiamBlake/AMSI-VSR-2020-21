function [M,k] = contour_membership(xi,yi,bnd)
%CONTOUR_MEMBERSHIP Summary of this function goes here
%   Detailed explanation goes here
%
%   arguments:
%       xi - meshgrid of initial x-components
%       yi - meshgrid of initial y-components

% No. closed contour regions - assumes that all vortex regions are disjoint
% - CHECK THIS ASSUMPTION
nc = size(bnd.xc,2);
k = nc + 1;

M = zeros(size(xi,1), size(xi,2),k);
on = zeros(size(xi,1), size(xi,2),k);
for n = 1:nc
    % Interpret boundary as a polygon (with a very large number of sides)
    xv = bnd.xc{n};
    yv = bnd.yc{n};
    
    % Determine which points are inside the polygon
    [M(:,:,n),on(:,:,n)] = inpolygon(xi,yi,xv,yv);
    
end


% Incoherent background
incoh = repmat(sum(M,3) == 0, 1,1, nc+1);
incoh(:,:,1:nc) = 0;
M(incoh) = 1;


% Ensure no overlap


end

