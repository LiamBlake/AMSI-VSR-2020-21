function cbar = colourplot_land(x,y,C,stitle,xlab,ylab,clab,Znan,cmap)
%COLOURPLOT_LAND Summary of this function goes here
%   Detailed explanation goes here

% Default arguments
if nargin < 9
    cmap = parula;
end

Nx = length(x);
Ny = length(y);

% NAN as required
Z = reshape(C,Nx*Ny,[]);
Z(Znan) = NaN;

cbar = colourplot(x,y,reshape(Z,Nx,Ny,[]),stitle,xlab,ylab, clab,cmap);

end

