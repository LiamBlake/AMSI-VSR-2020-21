function [lavd,w] = LAVD(X,Y,u,v,tspan,step)
%LAVD Lagrangian-averaged vorticity deviation
%   Detailed explanation goes here

% Default step size
if nargin < 6
    step = 0.01;
end

% Calculate vorticity field 
w = NaN(size(X,1:3));

for i = 1:length(tspan)
    % Approximate derivatives with central difference
    xp = X(:,:,i);
    yp = Y(:,:,i);
    t = tspan(i)*ones(size(xp));
    dvdx = (v(t,xp + step,yp) - v(t,xp - step,yp))/(2*step);
    dudy = (u(t,xp,yp + step) - u(t,xp,yp - step))/(2*step);
    
    w(:,:,i) = dvdx - dudy;
end



% Integrate for LAVD
Nx = size(X,1); Ny = size(X,2);
wrs = reshape(w, Nx*Ny, length(tspan));

% Spatial mean of vorticity
wbar = mean(wrs', 2); 

% Integrate
lavd = trapz(tspan, abs(bsxfun(@minus, wrs', wbar)), 1);
lavd = reshape(lavd, Nx, Ny);

end

