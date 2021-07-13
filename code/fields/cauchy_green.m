function C = cauchy_green(X,Y,i,j)%(X,Y,grid, i,j)
%CAUCHY_GREEN Calculates finite difference approximation of the
%Cauchy-Green tensor for a fluid flow. Assume uniform grid of initial
%conditions

% Use smallest distance to next elements in discrete grid
%dx = abs(grid(i+1,j,1) - grid(i,j,1));
%dy = abs(grid(i,j+1,2) - grid(i,j,2));
dx = abs(X(i+1,j,1) - X(i,j,1));
dy = abs(Y(i,j+1,1) - Y(i,j,1));

% Approximate Jacobian derivative with finite difference
DF = [(X(i+1,j,end) - X(i-1,j,end))/(2*dx) (X(i,j+1,end) - X(i,j-1,end))/(2*dy);
      (Y(i+1,j,end) - Y(i-1,j,end))/(2*dx) (Y(i,j+1,end) - Y(i,j-1,end))/(2*dy)];

% Calculate Cauchy-Green tensor
C = DF'*DF;

end