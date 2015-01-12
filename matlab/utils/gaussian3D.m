function y = gaussian3D(N)
% Compute a 3D Gaussian function 

[X,Y,Z] = ndgrid(linspace(-1,1,N));

y = exp(-(X.^2+Y.^2+Z.^2));
