% Testing kron for partial derivs


n = 100;

f = @(x,y)sin(x).*exp(x).*cos(y);
fx_e = @(x,y)(cos(y).*(cos(x).*exp(x)+sin(x).*exp(x)));


h = 2/(n-1);

x = -1:h:1;
[xx,yy] = ndgrid(x);
u = f(xx,yy);

e = ones(n,1);
z = zeros(n,1);
A1 = spdiags([-e z e],0:2,n-2,n)/(2*h);

A = kron(eye(n),A1);

tic
ux1 = (u(3:end,:)-u(1:end-2,:))/(2*h);
toc

tic
ux2 = reshape(A*u(:),[n-2,n]);
toc

ux_e = fx_e(xx,yy);

norm(ux2-ux1)
norm(ux2-ux_e(2:end-1,:))/norm(ux_e)