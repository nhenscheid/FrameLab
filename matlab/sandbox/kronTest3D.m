% Testing kron for partial derivs

clear all;
n = 128;
f = @(x,y,z)sin(x).*exp(z).*cos(y);
fx_e = @(x,y,z)cos(x).*exp(z).*cos(y);
fy_e = @(x,y,z)(-sin(y).*sin(x).*exp(z));
fxx_e = @(x,y,z)(-sin(x).*exp(z).*cos(y));
fyy_e = @(x,y,z)(-sin(x).*sin(y).*exp(z));
fxy_e = @(x,y,z)-cos(x).*sin(y).*exp(z);
fxyz_e = @(x,y,z)-cos(x).*exp(z).*sin(y);
fyz_e = @(x,y,z)-sin(y).*sin(x).*exp(z);
h = 2/(n-1);
x = -1:h:1;
[xx,yy,zz] = ndgrid(x);
u = f(xx,yy,zz);

%%
% first partials
e = ones(n,1);
z = zeros(n,1);
A1 = spdiags([-e z e],0:2,n-2,n)/(2*h);
A2 = kron(eye(n),A1);
A = kron(eye(n),A2);

disp('computing partial-x using array notation')
tic
ux1 = (u(3:end,:,:)-u(1:end-2,:,:))/(2*h);
toc

disp('computing partial-x using kron matrix')
tic
ux2 = reshape(A*u(:),[n-2,n,n]);
toc

disp('Computing exact solution')
ux_e = fx_e(xx,yy,zz);

disp('max error between kron and array')
max(abs(ux2(:)-ux1(:)))
ux_e2 = ux_e(2:end-1,:,:);
disp('2-norm error from exact')
norm(ux2(:)-ux_e2(:))/norm(ux_e2(:))


disp('computing partial y')

up = permute(u,[2,1,3]);

disp('computing partial-y using array notation')
tic
uy1 = (u(:,3:end,:)-u(:,1:end-2,:))/(2*h);
toc

disp('computing partial-y with permuted array and kron')
tic
uy2 = ipermute(reshape(A*up(:),[n-2,n,n]),[2,1,3]);
toc

uy_e = fy_e(xx,yy,zz);

disp('max error between kron and array')
max(abs(uy2(:)-uy1(:)))
uy_e2 = uy_e(:,2:end-1,:);
disp('2-norm error from exact')
norm(uy2(:)-uy_e2(:))/norm(uy_e2(:))

% Testing second derivatives using kron 

disp('Testing second derivatives')
e = ones(n,1);
z = zeros(n,1);
A1 = spdiags([e,z,-2*e,z,e],0:4,n-4,n);
A2 = kron(eye(n),A1);
A = kron(eye(n),A2)/(4*h^2);

disp('Array method')
tic
uxx1 = (u(5:end,:,:)+u(1:end-4,:,:)-2*u(3:end-2,:,:))/(4*h^2);
toc

disp('Kron method')
tic
uxx2 = reshape(A*u(:),[n-4,n,n]);
toc


uxx_e = fxx_e(xx,yy,zz);
uxx_e2 = uxx_e(3:end-2,:,:);

disp('max error between kron and array')
max(abs(uxx2(:)-uxx1(:)))
disp('L2 error from exact')
norm(uxx_e2(:)-uxx2(:))/norm(uxx_e2(:))


%%
disp('Testing second y derivative')
e = ones(n,1);
z = zeros(n,1);
A1 = spdiags([e,z,-2*e,z,e],0:4,n-4,n);
A2 = kron(A1,eye(n));
A = kron(eye(n),A2)/(4*h^2);

disp('Array method')
tic
uyy1 = (u(:,5:end,:)+u(:,1:end-4,:)-2*u(:,3:end-2,:))/(4*h^2);
toc

disp('Kron method')
tic
uyy2 = reshape(A*u(:),[n,n-4,n]);
toc


uyy_e = fxx_e(xx,yy,zz);
uyy_e2 = uyy_e(:,3:end-2,:);

disp('max error between kron and array')
max(abs(uyy2(:)-uyy1(:)))
disp('L2 error from exact')
norm(uyy_e2(:)-uyy1(:))/norm(uyy_e2(:))


%%
disp('Testing xy derivative')
e = ones(n,1);
z = zeros(n,1);
A1 = spdiags([-e,z,e],0:2,n-2,n);
A2 = kron(A1,A1);
A = kron(eye(n),A2)/(4*h^2);

disp('Array method')
tic
uxy1 = (u(3:end,3:end,:)+u(1:end-2,1:end-2,:)-u(3:end,1:end-2,:)-u(1:end-2,3:end,:))/(4*h^2);
toc

disp('Kron method')
tic
uxy2 = reshape(A*u(:),[n-2,n-2,n]);
toc


uxy_e = fxy_e(xx,yy,zz);
uxy_e2 = uxy_e(2:end-1,2:end-1,:);

disp('max error between kron and array')
max(abs(uxy2(:)-uxy1(:)))
disp('L2 error from exact')
norm(uxy_e2(:)-uxy1(:))/norm(uxy_e2(:))


%%
disp('Testing xyz derivative')
e = ones(n,1);
z = zeros(n,1);
A1 = spdiags([-e,z,e],0:2,n-2,n);
A2 = kron(A1,A1);
A = kron(A1,A2)/(8*h^3);

disp('Kron method')
tic
uxyz2 = reshape(A*u(:),[n-2,n-2,n-2]);
toc


uxyz_e = fxyz_e(xx,yy,zz);
uxyz_e2 = uxyz_e(2:end-1,2:end-1,2:end-1);

disp('L2 error from exact')
norm(uxyz_e2(:)-uxyz2(:))/norm(uxyz_e2(:))


%%
disp('Testing yz derivative')
e = ones(n,1);
z = zeros(n,1);
A1 = spdiags([-e,z,e],0:2,n-2,n);
A2 = kron(A1,eye(n));
A = kron(A1,A2);

disp('Kron method')
tic
uyz2 = reshape(A*u(:),[n,n-2,n-2])/(4*h^2);
toc


uyz_e = fyz_e(xx,yy,zz);
uyz_e2 = uyz_e(:,2:end-1,2:end-1);

disp('L2 error from exact')
norm(uyz_e2(:)-uyz2(:))/norm(uyz_e2(:))
