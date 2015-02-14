% Solving a transport equation using Augmented Lagrangian?

clear all;
close all;
Nx = 100;
dx = 1/(Nx-1);

Nt = 100;
dt = 1/(Nt-1);

x = 0:dx:1;
t = 0:dt:1;

f = @(x)sin(2*pi*x);

u = zeros(Nt,Nx);

u(1,:) = f(x);

plot(u(1,:))


for i = 2:Nt
    u(i,1:end-1) = u(i-1,1:end-1) + (dt/dx)*(u(i-1,2:end)-u(i-1,1:end-1));
    u(i,end) = 0;
    plot(u(i,:)),drawnow,pause(0.001)
end


ex = ones(Nx,1);
et = ones(Nt,1);
Dx = spdiags([-ex ex],0:1,Nx,Nx)/dx;
Dx(end,end) = 1;
Dt = spdiags([-et et],0:1,Nt,Nt)/dt;

Dx = kron(Dx,speye(Nt));
Dt = kron(speye(Nx),Dt);


Dxu = reshape(Dx*u(:),[Nt,Nx]);

Dtu = reshape(Dt*u(:),[Nt,Nx]);

norm(Dtu(:)-Dxu(:))


v=reshape((Dt-Dx)\zeros(Nt*Nx,1),[Nt,Nx]);