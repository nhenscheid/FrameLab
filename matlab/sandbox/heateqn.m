clear all;

Nx = 40;
Nt = 40;
dx = 1/(Nx-1);
dt = 1/(Nt-1);


x = (0:dx:1)';
t = (0:dt:1)';
[T,X] = ndgrid(t,x);

u = @(t,x)(exp(-t).*sin(x));

f0 = @(x)(sin(x));


Dxx = -spdiags([ones(Nx,1) -2*ones(Nx,1) ones(Nx,1)],0:2,Nx-2,Nx)/dx^2;
Dt = spdiags([-ones(Nt,1) zeros(Nt,1) ones(Nt,1)],0:2,Nt-2,Nt)/(2*dt);

trimx = speye(Nx); trimx = trimx(2:end-1,:);
trimx = kron(trimx,speye(Nt-2));
trimt = speye(Nt); trimt = trimt(2:end-1,:);
trimt = kron(speye(Nx-2),trimt);

D = trimt*kron(Dxx,speye(Nt)) + trimx*kron(speye(Nx),Dt);

y = f0(x);

yp = Dxx*y;

U = u(T,X);
Up = D*U(:);