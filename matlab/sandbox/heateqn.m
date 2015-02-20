Nx = 50;
Nt = 50;
dx = 1/(Nx-1);
dt = 1/(Nt-1);


x = (0:dx:1)';
t = (0:dt:1)';
[X,Y] = ndgrid(x,y);

u = @(x,t)(exp(-t).*sin(x));

f0 = @(x)(sin(x));


Dxx = -spdiags([ones(Nx,1) -2*ones(Nx,1) ones(Nx,1)],0:2,Nx-2,Nx)/dx^2;
Dt = spdiags([-ones(Nt,1) ones(Nt,1)],0:1,Nt-1,Nt)/dt;

trimx = speye(Nx); trimx = trimx(2:end-1,:);


D = kron(Dxx,speye(Nt))+kron(speye(Nx),Dt);

y = f0(x);

yp = Dxx*y;

U = u(X,Y);
Up = D*U(:);