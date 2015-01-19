% Testing John's equation, helix parametrization

clear all;

% Generate object
N = 64;
u0 = DataTypes.ObjectData(3,single(unitBall(N,3)),[2,2,2]); %unit ball, diameter = 2
u1 = DataTypes.ObjectData(3,single(gaussian3D(N)),[7,7,7]);
%u0 = DataTypes.ObjectData(3,single(phantom3d(N)),[2,2,2]);

% Set up the scanner
nd = 256;
rps = 1;
fps = 256;
zmax = 0.25;
vtab = 0.25;
nHelix = 2;
dphi = 2*pi*rps/fps;
phaseShift = [0,dphi];
cbct.verbose = true;
%cbct.GPU = 1;
cbct = Operators.ConeBeamScanner('multiHelix',nd,nd,[],zmax,rps,vtab,fps,nHelix,phaseShift);

% Compute forward transform of u0
disp('Computing multihelix scan for u0');
f = cbct.apply(u0);
disp('Computing multihelix scan for u1');
g = cbct.apply(u1);

% Apply John's equation
D = f.applyJohn();

Dt = g.applyJohnAdjoint();


y = g.dataArrayNorm(2:end-1,2:end-1,2:end-1,1);
x = f.dataArrayNorm(2:end-1,2:end-1,2:end-1,1);

abs(dot(D(:),y(:))-dot(Dt(:),x(:)))/abs(dot(D(:),y(:)))