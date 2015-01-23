% Testing John's equation, helix parametrization

clear all;

% Generate object
N = 512;
u0 = DataTypes.ObjectData(3,single(unitBall(N,3)),[2,2,2]); %unit ball, diameter = 2
u1 = DataTypes.ObjectData(3,single(gaussian3D(N)),[5,5,5]);
%u0 = DataTypes.ObjectData(3,single(phantom3d(N)),[2,2,2]);

% Set up the scanner
nd = 256;
rps = 1;
fps = 128;
zmax = 1;
vtab = 1;
nHelix = 3;
dphi = 2*pi*rps/fps;
phaseShift = dphi*(0:nHelix-1);
cbct = Operators.ConeBeamScanner('multiHelix',nd,nd,[],zmax,rps,vtab,fps,nHelix,phaseShift);
cbct.verbose = true;
cbct.GPU = 1;

% Compute forward transform of u0
disp('Computing multihelix scan for u0');
f = cbct.apply(u0);
disp('Computing multihelix scan for u1');
g = cbct.apply(u1);

% Apply John's equation
D = f.applyJohn(0);

Dt = g.applyJohnAdjoint(0);


y = g.dataArrayNorm(2:end-1,2:end-1,2:end-1,1);
x = f.dataArrayNorm(2:end-1,2:end-1,2:end-1,1);

abs(dot(D(:),y(:))-dot(Dt(:),x(:)))

disp(sprintf('%s%d\n','max(D)=',max(D(:))))
