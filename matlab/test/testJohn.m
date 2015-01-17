% Testing John's equation, helix parametrization

clear all;

% Generate object
N = 64;
u0 = DataTypes.ObjectData(3,single(unitBall(N,3)),[2,2,2]); %unit ball, diameter = 2
%u0 = DataTypes.ObjectData(3,single(gaussian3D(N)),[7,7,7]);
%u0 = DataTypes.ObjectData(3,single(phantom3d(N)),[2,2,2]);

% Set up the scanner
nd = 128;
rps = 1;
fps = 64;
zmax = 1;
vtab = 1;
nHelix = 2;
dphi = 2*pi*rps/fps;
phaseShift = [0,dphi];
cbct.verbose = true;
%cbct.GPU = 1;
cbct = Operators.ConeBeamScanner('multiHelix',nd,nd,[],zmax,rps,vtab,fps,nHelix,phaseShift);

% Compute forward transform of u0
disp('Computing multihelix scan');
f = cbct.apply(u0);

% Apply John's equation
D = f.applyJohn();