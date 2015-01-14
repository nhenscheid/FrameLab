% Test double helix Geometry 

clear all;clc
N = 64;
u0 = DataTypes.ObjectData(3,single(phantom3d(N)),[2,2,2]);
nd = 64;
nv = 64;
phase = pi/4;


disp('Testing double helix scan');
cbct = Operators.ConeBeamScanner('doubleHelix',nd,nd,[],2,4,1,nv,phase);
cbct.verbose = true;
cbct.GPU = 1;

yhelix = cbct.apply(u0);


%disp('Computing AtAu');

%AtAuHelix = cbct.applyAdjoint(yhelix);