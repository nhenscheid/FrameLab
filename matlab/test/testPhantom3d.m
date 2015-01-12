% Test the cone beam scanner with a unit ball
clear all;clc
N = 128;
u0 = DataTypes.ObjectData(3,single(phantom3d(N)),[2,2,2]);
nd = 128;
nv = 64;

disp('Testing circular scan');
cbct = Operators.ConeBeamScanner('circle',nd,nd,nv);

cbct.verbose = true;
cbct.GPU = 1;

ycircle = cbct.apply(u0);

AtAuCircle = cbct.applyAdjoint(ycircle);


disp('Testing helical scan');
cbct = Operators.ConeBeamScanner('helix',nd,nd,[],2,4,1,128);
cbct.verbose = true;
cbct.GPU = 1;

yhelix = cbct.apply(u0);



disp('Computing AtAu');

AtAuHelix = cbct.applyAdjoint(yhelix);