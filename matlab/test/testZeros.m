% Test the cone beam scanner with a matrix of zeros (useful for making sure
% the coordinates are correct!

clear all;clc
N = 128;
u0 = DataTypes.ObjectData(3,single(zeros(N,N,N)),[2,2,2]);
nd = 16;

cbct = Operators.ConeBeamScanner('helix',nd,nd,[],2,2,1,10);

cbct.verbose = true;
cbct.GPU = 0;


y = cbct.apply(u0);




