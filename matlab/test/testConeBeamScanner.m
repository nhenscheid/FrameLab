clear u cbct y Aty;

u = DataTypes.ObjectData(3,single(phantom3d(128)),[10,10,10]);
cbct = Operators.ConeBeamScanner(256,256,128);
disp('computing forward transform')
y = cbct.apply(u);
disp('computing adjoint transform')
Aty = cbct.applyAdjoint(y);

%% Test on a very small dataset
clc; clear u cbct y Aty;
u = DataTypes.ObjectData(3,single(zeros(10,10,10)),[1,1,1]);
cbct = Operators.ConeBeamScanner(4,4,4);
cbct.GPU = 0;
cbct.verbose = true;
disp('computing forward transform')
y = cbct.apply(u);
disp('computing adjoint transform')
Aty = cbct.applyAdjoint(y);

