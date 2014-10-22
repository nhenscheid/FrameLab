clear u cbct y Aty;

u = DataTypes.ObjectData(3,single(phantom3d(128)),[12,12,12]);
cbct = Operators.ConeBeamScanner(256,256,128);
disp('computing forward transform')
y = cbct.apply(u);
disp('computing adjoint transform')
Aty = cbct.applyAdjoint(y);
