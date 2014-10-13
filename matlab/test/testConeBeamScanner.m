clear u cbct y Aty;

u = DataTypes.ObjectData(3,single(phantom3d(256)),[10,10,10]);
cbct = Operators.ConeBeamScanner();
y = cbct.apply(u);
Aty = cbct.applyAdjoint(y);
