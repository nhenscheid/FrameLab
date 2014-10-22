u = DataTypes.ObjectData(2,single(phantom(256)),[10,10]);
fbct = Operators.FanBeamScanner(256,16)
y = fbct.apply(u)

Aty = fbct.applyAdjoint(y)
