u = DataTypes.ObjectData(2,single(phantom(256)),[10,10]);
fbct = Operators.FanBeamScanner(256,256);
y = fbct.apply(u);

%AtAu = fbct.applyAdjoint(y);