% Testing Cone Beam reconstruction with cgls
clear all;

lam = 0.001;



u = DataTypes.ObjectData(3,single(phantom3d(128)),[10,10,10]);
cbct = Operators.ConeBeamScanner(256,256,64);
A = @(x)cbct.applyAdjoint(cbct.apply(x))+lam*x;
f0 = cbct.apply(u);

rhs = cbct.applyAdjoint(f0); 

rk = -rhs;

pk = rhs;

uk = DataTypes.ObjectData(3,single(zeros(128,128,128)),[10,10,10]);

gammaold = rk.objdot(rk);

icg = 0;
cgiter = 50;
cgtol = 1e-16;

while icg<cgiter
    icg = icg+1
    qk = A(pk);
    alphak = gammaold/pk.objdot(qk);
    uk = uk + alphak*pk
    rk = rk + alphak*qk
    gammanew = rk.objdot(rk)
    betak = gammanew/gammaold;
    pk = -rk + betak*pk;
    gammaold = gammanew;
    if(gammaold<cgtol)
        icg = cgiter;
    end
    uk.plotData(3,3,10);
end