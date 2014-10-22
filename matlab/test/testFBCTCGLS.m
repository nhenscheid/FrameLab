% Testing Cone Beam reconstruction with cgls
clear all;close all;

lam = 1;



u = DataTypes.ObjectData(2,single(phantom(256)),[10,10]);
fbct = Operators.FanBeamScanner(256,64);
A = @(x)(fbct.applyAdjoint(fbct.apply(x))+lam*x);
f0 = fbct.apply(u);

rhs = fbct.applyAdjoint(f0); 

rk = rhs;

pk = rhs;

uk = DataTypes.ObjectData(2,single(zeros(256)),[10,10]);

gammaold = rk.objdot(rk);

icg = 0;
cgiter = 50;
cgtol = 1e-16;

while icg<cgiter
    icg = icg+1;
    qk = A(pk);
    alphak = gammaold/(pk.objdot(qk));
    uk = uk + alphak*pk;
    rk = rk + (-alphak*qk);
    gammanew = rk.objdot(rk)
    betak = gammanew/gammaold;
    pk = rk + betak*pk;
    gammaold = gammanew;
    if(gammaold<cgtol)
        icg = cgiter;
    end
    imshow(uk.dataArray,[]);
end