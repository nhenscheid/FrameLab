% Testing Cone Beam reconstruction with cgls
clear all;

lam = 0.01;


nd = 64;
rps = 1;
fps = 128;
nv = 64;
zmax = 1;
vtab = 1;
NTrue = 64;
NRecon = 64;

%cbct = Operators.ConeBeamScanner('helix',nd,nd,[],zmax,rps,vtab,fps);
cbct = Operators.ConeBeamScanner('circle',nd,nd,nv);

u = DataTypes.ObjectData(3,single(phantom3d(NTrue)),[2,2,2]);


A = @(x)(cbct.applyAdjoint(cbct.apply(x),[NRecon,NRecon,NRecon])+lam*x);
f0 = cbct.apply(u);

rhs = cbct.applyAdjoint(f0,[NRecon,NRecon,NRecon]); 

rk = -rhs;

pk = rhs;

uk = DataTypes.ObjectData(3,single(zeros(NRecon,NRecon,NRecon)),[2,2,2]);

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