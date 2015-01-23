% TWF regularized data interpolation with John's equation

clear all;

% Generate object
N = 256;
u0 = DataTypes.ObjectData(3,single(unitBall(N,3)),[2,2,2]); %unit ball, diameter = 2
lam = 1;
cgiter = 10;
cgtol = 1e-8;
%u1 = DataTypes.ObjectData(3,single(gaussian3D(N)),[5,5,5]);
%u0 = DataTypes.ObjectData(3,single(phantom3d(N)),[2,2,2]);

% Set up the scanner
nd = 256;
rps = 1;
fps = 128;
zmax = 1;
vtab = 1;
nHelix = 1;
dphi = 2*pi*rps/fps;
phaseShift = dphi*(0:nHelix-1);
cbct = Operators.ConeBeamScanner('multiHelix',nd,nd,[],zmax,rps,vtab,fps,nHelix,phaseShift);
cbct.verbose = true;
cbct.GPU = 1;

% Define local handles for operators
RL = @(x)changeScanner(x,u0,'down');
RLt = @(x)changeScanner(x,u0,'up');
D = @(x)(x.applyJohn(0));
Dt = @(x)(x.applyJohnAdjoint(0));
A = @(x)(RLt(RL(x))+lam*Dt(D(x)));

% Compute forward transform of u0
disp('Computing single helix scan for u0');
f0 = cbct.apply(u0);

% Interpolated f 
f = RLt(f0);
b = RLt(f0);

% CGLS
    %rhs = obj.At(obj.b);
    rk = A(f)-b;
    pk = -rk;
    uk = f;
    gammaold = rk.objdot(rk); %object rk must have dot product method
    icg = 0;
    cgiter = obj.globalIter;
    cgtol = obj.convergenceTol;
    % Main conjugate gradient loop
    while icg<cgiter
        icg = icg+1;
        qk = A(pk);
        alphak = gammaold/(pk.objdot(qk));
        uk = uk + alphak*pk;
        rk = rk + alphak*qk;
        gammanew = rk.objdot(rk)
        betak = gammanew/gammaold;
        pk = -rk + betak*pk;
        gammaold = gammanew;
        if(gammaold<cgtol)
            icg = cgiter;
        end
    end
    y = uk;