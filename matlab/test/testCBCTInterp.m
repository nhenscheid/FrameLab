% TWF regularized data interpolation with John's equation

clear all;

% Generate object
N = 256;
u0 = DataTypes.ObjectData(3,single(unitBall(N,3)),[2,2,2]); %unit ball, diameter = 2
lam = 0.125;
cgiter = 10;
cgtol = 1e-8;
%u1 = DataTypes.ObjectData(3,single(gaussian3D(N)),[5,5,5]);
%u0 = DataTypes.ObjectData(3,single(phantom3d(N)),[2,2,2]);

% Set up the scanner
nd = 64;
rps = 1;
fps = 16;
zmax = 1;
vtab = 1;
nHelix = 1;
dphi = 2*pi*rps/fps;
phaseShift = dphi*(0:nHelix-1);
cbct = Operators.ConeBeamScanner('multiHelix',nd,nd,[],zmax,rps,vtab,fps,nHelix,phaseShift);
nHelix2 = 2;
phaseShift2 = dphi*(0:nHelix2-1);
cbct2 = Operators.ConeBeamScanner('multiHelix',nd,nd,[],zmax,rps,vtab,fps,nHelix2,phaseShift2);
cbct.verbose = true;
cbct.GPU = 1;

% Define local handles for operators
RL = @(x)changeScanner(x,u0,'down');
RLt = @(x)changeScanner(x,u0,'up');
DtD = @(x)(x.applyJohnNormal(0));
A = @(x)(lam*RLt(RL(x))+DtD(x));

% Compute forward transform of u0
disp('Computing single helix scan for u0');
f0 = cbct.apply(u0);
f_exact = cbct2.apply(u0);

% Interpolated f 
f = RLt(f0);
b = RLt(f0);

iMain = 0;
mainIter = 10;
% CGLS
    while iMain < mainIter
        rk = A(f)-b;
        pk = -rk;
        %uk = f;
        gammaold = rk.objdot(rk); %object rk must have dot product method
        icg = 0;
        %cgiter = obj.globalIter;
        %cgtol = obj.convergenceTol;
        % Main conjugate gradient loop
        while icg<cgiter
            icg = icg+1;
            qk = A(pk);
            alphak = gammaold/(pk.objdot(qk));
            f = f + alphak*pk;
            rk = rk + alphak*qk;
            gammanew = rk.objdot(rk)
            betak = gammanew/gammaold;
            pk = -rk + betak*pk;
            gammaold = gammanew;
            if(gammaold<cgtol)
                icg = cgiter;
            end
            disp('Error from true')
            norm(f.dataArray(:)-f_exact.dataArray(:))/norm(f_exact.dataArray(:))
        end
        lam = lam*2
    end
    
