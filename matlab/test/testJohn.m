% Testing John's equation, helix parametrization

clear all;

% Generate object
N = 384;
u0 = DataTypes.ObjectData(3,single(unitBall(N,3)),[2,2,2]); %unit ball, diameter = 2
%u0 = DataTypes.ObjectData(3,single(gaussian3D(N)),[5,5,5]);
%u0 = DataTypes.ObjectData(3,single(phantom3d(N)),[2,2,2]);

% Set up the scanner
nd = 128;
rps = 0.25;
fps = 128;
zmax = 0.5;
vtab = 0.5;
nHelix = 2;
dphi = 2*pi*rps/fps;
phaseShift = dphi*(0:nHelix-1);
%phaseShift = [0 dphi];
cbct = Operators.ConeBeamScanner('multiHelix',nd,nd,[],zmax,rps,vtab,fps,nHelix,phaseShift);
cbct.verbose = true;
cbct.GPU = 1;

% Compute forward transform of u0
%disp('Computing multihelix scan for u0');
f = cbct.apply(u0);
%disp('Computing multihelix scan for u1');
%g = cbct.apply(u1);

% Apply John's equation
% D = f.applyJohn(0);
% 
% Dt = g.applyJohnAdjoint(0);
% 
% 
% y = g.dataArrayNorm(2:end-1,2:end-1,2:end-1,1);
% x = f.dataArrayNorm(2:end-1,2:end-1,2:end-1,1);
% 
% abs(dot(D(:),y(:))-dot(Dt(:),x(:)))
% 
%disp(sprintf('%s%d\n','max(D)=',max(D(:))))
y_det = cbct.para.y_det*cbct.para.scale;
z_det = cbct.para.z_det*cbct.para.scale;
sd_phi = cbct.para.sd_phi;
h = cbct.vtab/(2*pi*cbct.rps); 
nv = cbct.nv/cbct.nHelix;
na = cbct.na;
nb = cbct.nb;
yExact = zeros(na,nb,nv,nHelix);
yExactNorm = zeros(na,nb,nv,nHelix);
for i=1:nHelix
    [A,B,T] = ndgrid(y_det,z_det,sd_phi((i-1)*nv+1:i*nv));
    zeta = -phaseShift(i)*h;
    yE = coneBeamSphere(T,A,B,zeta,h,cbct.SO,cbct.OD,0);
    yExact(:,:,:,i) = yE;
    yEN = coneBeamSphere(T,A,B,zeta,h,cbct.SO,cbct.OD,1);
    yExactNorm(:,:,:,i) = yEN;
    yApprox = f.dataArray(:,:,:,i);
    yApproxNorm = f.dataArrayNorm(:,:,:,i);
    disp(sprintf('%s%d\n','L2 error for unnormalized transform is ',sqrt(sum(abs(yE(:)-yApprox(:)).^2))/sqrt(sum(abs(yE(:))))))
    disp(sprintf('%s%d\n','L2 error for normalized transform is ',sqrt(sum(abs(yEN(:)-yApproxNorm(:)).^2))/sqrt(sum(abs(yE(:))))))
    disp(sprintf('%s%d\n','max error for unnormalized transform is ',max(abs(yE(:)-yApprox(:)))/max(abs(yE(:)))))
    disp(sprintf('%s%d\n','max error for normalized transform is ',max(abs(yEN(:)-yApproxNorm(:)))/max(abs(yEN(:)))))
end

% 
y = DataTypes.CTData(cbct,yExact,yExactNorm,[2,2,2]);
%D = g.applyJohnNormal(0);
D1 = y.applyJohn(1);
%D2 = g.applyJohn(0);
%size(f.dataArrayNorm)
D2 = f.applyJohn(1);