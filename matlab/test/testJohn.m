% Testing John's equation, helix parametrization

clear all;

% Generate object
N = 256;
%u0 = DataTypes.ObjectData(3,single(unitBall(N,3)),[2,2,2]); %unit ball, diameter = 2
%u0 = DataTypes.ObjectData(3,single(gaussian3D(N)),[8,8,8]);
u0 = DataTypes.ObjectData(3,single(phantom3d(N)),[2,2,2]);

% Set up the scanner
nd = 128;
rps = 1;
fps = 128;
zmax = 2;
vtab = 1;
nHelix = 2;
dphi = 2*pi*rps/fps;
phaseShift = [0,dphi];
cbct.verbose = true;
cbct.GPU = 1;
cbct = Operators.ConeBeamScanner('multiHelix',nd,nd,[],zmax,rps,vtab,fps,nHelix,phaseShift);

% Compute forward transform of u0
disp('Computing multihelix scan');
f = cbct.apply(u0);


% Extract some scanner parameters 
nv = cbct.nv/nHelix;
na = cbct.na;
nb = cbct.nb;
SO = cbct.SO;
OD = cbct.OD;
scale = cbct.para.scale;
cos_phi = cbct.para.cos_phi;
sin_phi = cbct.para.sin_phi;
sd_z = scale*cbct.para.sd_z;
h = cbct.vtab/(2*pi*cbct.rps); 

% Compute detector positions
y_det = scale*(cbct.para.y_det);
z_det = scale*(cbct.para.z_det);

% Verify that f is correct

g1 = f.dataArrayNorm(:,:,:,1);
g2 = f.dataArrayNorm(:,:,:,2);

framelet = Transforms.FrameletSystem(3,'linear',2);
disp('Computing framelet expansion of f')
alpha1 = framelet.forwardTransform(g1);
alpha2 = framelet.forwardTransform(g2);

g1 = alpha1.frameletArray{2}{1,1};
g2 = alpha2.frameletArray{2}{1,1};

% Compute approximate John's equation 
disp('approximating Johns Equation')
R = SO+OD;
dzeta = -phaseShift(2)*h
dphi = 2*pi*cbct.rps/cbct.fps
da = scale*cbct.para.dy_det
db = scale*cbct.para.dz_det


% gtb
gtb = (g1(:,3:end,3:end)+g1(:,1:end-2,1:end-2)-g1(:,1:end-2,3:end)-g1(:,3:end,1:end-2))/(4*db*dphi);
% gaz
gaz = (g2(3:end,:,:)+g1(1:end-2,:,:)-g2(1:end-2,:,:)-g1(3:end,:,:))/(4*da*dzeta);
% gbz
gbz = (g2(:,3:end,:)+g1(:,1:end-2,:)-g2(:,1:end-2,:)-g1(:,3:end,:))/(4*db*dzeta);
% gb
gb = (g1(:,3:end,:)-g1(:,1:end-2,:))/(2*db);
% gbb
gbb = (g1(:,3:end,:)+g1(:,1:end-2,:)-2*g1(:,2:end-1,:))/(db*db);
% gab
gab = (g1(3:end,3:end,:)+g1(1:end-2,1:end-2,:)-g1(3:end,1:end-2,:)-g1(1:end-2,3:end,:))/(4*da*db);
% a
[a,b,NULL] = ndgrid(y_det,z_det,zeros(nv,1));
% b 


D = R*gtb(2:end-1,:,:) - R*SO*gaz(:,2:end-1,2:end-1) -...
    R*h*gbz(2:end-1,:,2:end-1) +...
    2.*a(2:end-1,2:end-1,2:end-1).*gb(2:end-1,:,2:end-1) +...
    a(2:end-1,2:end-1,2:end-1).*b(2:end-1,2:end-1,2:end-1).*gbb(2:end-1,:,2:end-1) +...
    (a(2:end-1,2:end-1,2:end-1).^2+R^2).*gab(:,:,2:end-1);







