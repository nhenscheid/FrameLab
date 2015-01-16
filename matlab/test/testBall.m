% Test the fan and cone beam scanners with a unit disc/ball
%% Fan beam geometry
N = 1000;
u0 = DataTypes.ObjectData(2,single(unitBall(N,2)),[2,2]);
nd = 10;
nv = 1;
fbct = Operators.FanBeamScanner(nd,nv)
y = fbct.apply(u0);
y = y.dataArray;
y_det = fbct.para.y_det*fbct.para.scale;

disp('Computing an exact solution.  Should make this a function.');
l = @(a1)4*sqrt(-6+625/(100+a1^2));

yExact = zeros(nd,1);
for i=1:nd
    if imag(l(y_det(i)))==0
        yExact(i) = l(y_det(i));
    end
end
disp(sprintf('%s%d\n','Max error is ',max(abs(y-yExact(:)))/max(abs(yExact))))
disp(sprintf('%s%d\n','L2 error is ',norm(y-yExact)/norm(yExact)))
%% Circular source path
clear all;
N = 256;
u0 = DataTypes.ObjectData(3,single(unitBall(N,3)),[2,2,2]); %unit ball, diameter = 2
nd = 128;
nv = 10;

disp('Testing circular scan');
cbct = Operators.ConeBeamScanner('circle',nd,nd,nv);

cbct.verbose = true;
cbct.GPU = 1;

y = cbct.apply(u0);

y_det = cbct.para.y_det*cbct.para.scale;
z_det = cbct.para.z_det*cbct.para.scale;

% Exact solution
[Y,Z] = ndgrid(y_det,z_det);
yExact = coneBeamSphere(0,Y,Z,0,0,cbct.SO,cbct.OD,0);

for i = 1:nv
    disp('For view ');
    i
    y1 = y.dataArray(:,:,i);
    disp(sprintf('%s%d\n','Max error is ',max(abs(yExact(:)-y1(:)))))
    disp(sprintf('%s%d\n','L2 error is ',sqrt(sum(abs(yExact(:)-y1(:)).^2))/sqrt(sum(abs(yExact(:))))))
end

%% Helix 
clear all;
N = 256;
u0 = DataTypes.ObjectData(3,single(unitBall(N,3)),[2,2,2]); %unit ball, diameter = 2
nd = 256;
nv = 32;
disp('Testing helical scan');
cbct = Operators.ConeBeamScanner('helix',nd,nd,[],2,4,2,nv);
cbct.verbose = true;
cbct.GPU = 1;

y = cbct.apply(u0);

disp('Computing exact solution')
y_det = cbct.para.y_det*cbct.para.scale;
z_det = cbct.para.z_det*cbct.para.scale;
sd_phi = cbct.para.sd_phi;
h = cbct.vtab/(2*pi*cbct.rps); 

[A,B,T] = ndgrid(y_det,z_det,sd_phi);
yExact = coneBeamSphere(T,A,B,0,h,cbct.SO,cbct.OD,0);
yExactNorm = coneBeamSphere(T,A,B,0,h,cbct.SO,cbct.OD,1);

disp(sprintf('%s%d\n','L2 error for unnormalized transform is ',sqrt(sum(abs(yExact(:)-y.dataArray(:)).^2))/sqrt(sum(abs(yExact(:))))))

disp(sprintf('%s%d\n','L2 error for normalized transform is ',sqrt(sum(abs(yExactNorm(:)-y.dataArrayNorm(:)).^2))/sqrt(sum(abs(yExactNorm(:))))))

%% Multihelix
clear all;
N = 256;
u0 = DataTypes.ObjectData(3,single(unitBall(N,3)),[2,2,2]); %unit ball, diameter = 2
nd = 128;
rps = 4;
fps = 64;
zmax = 2;
vtab = 2;
nHelix = 2;
dphi = 2*pi*rps/fps;
phaseShift = [0,dphi];
disp('Testing multihelix scan');
cbct = Operators.ConeBeamScanner('multiHelix',nd,nd,[],zmax,rps,vtab,fps,nHelix,phaseShift);

cbct.verbose = true;
cbct.GPU = 1;
y = cbct.apply(u0);

y_det = cbct.para.y_det*cbct.para.scale;
z_det = cbct.para.z_det*cbct.para.scale;
sd_phi = cbct.para.sd_phi;
h = cbct.vtab/(2*pi*cbct.rps); 
nv = cbct.nv/cbct.nHelix;

for i=1:nHelix
    disp(sprintf('%s%i\n','For helix # ',i))
    [A,B,T] = ndgrid(y_det,z_det,sd_phi((i-1)*nv+1:i*nv));
    zeta = -phaseShift(i)*h;
    yExact = coneBeamSphere(T,A,B,zeta,h,cbct.SO,cbct.OD,0);
    yExactNorm = coneBeamSphere(T,A,B,zeta,h,cbct.SO,cbct.OD,1);
    yApprox = y.dataArray(:,:,:,i);
    yApproxNorm = y.dataArrayNorm(:,:,:,i);
    disp(sprintf('%s%d\n','L2 error for unnormalized transform is ',sqrt(sum(abs(yExact(:)-yApprox(:)).^2))/sqrt(sum(abs(yExact(:))))))
    disp(sprintf('%s%d\n','L2 error for normalized transform is ',sqrt(sum(abs(yExactNorm(:)-yApproxNorm(:)).^2))/sqrt(sum(abs(yExactNorm(:))))))
end




