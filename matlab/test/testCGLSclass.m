%Testing Fan Beam reconstruction with cgls
clc;clear all;close all;

lam = 1;
fbct = Operators.FanBeamScanner(256,64);
A = @(x)fbct.apply(x);
At = @(x)fbct.applyAdjoint(x);
u0 = DataTypes.ObjectData(2,single(zeros(256)),[10,10]);
cgiter = 500;
cgtol = 1e-14;

%% Noise free version
u = DataTypes.ObjectData(2,single(phantom(256)),[10,10]);
f0 = fbct.apply(u);
solver = Optimizers.CGLS(A,At,f0,u0,lam,cgiter,cgtol);

disp('Solving a fan beam reconstruction problem using CGLS - noise free');
disp('Scanner settings: ');
disp(fbct.para)
disp('Optimization settings: ');
disp(solver)

y = solver.solveCGLS();

figure('Name',['Noise free reconstruction']);
imshow(y.dataArray,[]);

%% Noisy version
sigma = 0.2;
u = DataTypes.ObjectData(2,single(phantom(256))+single(sigma*randn(256)),[10,10]);
f0 = fbct.apply(u);
solver = Optimizers.CGLS(A,At,f0,u0,lam,cgiter,cgtol);

disp('Solving a fan beam reconstruction problem using CGLS');
disp(['Noise level sigma = ',num2str(sigma)]);
disp('Scanner settings: ');
disp(fbct.para)
disp('Optimization settings: ');
disp(solver)

y = solver.solveCGLS();

figure('Name',['Noisy reconstruction with sigma = ',num2str(sigma)]);
imshow(y.dataArray,[]);

%% Test cone beam reconstruction using CGLS

lam = 10;
cbct = Operators.ConeBeamScanner(256,256,128);
A = @(x)cbct.apply(x);
At = @(x)cbct.applyAdjoint(x);
u3d = DataTypes.ObjectData(3,single(phantom3d(256)),[10,10,10]);
u03d = DataTypes.ObjectData(3,single(zeros(256,256,256)),[10,10,10]);
cgiter = 10;
cgtol = 1e-14;

f03d = A(u3d);

solver = Optimizers.CGLS(A,At,f03d,u03d,lam,cgiter,cgtol);

disp('Solving a noise free cone beam reconstruction using CGLS');
disp('Scanner settings: ');
disp(cbct.para);
disp('Optimization settings: ');
disp(solver);

y3d = solver.solveCGLS();


