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