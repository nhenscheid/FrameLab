function [u,f0,err]=testTWFRecon(type,Niter,lam,mu,sigma)

%2D and 3D Fan/Cone beam CT inversion via framelet analysis model
% Usage: 
% [u,f0,err] = testTWFRecon(type,Niter,lam,mu,sigma)
%
% Parameters:
% type is either 'fan', 'circle' or 'helix'
% lam & mu are augmented Lagrangian parameters - mu is the 'Tikhonov'
% parameter and lam is the 'sparsity' parameter.
% sigma is the standard deviation of the additive Gaussian noise
% 
% Details:
% Optimization is via the augmented Lagrangian method (e.g. Nocedal &
% Wright 2006, Dong and Zhang 2012)
% Nick Henscheid, 2013
% Subroutines by:
% Jerome Gilles, (UCLA, Bregman Cookbook),
% Jian-Feng Cai (NUS, Framelet decompositions)
% Hao Gao (SJTU, Fast CT transforms) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
%*** NOTES ***%
% lam=0.001 and mu=0.1 seem to work best. 

% Define forward and adjoint cone beam operators
%cbct = Operators.ConeBeamScanner('circle',64,64,64);
nd = 128;
rps = 1;
fps = 32;
zmax = 1;
vtab = 1;
%nHelix = 3;
%dphi = 2*pi*rps/fps;
%phaseShift = dphi*(0:nHelix-1);
%cbct = Operators.ConeBeamScanner('multiHelix',nd,nd,[],zmax,rps,vtab,fps,nHelix,phaseShift);
cbct = Operators.ConeBeamScanner('helix',nd,nd,[],zmax,rps,vtab,fps);
A = @(x)cbct.apply(x);
At = @(x)cbct.applyAdjoint(x);

% Generate phantom, projection and initial guess
NTrue = 128;
%load MRbrain.mat;
%utrue=single(I)./max(max(max(I)));clear I;
uTrue = DataTypes.ObjectData(3,single(phantom3d(NTrue)),[2,2,2]);
NRecon = 128;
u0 = DataTypes.ObjectData(3,single(zeros(NRecon,NRecon,NRecon)),[2,2,2]);
f0 = A(uTrue)
% Add Gaussian noise with deviation sigma to data
%NOISE CURRENTLY DISABLED
%f0 = f0+sigma*randn(cbct.na,cbct.nb,cbct.nv);
Atf0 = At(f0);  % Backprojected data

% Create framelet transforms 
framelet = Transforms.FrameletSystem(3,'linear',1);
W = @(x)(x.frameletTransform(framelet));
Wt = @(x)(DataTypes.ObjectData(3,x.adjointFrameletTransform(framelet),[2,2,2]));
alpha = W(u0);
v = W(u0);  % v is the Lagrange multiplier

% Show the error interactively?
plots=1;

% Conjugate gradient variables

k=0;

% Prepare figures for updating
if (plots==1)
[fig3,h1,h2,h3,h4]=makeplots();
end
% Vector for relative errors from utrue 
err=zeros(Niter,1);

tic

% Create CGLS class to solve least squares
cgiter=50;
cgtol=1e-12;
b = Atf0+mu*Wt(alpha-v);
solver = Optimizers.CGLS(A,At,b,u0,mu,cgiter,cgtol);

while k<Niter     
    k=k+1;
    % Step 1: Approximate u(k+1) by solving (PtP+mu*I)u=Ptf0+muWt(alpha-v) 
    % using conjugate gradient.  
    u = solver.solveCGLS();
    Wu = W(u);
    % Step 2: find alpha(k+1) by performing framelet soft thresholding
    alpha = v + Wu;
    alpha.softThreshold(lam/mu);
    
    % Step 3: update Lagrange multiplier
    v = v + mu*(Wu-alpha);
    %v = v + (Wu-alpha);
    % Step 4: update CGLS solver
    b = Atf0 + mu*Wt(alpha-v);
    clear solver;
    solver = Optimizers.CGLS(A,At,b,u,mu,cgiter,cgtol);
    % Plot the image, display the error
    err(k)=norm(u-uTrue)/norm(uTrue);
    if (plots==1)
    set(h1,'String',num2str(k));
    set(h3,'String',num2str(err(k)));
    drawnow
    end
  
end
   

end

function [fig3,h1,h2,h3,h4]=makeplots()
fig3=figure('Position',[0 500 200 80],'ToolBar','none','MenuBar','none');
uicontrol('Style','text','Position',[0 60 100 20],'String','Iteration');
uicontrol('Style','text','Position',[100 60 100 20],'String','CG Iter');
uicontrol('Style','text','Position',[0 20 100 20],'String','Tot Error');
uicontrol('Style','text','Position',[100 20 100 20],'String','CG Error');
h1=uicontrol('Style','text','Position',[0 40 100 20]);
h2=uicontrol('Style','text','Position',[100 40 100 20]);
h3=uicontrol('Style','text','Position',[0 0 100 20]);
h4=uicontrol('Style','text','Position',[100 0 100 20]);
end