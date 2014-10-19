% *----------------------------------------------
% 	Author Contact Information:
% 	Hao Gao
% 	hao.gao@emory.edu || hao.gao.2012@gmail.com
% 	Department of Mathematics and Computer Science, Emory University
% 	Department of Radiology and Imaging Sciences, Emory University
% 
% 	Copyright (c) Hao Gao 2012
% ----------------------------------------------*/
%
% If you find this code useful, you may cite the following reference:
% H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
% The full source codes are available at https://sites.google.com/site/fastxraytransform

clc;clear all;close all;

% fan-beam geometry %
nx=256;ny=nx;% number of pixels: nx*ny
dx=single(250/nx);dy=dx;% pixel size: dx*dy
nv=668;% number of views (projections)
SO=single(1000);% the distance from source to isocenter
OD=single(500);% the distance from isocenter to detector
nd=512;% number of detector
dy_det=single(0.388*1024/nd);% size of detector
sd_phi=single(2*pi/nv*(0:nv-1));% view angles
y_os=single(0*dy_det);% isocenter offset with respect to the detector center
y_det=single(((-nd/2:nd/2-1)+0.5)*dy_det)+y_os;% detector coordinate
y_det2=single((-nd/2:nd/2)*dy_det)+y_os;
% normalize so that dx=dy=1 %
scale=dx;
SO=SO/scale;OD=OD/scale;y_det=y_det/scale;
dy_det=dy_det/scale;y_os=y_os/scale;
% parameters (single image frame) %
% note: this version works for multiple images as well by supplying id_v -- which projection views are for each frame%
%load phantom_2d.mat x0 % 2D phantom image
x0=single(phantom(256));
nt=1;
X0=x0(:);
Id_v=cell(1,nt);
Id_v{1}=uint32(0:nv-1);
id_X=[];Nv=zeros(1,nt);
for i=1:nt
    id_X=[id_X (i-1)*ones(1,numel(Id_v{i}))];
    Nv(i)=numel(Id_v{i});
end
tmp_size=max(Nv);
id_Y=zeros(tmp_size,nt);
for i=1:nt
    id_Y(1:Nv(i),i)=Id_v{i};
end
para=struct('version',[],'GPU',[],...
    'SO',single(SO),'OD',single(OD),'scale',single(scale),'nx',uint32(nx),'ny',uint32(ny),'nv',uint32(nv),...
    'sd_phi',sd_phi,'y_det',y_det,'id_X',uint32(id_X),'nt',uint32(nt),'dy_det',dy_det,'y_os',y_os,...
    'id_Y',uint32(id_Y),'Nv',uint32(Nv),'tmp_size',uint32(tmp_size),...
    'cos_phi',cos(sd_phi),'sin_phi',sin(sd_phi),'cos_det',[],'sin_det',[]);        
angle_det=atan2(y_det,SO+OD);para.cos_det=cos(angle_det);para.sin_det=sin(angle_det);
angle_det2=atan2(y_det2,SO+OD);para.cos_det2=single(cos(angle_det2));para.sin_det2=single(sin(angle_det2));% for finite-size beam

% X-ray Transform %
para.version=uint32(0);para.GPU=uint32(0); % Siddon's algorithm - CPU
tic;y_siddon_cpu=Ax_fan_mf(X0,para);toc;
para.version=uint32(0);para.GPU=uint32(1); % Siddon's algorithm - GPU
tic;y_siddon_gpu=Ax_fan_mf(X0,para);toc;
para.version=uint32(1);para.GPU=uint32(0); % new algorithm - CPU
tic;y_new_cpu=Ax_fan_mf(X0,para);toc;
para.version=uint32(1);para.GPU=uint32(1); % new algorithm - GPU
tic;y_new_gpu=Ax_fan_mf(X0,para);toc;
para.version=uint32(2);para.GPU=uint32(0); % new algorithm (finite beam) - CPU
tic;y_new_cpu_fb=Ax_fan_mf(X0,para);toc;
para.version=uint32(2);para.GPU=uint32(1); % new algorithm (finite beam) - GPU
tic;y_new_gpu_fb=Ax_fan_mf(X0,para);toc;

y_siddon_cpu=reshape(y_siddon_cpu,[nd nv]);
y_siddon_gpu=reshape(y_siddon_gpu,[nd nv]);
y_new_cpu=reshape(y_new_cpu,[nd nv]);
y_new_gpu=reshape(y_new_gpu,[nd nv]);
y_new_cpu_fb=reshape(y_new_cpu_fb,[nd nv]);
y_new_gpu_fb=reshape(y_new_gpu_fb,[nd nv]);
figure;imshow([y_siddon_cpu abs(y_siddon_gpu-y_siddon_cpu) abs(y_new_cpu-y_siddon_cpu) abs(y_new_gpu-y_siddon_cpu)],[]);
figure;imshow([y_new_cpu_fb abs(y_new_gpu_fb-y_new_cpu_fb)],[]);

% Adjoint X-ray Transform %
y=y_new_gpu;
para.version=uint32(1);para.GPU=uint32(0); % new algorithm (adjoint X-ray transform) - CPU
tic;x_new_cpu=Atx_fan_mf(y,para);toc;
para.version=uint32(1);para.GPU=uint32(1); % new algorithm (adjoint X-ray transform) - GPU
tic;x_new_gpu=Atx_fan_mf(y,para);toc;
para.version=uint32(2);para.GPU=uint32(0); % new algorithm (adjoint X-ray transform - finite beam) - CPU
tic;x_new_cpu_fb=Atx_fan_mf(y,para);toc;
para.version=uint32(2);para.GPU=uint32(1); % new algorithm (adjoint X-ray transform - finite beam) - GPU
tic;x_new_gpu_fb=Atx_fan_mf(y,para);toc;

x_new_cpu=reshape(x_new_cpu,[nx ny]);
x_new_gpu=reshape(x_new_gpu,[nx ny]);
x_new_cpu_fb=reshape(x_new_cpu_fb,[nx ny]);
x_new_gpu_fb=reshape(x_new_gpu_fb,[nx ny]);
figure;imshow([x_new_cpu abs(x_new_gpu-x_new_cpu)],[]);
figure;imshow([x_new_cpu_fb abs(x_new_gpu_fb-x_new_cpu_fb)],[]);
