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

% cone-beam geometry %
nx=256;ny=nx;nz=256;
dx=single(250/nx);dy=dx;dz=single(250/nz);
% nv=4;
nv=668;
SO=single(1000);OD=single(500);
na=512;dy_det=single(0.388*1024/na);
nb=384;dz_det=dy_det;
sd_phi=single(2*pi/nv*(0:nv-1));    
sd_z=single(zeros(1,nv)); % for helical scan
y_os=single(0*dy_det);
y_det=single(((-na/2:na/2-1)+0.5)*dy_det+y_os);
z_det=single(((-nb/2:nb/2-1)+0.5)*dz_det);
% normalize so that dx=dy=1 %
scale=dx;
sd_z=sd_z/scale;SO=SO/scale;OD=OD/scale;
y_det=y_det/scale;z_det=z_det/scale;
dy_det=dy_det/scale;dz_det=dz_det/scale;
y_os=y_os/scale;dz=dz/scale;
% parameters (single image frame) %
% note: this version works for multiple images as well by supplying id_v -- which projection views are for each frame%
%load phantom_3d.mat x0 % 3D phantom image
x0 = phantom3d(nx);
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
    'SO',single(SO),'OD',single(OD),'scale',single(scale),'nx',uint32(nx),'ny',uint32(ny),'nz',uint32(nz),'nt',uint32(nt),...
    'sd_phi',sd_phi,'sd_z',sd_z,'y_det',y_det,'z_det',z_det,'cos_phi',single(cos(sd_phi)),'sin_phi',single(sin(sd_phi)),...
    'dy_det',single(dy_det),'y_os',single(y_os/dy_det),'dz_det',single(dz_det),'dz',single(dz),...
    'id_X',uint32(id_X),'id_Y',uint32(id_Y),'Nv',uint32(Nv),'tmp_size',uint32(tmp_size),'nv_block',uint32(4));

% X-ray Transform %
% para.version=uint32(0);para.GPU=uint32(0); % Siddon's algorithm - CPU
% tic;y_siddon_cpu=Ax_cone_mf(X0,para);toc;
para.version=uint32(0);para.GPU=uint32(1); % Siddon's algorithm - GPU
tic;y_siddon_gpu=Ax_cone_mf(X0,para);toc;
% para.version=uint32(1);para.GPU=uint32(0); % new algorithm - CPU
% tic;y_new_cpu=Ax_cone_mf(X0,para);toc;
para.version=uint32(1);para.GPU=uint32(1); % new algorithm - GPU
tic;y_new_gpu=Ax_cone_mf(X0,para);toc;

% [mean(abs(y_siddon_gpu-y_siddon_cpu)) mean(abs(y_new_cpu-y_siddon_cpu)) mean(abs(y_new_gpu-y_siddon_cpu))]/mean(abs(y_siddon_cpu))
% y_siddon_cpu=reshape(y_siddon_cpu,[na nb nv]);
% y_siddon_gpu=reshape(y_siddon_gpu,[na nb nv]);
% y_new_cpu=reshape(y_new_cpu,[na nb nv]);
% y_new_gpu=reshape(y_new_gpu,[na nb nv]);
% iv=nv/2;figure;imshow([y_siddon_cpu(:,:,iv) abs(y_siddon_gpu(:,:,iv)-y_siddon_cpu(:,:,iv))...
%     abs(y_new_cpu(:,:,iv)-y_siddon_cpu(:,:,iv)) abs(y_new_gpu(:,:,iv)-y_siddon_cpu(:,:,iv))],[]);

% Adjoint X-ray Transform %
y=y_new_gpu;
% para.version=uint32(1);para.GPU=uint32(0); % new algorithm (adjoint X-ray transform) - CPU
% tic;x_new_cpu=Atx_cone_mf(y,para);toc;
para.version=uint32(1);para.GPU=uint32(1); % new algorithm (adjoint X-ray transform) - GPU
tic;x_new_gpu=Atx_cone_mf(y,para);toc;
return
[mean(abs(x_new_gpu-x_new_cpu))]/mean(abs(x_new_cpu))
x_new_cpu=reshape(x_new_cpu,[nx ny nz]);
x_new_gpu=reshape(x_new_gpu,[nx ny nz]);
iz=nz/2;figure;imshow([x_new_cpu(:,:,iz) abs(x_new_gpu(:,:,iz)-x_new_cpu(:,:,iz))],[]);
