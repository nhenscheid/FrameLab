// *----------------------------------------------
// 	Author Contact Information:
// 	Hao Gao
// 	hao.gao@emory.edu || hao.gao.2012@gmail.com
// 	Department of Mathematics and Computer Science, Emory University
// 	Department of Radiology and Imaging Sciences, Emory University
//
// 	Copyright (c) Hao Gao 2012
// ----------------------------------------------*/
//
// If you find this code useful, you may cite the following reference:
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
// The full source codes are available at https://sites.google.com/site/fastxraytransform

#include <stdio.h>
#include "mex.h"
#include "Atx_cone_mf.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//assuming dx=dy; during the computation, dx=dy=1; after the computation, the "dx=dy=scale" is taken into account.
//the output can be multiple 3D images, e.g., temporally-resolved images from 4DCT.
//the user needs to supply
//      "id_Y" -- the projection index for each image
//      "nv_block" -- to deal with the GPU time limit
{
    int *id_Y,*Nv,nx,ny,nz,nt,nv,na,nb,tmp_size,nv_block,version,GPU;
    float *X,*y,*sd_z,*y_det,*z_det,*cos_phi,*sin_phi,SO,OD,scale,dy_det,y_os,dz_det,dz;

    y=(float*)mxGetData(prhs[0]);

    version=(int)mxGetScalar(mxGetField(prhs[1],0,"version"));
    GPU=(int)mxGetScalar(mxGetField(prhs[1],0,"GPU"));
    SO=(float)mxGetScalar(mxGetField(prhs[1],0,"SO"));
    OD=(float)mxGetScalar(mxGetField(prhs[1],0,"OD"));
    scale=(float)mxGetScalar(mxGetField(prhs[1],0,"scale"));
    nx=(int)mxGetScalar(mxGetField(prhs[1],0,"nx"));
    ny=(int)mxGetScalar(mxGetField(prhs[1],0,"ny"));
    nz=(int)mxGetScalar(mxGetField(prhs[1],0,"nz"));
    nt=(int)mxGetScalar(mxGetField(prhs[1],0,"nt"));
    cos_phi=(float*)mxGetData(mxGetField(prhs[1],0,"cos_phi"));
    sin_phi=(float*)mxGetData(mxGetField(prhs[1],0,"sin_phi"));
    nv=(int)mxGetNumberOfElements(mxGetField(prhs[1],0,"sin_phi"));
    sd_z=(float*)mxGetData(mxGetField(prhs[1],0,"sd_z"));
    y_det=(float*)mxGetData(mxGetField(prhs[1],0,"y_det"));//initial y position of the detector
    na=(int)mxGetNumberOfElements(mxGetField(prhs[1],0,"y_det"));
    z_det=(float*)mxGetData(mxGetField(prhs[1],0,"z_det"));//initial z position of the detector
    nb=(int)mxGetNumberOfElements(mxGetField(prhs[1],0,"z_det"));
    y_os=(float)mxGetScalar(mxGetField(prhs[1],0,"y_os"));
    dy_det=(float)mxGetScalar(mxGetField(prhs[1],0,"dy_det"));
    dz_det=(float)mxGetScalar(mxGetField(prhs[1],0,"dz_det"));
    dz=(float)mxGetScalar(mxGetField(prhs[1],0,"dz"));
    id_Y=(int*)mxGetData(mxGetField(prhs[1],0,"id_Y"));
    Nv=(int*)mxGetData(mxGetField(prhs[1],0,"Nv"));
    tmp_size=(int)mxGetScalar(mxGetField(prhs[1],0,"tmp_size"));
    nv_block=(int)mxGetScalar(mxGetField(prhs[1],0,"nv_block"));

    plhs[0]=mxCreateNumericMatrix(nx*ny*nz*nt,1,mxSINGLE_CLASS,mxREAL);
    X=(float*)mxGetData(plhs[0]);

    switch(version)
    {   case 1:
            if(GPU==0)
            {Atx_cone_mf_cpu_new(X,y,SO,OD,scale,nx,ny,nz,sd_z,y_det,z_det,dy_det,y_os,dz_det,dz,nt,id_Y,Nv,tmp_size,nv,cos_phi,sin_phi,na,nb);}
            else
            {Atx_cone_mf_gpu_new(X,y,cos_phi,sin_phi,sd_z,y_det,z_det,id_Y,Nv,SO,OD,scale,dy_det,y_os,dz_det,dz,nx,ny,nz,nt,na,nb,nv,tmp_size,nv_block);}
        break;
    }
}

