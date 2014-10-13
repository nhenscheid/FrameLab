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

#include "mex.h"
#include "Ax_cone_mf.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//assuming dx=dy; during the computation, dx=dy=1; after the computation, the "dx=dy=scale" is taken into account.
//the input can be multiple 3D images, e.g., temporally-resolved images from 4DCT.
//the user needs to supply
//      "id_X" -- the image index for each projection
//      "id_Y" -- the projection index for each image
//      "nv_block" -- to deal with the GPU time limit
{
    int *id_X,*id_Y,*Nv,nx,ny,nz,nt,nv,na,nb,tmp_size,nv_block,version,GPU;
    float *X,*y,*sd_phi,*sd_z,*y_det,*z_det,SO,OD,scale,dz;

    X=(float*)mxGetData(prhs[0]);

    version=(int)mxGetScalar(mxGetField(prhs[1],0,"version"));
    GPU=(int)mxGetScalar(mxGetField(prhs[1],0,"GPU"));
    SO=(float)mxGetScalar(mxGetField(prhs[1],0,"SO"));
    OD=(float)mxGetScalar(mxGetField(prhs[1],0,"OD"));
    scale=(float)mxGetScalar(mxGetField(prhs[1],0,"scale"));
    nx=(int)mxGetScalar(mxGetField(prhs[1],0,"nx"));
    ny=(int)mxGetScalar(mxGetField(prhs[1],0,"ny"));
    nz=(int)mxGetScalar(mxGetField(prhs[1],0,"nz"));
    nt=(int)mxGetScalar(mxGetField(prhs[1],0,"nt"));
    sd_phi=(float*)mxGetData(mxGetField(prhs[1],0,"sd_phi"));
    sd_z=(float*)mxGetData(mxGetField(prhs[1],0,"sd_z"));
    nv=(int)mxGetNumberOfElements(mxGetField(prhs[1],0,"sd_phi"));
    y_det=(float*)mxGetData(mxGetField(prhs[1],0,"y_det"));
    na=(int)mxGetNumberOfElements(mxGetField(prhs[1],0,"y_det"));
    z_det=(float*)mxGetData(mxGetField(prhs[1],0,"z_det"));
    nb=(int)mxGetNumberOfElements(mxGetField(prhs[1],0,"z_det"));
    dz=(float)mxGetScalar(mxGetField(prhs[1],0,"dz"));
    id_Y=(int*)mxGetData(mxGetField(prhs[1],0,"id_Y"));
    Nv=(int*)mxGetData(mxGetField(prhs[1],0,"Nv"));
    tmp_size=(int)mxGetScalar(mxGetField(prhs[1],0,"tmp_size"));
    nv_block=(int)mxGetScalar(mxGetField(prhs[1],0,"nv_block"));
    id_X=(int*)mxGetData(mxGetField(prhs[1],0,"id_X"));

    plhs[0]=mxCreateNumericMatrix(nv*na*nb,1,mxSINGLE_CLASS,mxREAL);
    y=(float*)mxGetData(plhs[0]);

    switch(version)
    {   case 0:
            if(GPU==0)
            {Ax_cone_mf_cpu_siddon(X,y,SO,OD,scale,nx,ny,nz,nv,sd_phi,sd_z,na,nb,y_det,z_det,id_X,dz);}
            else
            {Ax_cone_mf_gpu_siddon(X,y,sd_phi,sd_z,y_det,z_det,id_Y,Nv,SO,OD,scale,dz,nx,ny,nz,nt,na,nb,nv,tmp_size,nv_block);}
        break;
        case 1:
            if(GPU==0)
            {Ax_cone_mf_cpu_new(X,y,SO,OD,scale,nx,ny,nz,nv,sd_phi,sd_z,na,nb,y_det,z_det,id_X,dz);}
            else
            {Ax_cone_mf_gpu_new(X,y,sd_phi,sd_z,y_det,z_det,id_Y,Nv,SO,OD,scale,dz,nx,ny,nz,nt,na,nb,nv,tmp_size,nv_block);}
        break;
    }
}

