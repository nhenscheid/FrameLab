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
#include "Atx_fan_mf.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//assuming dx=dy; during the computation, dx=dy=1; after the computation, the "dx=dy=scale" is taken into account.
//the output can be multiple 2D images, e.g., temporally-resolved images from 4DCT.
//the user needs to supply "id_Y" -- the projection index for each image.
{
    int *id_Y,*Nv,nx,ny,nt,nd,nv,tmp_size,version,GPU;
    float *x,*Y,*sd_phi,*y_det,*cos_phi,*sin_phi,*cos_det,*sin_det,SO,OD,scale,dy_det,y_os;

    Y=(float*)mxGetData(prhs[0]);

    version=(int)mxGetScalar(mxGetField(prhs[1],0,"version"));
    GPU=(int)mxGetScalar(mxGetField(prhs[1],0,"GPU"));
    SO=(float)mxGetScalar(mxGetField(prhs[1],0,"SO"));
    OD=(float)mxGetScalar(mxGetField(prhs[1],0,"OD"));
    scale=(float)mxGetScalar(mxGetField(prhs[1],0,"scale"));
    nx=(int)mxGetScalar(mxGetField(prhs[1],0,"nx"));
    ny=(int)mxGetScalar(mxGetField(prhs[1],0,"ny"));
    sd_phi=(float*)mxGetData(mxGetField(prhs[1],0,"sd_phi"));
    y_det=(float*)mxGetData(mxGetField(prhs[1],0,"y_det"));
    nd=(int)mxGetNumberOfElements(mxGetField(prhs[1],0,"y_det"));
    dy_det=(float)mxGetScalar(mxGetField(prhs[1],0,"dy_det"));
    nt=(int)mxGetScalar(mxGetField(prhs[1],0,"nt"));
    id_Y=(int*)mxGetData(mxGetField(prhs[1],0,"id_Y"));
    Nv=(int*)mxGetData(mxGetField(prhs[1],0,"Nv"));
    tmp_size=(int)mxGetScalar(mxGetField(prhs[1],0,"tmp_size"));
    nv=(int)mxGetScalar(mxGetField(prhs[1],0,"nv"));
    cos_phi=(float*)mxGetData(mxGetField(prhs[1],0,"cos_phi"));
    sin_phi=(float*)mxGetData(mxGetField(prhs[1],0,"sin_phi"));

    plhs[0]=mxCreateNumericMatrix(nx*ny*nt,1,mxSINGLE_CLASS,mxREAL);
    x=(float*)mxGetData(plhs[0]);

    switch(version)
    {   case 1:
            cos_det=(float*)mxGetData(mxGetField(prhs[1],0,"cos_det"));
            sin_det=(float*)mxGetData(mxGetField(prhs[1],0,"sin_det"));
            y_os=(float)mxGetScalar(mxGetField(prhs[1],0,"y_os"));
            if(GPU==0)
            {Atx_fan_mf_cpu_new(x,Y,SO,OD,scale,nx,ny,sd_phi,nd,y_det,dy_det,y_os,nt,id_Y,Nv,tmp_size,nv,cos_phi,sin_phi,cos_det,sin_det);}
            else
            {Atx_fan_mf_gpu_new(x,Y,SO,OD,scale,nx,ny,sd_phi,nd,y_det,dy_det,y_os,nt,id_Y,Nv,tmp_size,nv,cos_phi,sin_phi,cos_det,sin_det);}
        break;
        case 2:
            cos_det=(float*)mxGetData(mxGetField(prhs[1],0,"cos_det2"));
            sin_det=(float*)mxGetData(mxGetField(prhs[1],0,"sin_det2"));
            if(GPU==0)
            {Atx_fan_mf_cpu_new_fb(x,Y,SO,OD,scale,nx,ny,sd_phi,nd,y_det,dy_det,nt,id_Y,Nv,tmp_size,nv,cos_phi,sin_phi,cos_det,sin_det);}
            else
            {Atx_fan_mf_gpu_new_fb(x,Y,SO,OD,scale,nx,ny,sd_phi,nd,y_det,dy_det,nt,id_Y,Nv,tmp_size,nv,cos_phi,sin_phi,cos_det,sin_det);}
        break;
    }
}

