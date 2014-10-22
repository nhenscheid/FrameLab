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
#include "Ax_fan_mf.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//assuming dx=dy; during the computation, dx=dy=1; after the computation, the "dx=dy=scale" is taken into account.
//the input can be multiple 2D images, e.g., temporally-resolved images from 4DCT.
//the user needs to supply "id_X" -- the image index for each projection.
{
    int *id_X,nx,ny,nv,nt,nd,version,GPU;
    float *X,*y,*sd_phi,*y_det,SO,OD,scale,dy_det;

    X=(float*)mxGetData(prhs[0]);

    version=(int)mxGetScalar(mxGetField(prhs[1],0,"version"));
    GPU=(int)mxGetScalar(mxGetField(prhs[1],0,"GPU"));
    SO=(float)mxGetScalar(mxGetField(prhs[1],0,"SO"));
    OD=(float)mxGetScalar(mxGetField(prhs[1],0,"OD"));
    scale=(float)mxGetScalar(mxGetField(prhs[1],0,"scale"));
    nx=(int)mxGetScalar(mxGetField(prhs[1],0,"nx"));
    ny=(int)mxGetScalar(mxGetField(prhs[1],0,"ny"));
    sd_phi=(float*)mxGetData(mxGetField(prhs[1],0,"sd_phi"));
    nv=(int)mxGetNumberOfElements(mxGetField(prhs[1],0,"sd_phi"));
    y_det=(float*)mxGetData(mxGetField(prhs[1],0,"y_det"));
    nd=(int)mxGetNumberOfElements(mxGetField(prhs[1],0,"y_det"));
    id_X=(int*)mxGetData(mxGetField(prhs[1],0,"id_X"));
    nt=(int)mxGetScalar(mxGetField(prhs[1],0,"nt"));
    dy_det=(float)mxGetScalar(mxGetField(prhs[1],0,"dy_det"));

    plhs[0]=mxCreateNumericMatrix(nv*nd,1,mxSINGLE_CLASS,mxREAL);
    y=(float*)mxGetData(plhs[0]);
    
    Ax_fan_mf_cpu_new(X,y,SO,OD,scale,nx,ny,nv,sd_phi,nd,y_det,id_X);
}

