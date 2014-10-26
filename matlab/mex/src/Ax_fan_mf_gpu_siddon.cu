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

#include <math.h>
#include <stdlib.h>
#include "sort_alpha_d.h"

#define BLOCK_SIZE_x 16
#define BLOCK_SIZE_y 16
#define Nax 512+1 // Please note that you need increase the value of Nax manually, if Nx>512, since the thread doesn't support malloc.
#define Nay 512+1 // Please note that you need increase the value of Nay manually, if Nx>512, since the thread doesn't support malloc.
#define Na 512+512+2 // Please note that you need increase the value of Na manually, if Nx+Ny>1024, since the thread doesn't support malloc.

extern "C" void Ax_fan_mf_gpu_siddon(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X,int nt);

__global__ void Ax_fan_mf_gpu_kernel_siddon(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx0=threadIdx.x;
	int ty0=threadIdx.y;

	int iv=bx*BLOCK_SIZE_x+tx0;
	int id=by*BLOCK_SIZE_y+ty0;

    if(iv<nv&&id<nd)
    {
		int Ji,n,nx2,ny2,iy,i,j,n_alpha,tx,ty;
		float *x,alpha_x[Nax],alpha_y[Nay],alpha[Na];
		float Si,cos_phi,sin_phi,x1,y1,x2,y2,dx,dy,totallength,alpha_c;

	    iy=iv*nd+id;
		nx2=nx/2;ny2=ny/2;
		n=nx*ny;

        x=&X[id_X[iv]*n];
        cos_phi=cosf(sd_phi[iv]);sin_phi=sinf(sd_phi[iv]);
        x1=cos_phi*(-SO);
        y1=sin_phi*(-SO);
		x2=cos_phi*OD-sin_phi*y_det[id];
        y2=sin_phi*OD+cos_phi*y_det[id];
        // step 1: compute alpha_x and alpha_y
        dx=x2-x1;dy=y2-y1;
        if(x1<x2)
        {   for(i=-nx2,j=0;i<=nx2;i++,j++)
            {alpha_x[j]=(i-x1)/dx;}
        }
        else if(x1>x2)
        {   for(i=nx2,j=0;i>=-nx2;i--,j++)
            {alpha_x[j]=(i-x1)/dx;}
        }

        if(y1<y2)
        {   for(i=-ny2,j=0;i<=ny2;i++,j++)
            {alpha_y[j]=(i-y1)/dy;}
        }
        else if(y1>y2)
        {   for(i=ny2,j=0;i>=-ny2;i--,j++)
            {alpha_y[j]=(i-y1)/dy;}
        }
        // step 2: sort and bin alpha_x, alpha_y into alpha
        if(y1==y2)
        {   for(i=0;i<=nx;i++)
            {alpha[i]=alpha_x[i];}
            n_alpha=nx+1;
        }
        else
        {   if(x1==x2)
            {   for(i=0;i<=ny;i++)
                {alpha[i]=alpha_y[i];}
                n_alpha=ny+1;
            }
            else
            {n_alpha=sort_alpha(nx+1,ny+1,alpha_x,alpha_y,alpha);}
        }
        // step 3: compute x-ray transform y[iy]
        totallength=sqrtf(dx*dx+dy*dy);
        y[iy]=0;
        for(i=0;i<n_alpha-1;i++)
        {   alpha_c=(float)0.5*(alpha[i]+alpha[i+1]);
            tx=(int)floor(x1+alpha_c*dx+nx2);
            if(tx>=0&&tx<nx)
            {   ty=(int)floorf(y1+alpha_c*dy+ny2);
                if(ty>=0&&ty<ny)
                {   Ji=ty*nx+tx;
                    Si=(alpha[i+1]-alpha[i])*totallength;
                    y[iy]+=Si*x[Ji];
                }
            }
        }
        y[iy]*=scale;
    }
}

void Ax_fan_mf_gpu_siddon(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X,int nt)
// Siddon's method to compute the X-ray transform
{   float *y_d,*X_d,*sd_phi_d,*y_det_d;
	int *id_X_d;

	cudaMalloc(&y_d,nv*nd*sizeof(float));
	cudaMalloc(&X_d,nx*ny*nt*sizeof(float));
	cudaMalloc(&sd_phi_d,nv*sizeof(float));
	cudaMalloc(&y_det_d,nd*sizeof(float));
	cudaMalloc(&id_X_d,nv*sizeof(int));

	cudaMemcpy(X_d,X,nx*ny*nt*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(sd_phi_d,sd_phi,nv*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(y_det_d,y_det,nd*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(id_X_d,id_X,nv*sizeof(int),cudaMemcpyHostToDevice);

	dim3 dimBlock(BLOCK_SIZE_x,BLOCK_SIZE_y);
	dim3 dimGrid_t((nv+dimBlock.x-1)/dimBlock.x,(nd+dimBlock.y-1)/dimBlock.y);
	Ax_fan_mf_gpu_kernel_siddon<<<dimGrid_t, dimBlock>>>(X_d,y_d,SO,OD,scale,nx,ny,nv,sd_phi_d,nd,y_det_d,id_X_d);

	cudaMemcpy(y,y_d,nv*nd*sizeof(float),cudaMemcpyDeviceToHost);

    cudaFree(y_d);cudaFree(X_d);cudaFree(sd_phi_d);cudaFree(y_det_d);cudaFree(id_X_d);
}
