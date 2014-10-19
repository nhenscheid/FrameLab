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
#include <malloc.h>
#include "sort_alpha_d.h"

#define BLOCK_SIZE_x 16
#define BLOCK_SIZE_y 16
#define Nax 256+1 // Please note that you need increase the value of Nax manually, if Nx>256, since the thread doesn't support malloc.
#define Nay 256+1 // Please note that you need increase the value of Nax manually, if Ny>256, since the thread doesn't support malloc.
#define Naz 192+1 // Please note that you need increase the value of Nax manually, if Nz>192, since the thread doesn't support malloc.
#define Naxy 256+256+2 // Please note that you need increase the value of Naxy manually, if Nx+Ny>512, since the thread doesn't support malloc.
#define Na 256+256+192+3 // Please note that you need increase the value of Na manually, if Nx+Ny+Nz>704, since the thread doesn't support malloc.

extern "C" void Ax_cone_mf_gpu_siddon(float *X,float *y,float *sd_phi,float *sd_z,float *y_det,float *z_det,int *id_Y,int *Nv,
float SO,float OD,float scale,float dz,int nx,int ny,int nz,int nt,int na,int nb,int nv,int tmp_size,int nv_block);

__global__ void Ax_cone_mf_gpu_kernel_siddon(float *x,float *y,float *sd_phi,float *sd_z,float *y_det,float *z_det,int *id_Y,
float SO,float OD,float scale,float dz0,int nx,int ny,int nz,int na,int nb,int nv2)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx0=threadIdx.x;
	int ty0=threadIdx.y;

	int ia=bx*BLOCK_SIZE_x+tx0;
	int ib2=by*BLOCK_SIZE_y+ty0;

    if(ia<na&&ib2<nb*nv2)
    {
		int Ji,nd,nx2,ny2,nz2,iv0,iv,ib,id,i,j,nxy,n_alpha,tx,ty,tz;
		float alpha_x[Nax],alpha_y[Nay],alpha_z[Naz],alpha_xy[Naxy],alpha[Na];
		float Si,cos_phi,sin_phi,x1,y1,z1,x2,y2,z2,dx,dy,dz,totallength,alpha_c;

		nx2=nx/2;ny2=ny/2;nz2=nz/2;
		nd=na*nb;

		iv0=(int)floor((float)ib2/(float)nb);
		ib=ib2-iv0*nb;
		iv=id_Y[iv0];
		id=iv0*nd+ib*na+ia;

        cos_phi=cosf(sd_phi[iv]);sin_phi=sinf(sd_phi[iv]);
        x1=cos_phi*(-SO);
        y1=sin_phi*(-SO);
        z1=sd_z[iv];
		x2=cos_phi*OD-sin_phi*y_det[ia];
        y2=sin_phi*OD+cos_phi*y_det[ia];
		z2=z_det[ib]+sd_z[iv];

        // step 1: compute alpha_x, alpha_y and alpha_z
        dx=x2-x1;dy=y2-y1;dz=z2-z1;
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

        if(z1<z2)
        {   for(i=-nz2,j=0;i<=nz2;i++,j++)
            {alpha_z[j]=(i*dz0-z1)/dz0;}
        }
        else if(z1>z2)
        {   for(i=nz2,j=0;i>=-nz2;i--,j++)
            {alpha_z[j]=(i*dz0-z1)/dz0;}
        }
        // step 2: sort and bin alpha_x, alpha_y, alpha_z into alpha
        if(z1==z2)
        {   if(y1==y2)
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
        }
        else
        {   if(y1==y2)
            {   if(x1==x2)
                {   for(i=0;i<=nz;i++)
                    {alpha[i]=alpha_z[i];}
                    n_alpha=nz+1;
                }
                else
                {n_alpha=sort_alpha(nx+1,nz+1,alpha_x,alpha_z,alpha);}
            }
            else
            {   if(x1==x2)
                {n_alpha=sort_alpha(ny+1,nz+1,alpha_y,alpha_z,alpha);}
                else
                {   nxy=sort_alpha(nx+1,ny+1,alpha_x,alpha_y,alpha_xy);
                    n_alpha=sort_alpha(nxy,nz+1,alpha_xy,alpha_z,alpha);
                }
            }
        }
        // step 3: compute x-ray transform y[iy]
        totallength=sqrtf(dx*dx+dy*dy+dz*dz);
        y[id]=0;
        for(i=0;i<n_alpha-1;i++)
        {   alpha_c=(float)0.5*(alpha[i]+alpha[i+1]);
            tx=(int)floorf(x1+alpha_c*dx+nx2);
            if(tx>=0&&tx<nx)
            {   ty=(int)floorf(y1+alpha_c*dy+ny2);
                if(ty>=0&&ty<ny)
                {   tz=(int)floorf((z1+alpha_c*dz)/dz0+nz2);
                    if(tz>=0&&tz<nz)
                    {   Ji=tz*ny*nx+ty*nx+tx;
                        Si=(alpha[i+1]-alpha[i])*totallength;
                        y[id]+=Si*x[Ji];
                    }
                }
            }
        }
        y[id]*=scale;
    }
}

void Ax_cone_mf_gpu_siddon(float *X,float *y,float *sd_phi,float *sd_z,float *y_det,float *z_det,int *id_Y,int *Nv,
float SO,float OD,float scale,float dz,int nx,int ny,int nz,int nt,int na,int nb,int nv,int tmp_size,int nv_block)
// Siddon's method to compute the X-ray transform
{   float *y_d,*x_d,*sd_phi_d,*sd_z_d,*y_det_d,*z_det_d;
	int *id_Y_d,nd,N,id,v0,it,i,n,nv2;

	N=nx*ny*nz;
	nd=na*nb;

	cudaMalloc(&y_d,nv_block*nd*sizeof(float));
	cudaMalloc(&x_d,N*sizeof(float));

	cudaMalloc(&sd_phi_d,nv*sizeof(float));cudaMemcpy(sd_phi_d,sd_phi,nv*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&sd_z_d,nv*sizeof(float));cudaMemcpy(sd_z_d,sd_z,nv*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&y_det_d,na*sizeof(float));cudaMemcpy(y_det_d,y_det,na*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&z_det_d,nb*sizeof(float));cudaMemcpy(z_det_d,z_det,nb*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&id_Y_d,nt*tmp_size*sizeof(int));cudaMemcpy(id_Y_d,id_Y,nt*tmp_size*sizeof(int),cudaMemcpyHostToDevice);

	dim3 dimBlock(BLOCK_SIZE_x,BLOCK_SIZE_y);

	id=0;
	for(it=0;it<nt;it++)
	{	cudaMemcpy(x_d,&X[it*N],N*sizeof(float),cudaMemcpyHostToDevice);
		n=(Nv[it]+nv_block-1)/nv_block;
		v0=0;
		for(i=0;i<n;i++)
		{	if(i<n-1)
			{nv2=nv_block;}
			else
			{nv2=Nv[it]-nv_block*(n-1);}
			dim3 dimGrid_t((na+dimBlock.x-1)/dimBlock.x,(nv2*nb+dimBlock.y-1)/dimBlock.y);
			Ax_cone_mf_gpu_kernel_siddon<<<dimGrid_t, dimBlock>>>(x_d,y_d,sd_phi_d,sd_z_d,y_det_d,z_det_d,&id_Y_d[it*tmp_size+v0],
			SO,OD,scale,dz,nx,ny,nz,na,nb,nv2);
			cudaThreadSynchronize();
			cudaMemcpy(&y[id*nd],y_d,nv2*nd*sizeof(float),cudaMemcpyDeviceToHost);
			v0+=nv2;id+=nv2;
		}
	}
    cudaFree(y_d);cudaFree(x_d);cudaFree(sd_phi_d);cudaFree(sd_z_d);cudaFree(y_det_d);cudaFree(z_det_d);cudaFree(id_Y_d);
}
