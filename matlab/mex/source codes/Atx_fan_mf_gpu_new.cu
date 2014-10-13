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
#include "find_l_d.h"
#define MAX(a,b) (a>b?a:b)
#define MIN(a,b) (a<b?a:b)

#define BLOCK_SIZE_x 16
#define BLOCK_SIZE_y 16

extern "C" void Atx_fan_mf_gpu_new(float *x,float *Y,float SO,float OD,float scale,int nx,int ny,float *sd_phi,int nd,float *y_det,float dy_det,float y_os,
int nt,int *id_Y,int *Nv,int tmp_size,int nv,float *cos_phi,float *sin_phi,float *cos_det,float *sin_det);

__global__ void Atx_fan_mf_gpu_new_kernel(float *x,float *Y,float SO,float OD,float scale,int nx,int ny,float *sd_phi,int nd,float *y_det,float dy_det,float y_os,
int nt,int *id_Y,int *Nv,int tmp_size,float *cos_phi,float *sin_phi,float *cos_det,float *sin_det)
// Please note that this version has O(Nv) per thread, since GPU threads are already saturated.
// O(1) per thread can be achieved by parallelizing the "for" loop here, given sufficient number of GPU threads.
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;

	int ix=bx*BLOCK_SIZE_x+tx;
	int iy2=by*BLOCK_SIZE_y+ty;

    if(ix<nx&&iy2<ny*nt)
    {
		int nx2,ny2,nd2,it,iv,id,iy,nd_min,nd_max,j;
		float *y,xc,yc,xr,yr,SD,Px,Py,Pd,l,tmp;

		SD=SO+OD;

		nx2=nx/2;ny2=ny/2;nd2=nd/2;
		it=(int)floor((float)iy2/(float)ny);
		iy=iy2-it*ny;

		yc=(float)(iy+0.5-ny2);
		xc=(float)(ix+0.5-nx2);

        x[iy2*nx+ix]=0;
        for(j=0;j<Nv[it];j++)
        {   iv=id_Y[it*tmp_size+j];
            y=&Y[iv*nd];

            xr=cos_phi[iv]*xc+sin_phi[iv]*yc;
            yr=-sin_phi[iv]*xc+cos_phi[iv]*yc;
            tmp=SD/((xr+SO)*dy_det);
            nd_max=(int)floor((yr+1)*tmp-y_os+nd2);
            nd_min=(int)floor((yr-1)*tmp-y_os+nd2);
            for(id=MAX(0,nd_min);id<=MIN(nd_max,nd-1);id++)
            {   Px=-(sin_det[id]*cos_phi[iv]+cos_det[id]*sin_phi[iv]);
                Py=cos_det[id]*cos_phi[iv]-sin_det[id]*sin_phi[iv];
                Pd=SO*sin_det[id];
                l=find_l(Px,Py,Pd,xc,yc);
                if(l>0)
                {x[iy2*nx+ix]+=l*y[id];}
            }
        }
        x[iy2*nx+ix]*=scale;
    }
}

void Atx_fan_mf_gpu_new(float *x,float *Y,float SO,float OD,float scale,int nx,int ny,float *sd_phi,int nd,float *y_det,float dy_det,float y_os,
int nt,int *id_Y,int *Nv,int tmp_size,int nv,float *cos_phi,float *sin_phi,float *cos_det,float *sin_det)
// A new method for computing the adjoint X-ray transform (infinitely-narrow beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
{   float *x_d,*Y_d,*sd_phi_d,*y_det_d,*cos_phi_d,*sin_phi_d,*cos_det_d,*sin_det_d;
	int *id_Y_d,*Nv_d;

	cudaMalloc(&cos_phi_d,nv*sizeof(float));cudaMemcpy(cos_phi_d,cos_phi,nv*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&sin_phi_d,nv*sizeof(float));cudaMemcpy(sin_phi_d,sin_phi,nv*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&cos_det_d,nd*sizeof(float));cudaMemcpy(cos_det_d,cos_det,nd*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&sin_det_d,nd*sizeof(float));cudaMemcpy(sin_det_d,sin_det,nd*sizeof(float),cudaMemcpyHostToDevice);

	cudaMalloc(&x_d,nx*ny*nt*sizeof(float));
	cudaMalloc(&Y_d,nv*nd*sizeof(float));cudaMemcpy(Y_d,Y,nv*nd*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&sd_phi_d,nv*sizeof(float));cudaMemcpy(sd_phi_d,sd_phi,nv*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&y_det_d,nd*sizeof(float));cudaMemcpy(y_det_d,y_det,nd*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&id_Y_d,nt*tmp_size*sizeof(int));cudaMemcpy(id_Y_d,id_Y,nt*tmp_size*sizeof(int),cudaMemcpyHostToDevice);
	cudaMalloc(&Nv_d,nt*sizeof(int));cudaMemcpy(Nv_d,Nv,nt*sizeof(int),cudaMemcpyHostToDevice);

	dim3 dimBlock(BLOCK_SIZE_x,BLOCK_SIZE_y);
	dim3 dimGrid_t((nx+dimBlock.x-1)/dimBlock.x,(ny*nt+dimBlock.y-1)/dimBlock.y);
	Atx_fan_mf_gpu_new_kernel<<<dimGrid_t, dimBlock>>>(x_d,Y_d,SO,OD,scale,nx,ny,sd_phi_d,nd,y_det_d,dy_det,y_os,
	nt,id_Y_d,Nv_d,tmp_size,cos_phi_d,sin_phi_d,cos_det_d,sin_det_d);

	cudaMemcpy(x,x_d,nx*ny*nt*sizeof(float),cudaMemcpyDeviceToHost);

    cudaFree(x_d);cudaFree(Y_d);cudaFree(sd_phi_d);cudaFree(y_det_d);cudaFree(id_Y_d);cudaFree(Nv_d);
    cudaFree(cos_phi_d);cudaFree(sin_phi_d);cudaFree(cos_det_d);cudaFree(sin_det_d);
}
