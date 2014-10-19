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
#include "find_l_3d_d.h"
#define MAX(a,b) (a>b?a:b)
#define MIN(a,b) (a<b?a:b)

#define BLOCK_SIZE_x 16
#define BLOCK_SIZE_y 16

extern "C" void Atx_cone_mf_gpu_new(float *X,float *y,float *cos_phi,float *sin_phi,float *sd_z,float *y_det,float *z_det,int *id_Y,int *Nv,
float SO,float OD,float scale,float dy_det,float y_os,float dz_det,float dz,int nx,int ny,int nz,int nt,int na,int nb,int nv,int tmp_size,int nv_block);

__global__ void Atx_cone_mf_gpu_new_kernel(float *x,float *y,float *cos_phi,float *sin_phi,float *sd_z,float *y_det,float *z_det,int *id_Y,
float SO,float OD,float scale,float dy_det,float y_os,float dz_det,float dz,int nx,int ny,int nz,int na,int nb,int nv2)
// Please note that this version has O(Nv) per thread, since GPU threads are already saturated.
// O(1) per thread can be achieved by parallelizing the "for" loop here, given sufficient number of GPU threads.
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;

	int ix=bx*BLOCK_SIZE_x+tx;
	int iy2=by*BLOCK_SIZE_y+ty;

    if(ix<nx&&iy2<ny*nz)
    {   int nx2,ny2,nz2,na2,nb2,iv0,iv,ia,ib,iy,iz,na_min,na_max,nb_min,nb_max,idx;
		float xc,yc,zc,xr,yr,SD,l,tmp,x1,y1,z1,x2,y2,z2,d;

		SD=SO+OD;
		na2=na/2;nb2=nb/2;
		nx2=nx/2;ny2=ny/2;nz2=nz/2;
		d=(float)sqrt((1+dz*dz)/2);

		iz=(int)floor((float)iy2/(float)ny);
		iy=iy2-iz*ny;
		idx=iz*ny*nx+iy*nx+ix;

		zc=(float)(iz+0.5-nz2)*dz;
		yc=(float)(iy+0.5-ny2);
		xc=(float)(ix+0.5-nx2);

        for(iv0=0;iv0<nv2;iv0++)
        {   iv=id_Y[iv0];
            xr=cos_phi[iv]*xc+sin_phi[iv]*yc;
            yr=-sin_phi[iv]*xc+cos_phi[iv]*yc;

            tmp=SD/((xr+SO)*dy_det);
            na_max=(int)floor((yr+1)*tmp-y_os+na2);
            na_min=(int)floor((yr-1)*tmp-y_os+na2);

            tmp=SD/((xr+SO)*dz_det);
            nb_max=(int)floor((zc+d)*tmp+nb2);
            nb_min=(int)floor((zc-d)*tmp+nb2);

            for(ib=MAX(0,nb_min);ib<=MIN(nb_max,nb-1);ib++)
            {   for(ia=MAX(0,na_min);ia<=MIN(na_max,na-1);ia++)
                {   x1=cos_phi[iv]*(-SO);
                    y1=sin_phi[iv]*(-SO);
                    z1=sd_z[iv];
                    x2=cos_phi[iv]*OD-sin_phi[iv]*y_det[ia];
                    y2=sin_phi[iv]*OD+cos_phi[iv]*y_det[ia];
                    z2=z_det[ib]+sd_z[iv];
                    l=find_l_3d(x1,y1,z1,x2,y2,z2,1.0,1.0,dz,xc,yc,zc);
                    x[idx]+=l*y[iv0*na*nb+ib*na+ia];
                }
            }
        }
    }
}

__global__ void set2zero(float *x,int nx,int nyz)
{	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;

	int ix=bx*BLOCK_SIZE_x+tx;
	int iy=by*BLOCK_SIZE_y+ty;

    if(ix<nx&&iy<nyz)
    {x[iy*nx+ix]=0;}
}

__global__ void scalex(float *x,int nx,int nyz,float scale)
{	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;

	int ix=bx*BLOCK_SIZE_x+tx;
	int iy=by*BLOCK_SIZE_y+ty;

    if(ix<nx&&iy<nyz)
    {x[iy*nx+ix]*=scale;}
}

void Atx_cone_mf_gpu_new(float *X,float *y,float *cos_phi,float *sin_phi,float *sd_z,float *y_det,float *z_det,int *id_Y,int *Nv,
float SO,float OD,float scale,float dy_det,float y_os,float dz_det,float dz,int nx,int ny,int nz,int nt,int na,int nb,int nv,int tmp_size,int nv_block)
// A new method for computing the adjoint X-ray transform (infinitely-narrow beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
{   float *x_d,*y_d,*sd_z_d,*y_det_d,*z_det_d,*cos_phi_d,*sin_phi_d;
	int *id_Y_d,nd,N,id,v0,it,i,n,nv2;

	N=nx*ny*nz;
	nd=na*nb;

	cudaMalloc(&y_d,nv_block*nd*sizeof(float));
	cudaMalloc(&x_d,N*sizeof(float));

	cudaMalloc(&cos_phi_d,nv*sizeof(float));cudaMemcpy(cos_phi_d,cos_phi,nv*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&sin_phi_d,nv*sizeof(float));cudaMemcpy(sin_phi_d,sin_phi,nv*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&sd_z_d,nv*sizeof(float));cudaMemcpy(sd_z_d,sd_z,nv*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&y_det_d,na*sizeof(float));cudaMemcpy(y_det_d,y_det,na*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&z_det_d,nb*sizeof(float));cudaMemcpy(z_det_d,z_det,nb*sizeof(float),cudaMemcpyHostToDevice);
	cudaMalloc(&id_Y_d,nt*tmp_size*sizeof(int));cudaMemcpy(id_Y_d,id_Y,nt*tmp_size*sizeof(int),cudaMemcpyHostToDevice);

	dim3 dimBlock(BLOCK_SIZE_x,BLOCK_SIZE_y);
	dim3 dimGrid_t((nx+dimBlock.x-1)/dimBlock.x,(ny*nz+dimBlock.y-1)/dimBlock.y);

	id=0;
	for(it=0;it<nt;it++)
	{	n=(Nv[it]+nv_block-1)/nv_block;
		set2zero<<<dimGrid_t, dimBlock>>>(x_d,nx,ny*nz);
		v0=0;
		for(i=0;i<n;i++)
		{	if(i<n-1)
			{nv2=nv_block;}
			else
			{nv2=Nv[it]-nv_block*(n-1);}
            cudaMemcpy(y_d,&y[id*nd],nv2*nd*sizeof(float),cudaMemcpyHostToDevice);
			Atx_cone_mf_gpu_new_kernel<<<dimGrid_t, dimBlock>>>(x_d,y_d,cos_phi_d,sin_phi_d,sd_z_d,y_det_d,z_det_d,&id_Y_d[it*tmp_size+v0],
			SO,OD,scale,dy_det,y_os,dz_det,dz,nx,ny,nz,na,nb,nv2);
			cudaThreadSynchronize();
			v0+=nv2;id+=nv2;
		}
		scalex<<<dimGrid_t, dimBlock>>>(x_d,nx,ny*nz,scale);
		cudaMemcpy(&X[it*N],x_d,N*sizeof(float),cudaMemcpyDeviceToHost);
	}

    cudaFree(x_d);cudaFree(y_d);cudaFree(cos_phi_d);cudaFree(sin_phi_d);cudaFree(sd_z_d);cudaFree(y_det_d);cudaFree(z_det_d);cudaFree(id_Y_d);
}
