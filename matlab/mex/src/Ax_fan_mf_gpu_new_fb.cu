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
#include "find_area_d.h"
#define ABS(a) (a>0?a:-(a))

#define BLOCK_SIZE_x 16
#define BLOCK_SIZE_y 16

extern "C" void Ax_fan_mf_gpu_new_fb(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X,int nt,float dy_det);

__global__ void Ax_fan_mf_gpu_new_fb_kernel(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X,float dy_det)
// Please note that this version has O(Nx) per thread, since GPU threads are already saturated.
// O(1) per thread can be achieved by parallelizing the "for" loop here, given sufficient number of GPU threads.
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx0=threadIdx.x;
	int ty0=threadIdx.y;

	int iv=bx*BLOCK_SIZE_x+tx0;
	int id=by*BLOCK_SIZE_y+ty0;

    if(iv<nv&&id<nd)
    {
		int n,nx2,ny2,ix,iy,ymin,ymax,xmin,xmax,i,xx[4],yy[4];
		float *x,cos_phi,sin_phi,x1,y1,x2[2],y2[2],xx1,yy1,xx2,yy2,xc,yc,slope[2],SD,Px[2],Py[2],Pd[2],angle_det,l;

		SD=SO+OD;
		nx2=nx/2;ny2=ny/2;
		n=nx*ny;

		x=&X[id_X[iv]*n];
        cos_phi=(float)cos(sd_phi[iv]);sin_phi=(float)sin(sd_phi[iv]);
        x1=cos_phi*(-SO);
        y1=sin_phi*(-SO);

		x2[0]=cos_phi*OD-sin_phi*(y_det[id]-dy_det/2);
        y2[0]=sin_phi*OD+cos_phi*(y_det[id]-dy_det/2);
        angle_det=(float)atan2(y_det[id]-dy_det/2,SD);
        Px[0]=(float)-sin(angle_det+sd_phi[iv]);Py[0]=(float)cos(angle_det+sd_phi[iv]);
        Pd[0]=SO*(float)sin(angle_det);

        x2[1]=cos_phi*OD-sin_phi*(y_det[id]+dy_det/2);
        y2[1]=sin_phi*OD+cos_phi*(y_det[id]+dy_det/2);
        angle_det=(float)atan2(y_det[id]+dy_det/2,SD);
        Px[1]=(float)-sin(angle_det+sd_phi[iv]);Py[1]=(float)cos(angle_det+sd_phi[iv]);
        Pd[1]=SO*(float)sin(angle_det);

        y[iv*nd+id]=0;
        if(ABS(x1-(x2[0]+x2[1])/2)>ABS(y1-(y2[0]+y2[1])/2))
        {   if(x2[0]!=x1&&x2[1]!=x1)
            {   slope[0]=(y2[0]-y1)/(x2[0]-x1);
                slope[1]=(y2[1]-y1)/(x2[1]-x1);
            }
            for(ix=0;ix<nx;ix++)
            {   xc=(float)(ix+0.5-nx2);
                xx1=(float)(ix-nx2);xx2=xx1+1;

                if(x2[0]==x1||x2[1]==x1)
                {   if(ABS(x1-xc)<=0.5)
                    {ymin=0;ymax=ny-1;}
                }
                else
                {   ymin=ny-1;ymax=0;
                    yy[0]=(int)floor(y1+slope[0]*(xx1-x1)+ny2);
                    yy[1]=(int)floor(y1+slope[0]*(xx2-x1)+ny2);
                    yy[2]=(int)floor(y1+slope[1]*(xx1-x1)+ny2);
                    yy[3]=(int)floor(y1+slope[1]*(xx2-x1)+ny2);
                    for(i=0;i<4;i++)
                    {   if(yy[i]<ymin&&yy[i]>=0)
                        {ymin=yy[i];}
                        if(yy[i]>ymax&&yy[i]<=ny-1)
                        {ymax=yy[i];}
                    }
                }

                for(iy=ymin;iy<=ymax;iy++)
                {   yc=(float)(iy+0.5-ny2);
                    l=find_area(Px[0],Py[0],Pd[0],Px[1],Py[1],Pd[1],xc,yc);
                    if(l>0)
                    {y[iv*nd+id]+=l*x[iy*nx+ix];}
                }
            }
        }
        else
        {   if(y2[0]!=y1&&y2[1]!=y1)
            {   slope[0]=(x2[0]-x1)/(y2[0]-y1);
                slope[1]=(x2[1]-x1)/(y2[1]-y1);
            }
            for(iy=0;iy<ny;iy++)
            {   yc=(float)(iy+0.5-ny2);
                yy1=(float)(iy-ny2);yy2=yy1+1;

                if(y2[0]==y1||y2[1]==y1)
                {   if(ABS(y1-yc)<=0.5)
                    {xmin=0;xmax=nx-1;}
                }
                else
                {   xmin=nx-1;xmax=0;
                    xx[0]=(int)floor(x1+slope[0]*(yy1-y1)+nx2);
                    xx[1]=(int)floor(x1+slope[0]*(yy2-y1)+nx2);
                    xx[2]=(int)floor(x1+slope[1]*(yy1-y1)+nx2);
                    xx[3]=(int)floor(x1+slope[1]*(yy2-y1)+nx2);
                    for(i=0;i<4;i++)
                    {   if(xx[i]<xmin&&xx[i]>=0)
                        {xmin=xx[i];}
                        if(xx[i]>xmax&&xx[i]<=nx-1)
                        {xmax=xx[i];}
                    }
                }

                for(ix=xmin;ix<=xmax;ix++)
                {   xc=(float)(ix+0.5-nx2);
                    l=find_area(Px[0],Py[0],Pd[0],Px[1],Py[1],Pd[1],xc,yc);
                    if(l>0)
                    {y[iv*nd+id]+=l*x[iy*nx+ix];}
                }
            }
        }
        y[iv*nd+id]*=scale;
    }
}

void Ax_fan_mf_gpu_new_fb(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X,int nt,float dy_det)
// A new method for computing the X-ray transform (finite-size beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
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

	Ax_fan_mf_gpu_new_fb_kernel<<<dimGrid_t, dimBlock>>>(X_d,y_d,SO,OD,scale,nx,ny,nv,sd_phi_d,nd,y_det_d,id_X_d,dy_det);

	cudaMemcpy(y,y_d,nv*nd*sizeof(float),cudaMemcpyDeviceToHost);

    cudaFree(y_d);cudaFree(X_d);cudaFree(sd_phi_d);cudaFree(y_det_d);cudaFree(id_X_d);
}
