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
#define ABS(a) (a>0?a:-(a))

#define BLOCK_SIZE_x 16
#define BLOCK_SIZE_y 16

extern "C" void Ax_fan_mf_gpu_new(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X,int nt);

__global__ void Ax_fan_mf_gpu_new_kernel(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X)
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
		int n,nx2,ny2,ix,iy,c1,c2;
		float *x,cos_phi,sin_phi,x1,y1,x2,y2,xx1,yy1,xx2,yy2,slope,l,d;

		nx2=nx/2;ny2=ny/2;
		n=nx*ny;

		x=&X[id_X[iv]*n];
        cos_phi=(float)cos(sd_phi[iv]);sin_phi=(float)sin(sd_phi[iv]);
        x1=cos_phi*(-SO);
        y1=sin_phi*(-SO);

        x2=cos_phi*OD-sin_phi*y_det[id];
        y2=sin_phi*OD+cos_phi*y_det[id];

        y[iv*nd+id]=0;
        if(ABS(x1-x2)>ABS(y1-y2))
        {   slope=(y2-y1)/(x2-x1);
            for(ix=0;ix<nx;ix++)
            {   xx1=(float)(ix-nx2);xx2=xx1+1;
                if(slope>=0)
                {   yy1=y1+slope*(xx1-x1)+ny2;
                    yy2=y1+slope*(xx2-x1)+ny2;
                }
                else
                {   yy1=y1+slope*(xx2-x1)+ny2;
                    yy2=y1+slope*(xx1-x1)+ny2;
                }
                c1=(int)floor(yy1);
                c2=(int)floor(yy2);

                if(c2==c1)// c1 and c2 differs less than 1
                {   if(c1>=0&&c1<=ny-1)
                    {   d=yy2-yy1;l=(float)sqrt(d*d+1);
                        iy=c1;y[iv*nd+id]+=l*x[iy*nx+ix];
                    }
                }
                else
                {   if(c2>0&&c2<ny)
                    {   d=yy2-yy1;l=(float)sqrt(d*d+1);
                        iy=c1;y[iv*nd+id]+=((c2-yy1)/d)*l*x[iy*nx+ix];
                        iy=c2;y[iv*nd+id]+=((yy2-c2)/d)*l*x[iy*nx+ix];
                    }
                    else
                    {   if(c2==0)
                        {   d=yy2-yy1;l=(float)sqrt(d*d+1);
                            iy=c2;y[iv*nd+id]+=((yy2-c2)/d)*l*x[iy*nx+ix];
                        }
                        if(c2==ny)
                        {   d=yy2-yy1;l=(float)sqrt(d*d+1);
                            iy=c1;y[iv*nd+id]+=((c2-yy1)/d)*l*x[iy*nx+ix];
                        }
                    }
                }
            }
        }
        else
        {   slope=(x2-x1)/(y2-y1);
            for(iy=0;iy<ny;iy++)
            {   yy1=(float)(iy-ny2);yy2=yy1+1;

                if(slope>=0)
                {   xx1=x1+slope*(yy1-y1)+nx2;
                    xx2=x1+slope*(yy2-y1)+nx2;
                }
                else
                {   xx1=x1+slope*(yy2-y1)+nx2;
                    xx2=x1+slope*(yy1-y1)+nx2;
                }
                c1=(int)floor(xx1);
                c2=(int)floor(xx2);

                if(c2==c1)// c1 and c2 differs less than 1
                {   if(c1>=0&&c1<=nx-1)
                    {   d=xx2-xx1;l=(float)sqrt(d*d+1);
                        ix=c1;y[iv*nd+id]+=l*x[iy*nx+ix];
                    }
                }
                else
                {   if(c2>0&&c2<nx)
                    {   d=xx2-xx1;l=(float)sqrt(d*d+1);
                        ix=c1;y[iv*nd+id]+=((c2-xx1)/d)*l*x[iy*nx+ix];
                        ix=c2;y[iv*nd+id]+=((xx2-c2)/d)*l*x[iy*nx+ix];
                    }
                    else
                    {   if(c2==0)
                        {   d=xx2-xx1;l=(float)sqrt(d*d+1);
                            ix=c2;y[iv*nd+id]+=((xx2-c2)/d)*l*x[iy*nx+ix];
                        }
                        if(c2==ny)
                        {   d=xx2-xx1;l=(float)sqrt(d*d+1);
                            ix=c1;y[iv*nd+id]+=((c2-xx1)/d)*l*x[iy*nx+ix];
                        }
                    }
                }
            }
        }
        y[iv*nd+id]*=scale;
    }
}

void Ax_fan_mf_gpu_new(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X,int nt)
// A new method for computing the X-ray transform (infinitely-narrow beam)
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

	Ax_fan_mf_gpu_new_kernel<<<dimGrid_t, dimBlock>>>(X_d,y_d,SO,OD,scale,nx,ny,nv,sd_phi_d,nd,y_det_d,id_X_d);

	cudaMemcpy(y,y_d,nv*nd*sizeof(float),cudaMemcpyDeviceToHost);

    cudaFree(y_d);cudaFree(X_d);cudaFree(sd_phi_d);cudaFree(y_det_d);cudaFree(id_X_d);
}
