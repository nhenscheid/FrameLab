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

extern "C" void Ax_cone_mf_gpu_new(float *X,float *y,float *sd_phi,float *sd_z,float *y_det,float *z_det,int *id_Y,int *Nv,
float SO,float OD,float scale,float dz,int nx,int ny,int nz,int nt,int na,int nb,int nv,int tmp_size,int nv_block);

__global__ void Ax_cone_mf_gpu_kernel_new(float *x,float *y,float *sd_phi,float *sd_z,float *y_det,float *z_det,int *id_Y,
float SO,float OD,float scale,float dz,int nx,int ny,int nz,int na,int nb,int nv2)
// Please note that this version has O(Nx) per thread, since GPU threads are already saturated.
// O(1) per thread can be achieved by parallelizing the "for" loop here, given sufficient number of GPU threads.
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx0=threadIdx.x;
	int ty0=threadIdx.y;

	int ia=bx*BLOCK_SIZE_x+tx0;
	int ib2=by*BLOCK_SIZE_y+ty0;

    if(ia<na&&ib2<nb*nv2)
    {
		int nd,nx2,ny2,nz2,iv0,iv,ib,id,ix,iy,iz,cx1,cx2,cy1,cy2,cz1,cz2;
		float cos_phi,sin_phi,x1,y1,x2,y2,z1,z2,xx1,yy1,zz1,xx2,yy2,zz2,slope1,slope2,l,d,tmp,rx,ry,rz;

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

		y[id]=0;
        // assuming z1-z2 is small
        if(ABS(x1-x2)>ABS(y1-y2))
        {   slope1=(y2-y1)/(x2-x1);
            slope2=(z2-z1)/(x2-x1);
            for(ix=0;ix<nx;ix++)
            {   xx1=(float)(ix-nx2);xx2=xx1+1;
                if(slope1>=0)
                {   yy1=y1+slope1*(xx1-x1)+ny2;
                    yy2=y1+slope1*(xx2-x1)+ny2;
                }
                else
                {   yy1=y1+slope1*(xx2-x1)+ny2;
                    yy2=y1+slope1*(xx1-x1)+ny2;
                }
                cy1=(int)floor(yy1);
                cy2=(int)floor(yy2);
                if(slope2>=0)
                {   zz1=(z1+slope2*(xx1-x1))/dz+nz2;
                    zz2=(z1+slope2*(xx2-x1))/dz+nz2;
                }
                else
                {   zz1=(z1+slope2*(xx2-x1))/dz+nz2;
                    zz2=(z1+slope2*(xx1-x1))/dz+nz2;
                }
                cz1=(int)floor(zz1);
                cz2=(int)floor(zz2);

                if(cy2==cy1)
                {   if(cy1>=0&&cy1<=ny-1)
                    {   if(cz2==cz1)
                        {   if(cz1>=0&&cz1<=nz-1)// 11
                            {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                iy=cy1;iz=cz1;y[id]+=l*x[iz*ny*nx+iy*nx+ix];
                            }
                        }
                        else
                        {   if(cz2>0&&cz2<nz)// 12
                            {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                rz=(cz2-zz1)/(zz2-zz1);
                                iy=cy1;iz=cz1;y[id]+=rz*l*x[iz*ny*nx+iy*nx+ix];
                                iy=cy1;iz=cz2;y[id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];
                            }
                            else
                            {   if(cz2==0)// 13
                                {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    rz=(cz2-zz1)/(zz2-zz1);
                                    iy=cy1;iz=cz2;y[id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                }
                                if(cz2==nz)// 14
                                {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    rz=(cz2-zz1)/(zz2-zz1);
                                    iy=cy1;iz=cz1;y[id]+=rz*l*x[iz*ny*nx+iy*nx+ix];
                                }
                            }
                        }
                    }
                }
                else
                {   if(cy2>0&&cy2<ny)
                    {   if(cz2==cz1)
                        {   if(cz1>=0&&cz1<=nz-1)// 21
                            {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                ry=(cy2-yy1)/d;
                                iy=cy1;iz=cz1;y[id]+=ry*l*x[iz*ny*nx+iy*nx+ix];
                                iy=cy2;iz=cz1;y[id]+=(1-ry)*l*x[iz*ny*nx+iy*nx+ix];
                            }
                        }
                        else
                        {   if(cz2>0&&cz2<nz)// 22
                            {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                ry=(cy2-yy1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                if(ry>rz)
                                {   iy=cy1;iz=cz1;y[id]+=rz*l*x[iz*ny*nx+iy*nx+ix];
                                    iy=cy1;iz=cz2;y[id]+=(ry-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                    iy=cy2;iz=cz2;y[id]+=(1-ry)*l*x[iz*ny*nx+iy*nx+ix];
                                }
                                else
                                {   iy=cy1;iz=cz1;y[id]+=ry*l*x[iz*ny*nx+iy*nx+ix];
                                    iy=cy2;iz=cz1;y[id]+=(rz-ry)*l*x[iz*ny*nx+iy*nx+ix];
                                    iy=cy2;iz=cz2;y[id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                }
                            }
                            else
                            {   if(cz2==0)// 23
                                {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    ry=(cy2-yy1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                    if(ry>rz)
                                    {   iy=cy1;iz=cz2;y[id]+=(ry-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                        iy=cy2;iz=cz2;y[id]+=(1-ry)*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                    else
                                    {   iy=cy2;iz=cz2;y[id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                }
                                if(cz2==nz)// 24
                                {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    ry=(cy2-yy1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                    if(ry>rz)
                                    {   iy=cy1;iz=cz1;y[id]+=rz*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                    else
                                    {   iy=cy1;iz=cz1;y[id]+=ry*l*x[iz*ny*nx+iy*nx+ix];
                                        iy=cy2;iz=cz1;y[id]+=(rz-ry)*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                }
                            }
                        }
                    }
                    else
                    {   if(cy2==0)
                        {   if(cz2==cz1)
                            {   if(cz1>=0&&cz1<=nz-1)// 31
                                {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    ry=(cy2-yy1)/d;
                                    iy=cy2;iz=cz1;y[id]+=(1-ry)*l*x[iz*ny*nx+iy*nx+ix];
                                }
                            }
                            else
                            {   if(cz2>0&&cz2<nz)// 32
                                {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    ry=(cy2-yy1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                    if(ry>rz)
                                    {   iy=cy2;iz=cz2;y[id]+=(1-ry)*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                    else
                                    {   iy=cy2;iz=cz1;y[id]+=(rz-ry)*l*x[iz*ny*nx+iy*nx+ix];
                                        iy=cy2;iz=cz2;y[id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                }
                                else
                                {   if(cz2==0)// 33
                                    {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        ry=(cy2-yy1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                        if(ry>rz)
                                        {   iy=cy2;iz=cz2;y[id]+=(1-ry)*l*x[iz*ny*nx+iy*nx+ix];
                                        }
                                        else
                                        {   iy=cy2;iz=cz2;y[id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                        }
                                    }
                                    if(cz2==nz)// 34
                                    {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        ry=(cy2-yy1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                        if(ry>rz)
                                        {
                                        }
                                        else
                                        {   iy=cy2;iz=cz1;y[id]+=(rz-ry)*l*x[iz*ny*nx+iy*nx+ix];
                                        }
                                    }
                                }
                            }
                        }

                        if(cy2==ny)
                        {   if(cz2==cz1)
                            {   if(cz1>=0&&cz1<=nz-1)// 41
                                {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    ry=(cy2-yy1)/d;
                                    iy=cy1;iz=cz1;y[id]+=ry*l*x[iz*ny*nx+iy*nx+ix];
                                }
                            }
                            else
                            {   if(cz2>0&&cz2<nz)// 42
                                {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    ry=(cy2-yy1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                    if(ry>rz)
                                    {   iy=cy1;iz=cz1;y[id]+=rz*l*x[iz*ny*nx+iy*nx+ix];
                                        iy=cy1;iz=cz2;y[id]+=(ry-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                    else
                                    {   iy=cy1;iz=cz1;y[id]+=ry*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                }
                                else
                                {   if(cz2==0)// 43
                                    {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        ry=(cy2-yy1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                        if(ry>rz)
                                        {   iy=cy1;iz=cz2;y[id]+=(ry-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                        }
                                        else
                                        {
                                        }
                                    }
                                    if(cz2==nz)// 44
                                    {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        ry=(cy2-yy1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                        if(ry>rz)
                                        {   iy=cy1;iz=cz1;y[id]+=rz*l*x[iz*ny*nx+iy*nx+ix];
                                        }
                                        else
                                        {   iy=cy1;iz=cz1;y[id]+=ry*l*x[iz*ny*nx+iy*nx+ix];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {   slope1=(x2-x1)/(y2-y1);
            slope2=(z2-z1)/(y2-y1);
            for(iy=0;iy<ny;iy++)
            {   yy1=(float)(iy-ny2);yy2=yy1+1;
                if(slope1>=0)
                {   xx1=x1+slope1*(yy1-y1)+nx2;
                    xx2=x1+slope1*(yy2-y1)+nx2;
                }
                else
                {   xx1=x1+slope1*(yy2-y1)+nx2;
                    xx2=x1+slope1*(yy1-y1)+nx2;
                }
                cx1=(int)floor(xx1);
                cx2=(int)floor(xx2);
                if(slope2>=0)
                {   zz1=(z1+slope2*(yy1-y1))/dz+nz2;
                    zz2=(z1+slope2*(yy2-y1))/dz+nz2;
                }
                else
                {   zz1=(z1+slope2*(yy2-y1))/dz+nz2;
                    zz2=(z1+slope2*(yy1-y1))/dz+nz2;
                }
                cz1=(int)floor(zz1);
                cz2=(int)floor(zz2);

                if(cx2==cx1)
                {   if(cx1>=0&&cx1<=nx-1)
                    {   if(cz2==cz1)
                        {   if(cz1>=0&&cz1<=nz-1)// 11
                            {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                ix=cx1;iz=cz1;y[id]+=l*x[iz*ny*nx+iy*nx+ix];
                            }
                        }
                        else
                        {   if(cz2>0&&cz2<nz)// 12
                            {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                rz=(cz2-zz1)/(zz2-zz1);
                                ix=cx1;iz=cz1;y[id]+=rz*l*x[iz*ny*nx+iy*nx+ix];
                                ix=cx1;iz=cz2;y[id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];
                            }
                            else
                            {   if(cz2==0)// 13
                                {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    rz=(cz2-zz1)/(zz2-zz1);
                                    ix=cx1;iz=cz2;y[id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                }
                                if(cz2==nz)// 14
                                {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    rz=(cz2-zz1)/(zz2-zz1);
                                    ix=cx1;iz=cz1;y[id]+=rz*l*x[iz*ny*nx+iy*nx+ix];
                                }
                            }
                        }
                    }
                }
                else
                {   if(cx2>0&&cx2<nx)
                    {   if(cz2==cz1)
                        {   if(cz1>=0&&cz1<=nz-1)// 21
                            {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                rx=(cx2-xx1)/d;
                                ix=cx1;iz=cz1;y[id]+=rx*l*x[iz*ny*nx+iy*nx+ix];
                                ix=cx2;iz=cz1;y[id]+=(1-rx)*l*x[iz*ny*nx+iy*nx+ix];
                            }
                        }
                        else
                        {   if(cz2>0&&cz2<nz)// 22
                            {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                rx=(cx2-xx1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                if(rx>rz)
                                {   ix=cx1;iz=cz1;y[id]+=rz*l*x[iz*ny*nx+iy*nx+ix];
                                    ix=cx1;iz=cz2;y[id]+=(rx-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                    ix=cx2;iz=cz2;y[id]+=(1-rx)*l*x[iz*ny*nx+iy*nx+ix];
                                }
                                else
                                {   ix=cx1;iz=cz1;y[id]+=rx*l*x[iz*ny*nx+iy*nx+ix];
                                    ix=cx2;iz=cz1;y[id]+=(rz-rx)*l*x[iz*ny*nx+iy*nx+ix];
                                    ix=cx2;iz=cz2;y[id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                }
                            }
                            else
                            {   if(cz2==0)// 23
                                {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    rx=(cx2-xx1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                    if(rx>rz)
                                    {   ix=cx1;iz=cz2;y[id]+=(rx-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                        ix=cx2;iz=cz2;y[id]+=(1-rx)*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                    else
                                    {   ix=cx2;iz=cz2;y[id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                }
                                if(cz2==nz)// 24
                                {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    rx=(cx2-xx1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                    if(rx>rz)
                                    {   ix=cx1;iz=cz1;y[id]+=rz*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                    else
                                    {   ix=cx1;iz=cz1;y[id]+=rx*l*x[iz*ny*nx+iy*nx+ix];
                                        ix=cx2;iz=cz1;y[id]+=(rz-rx)*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                }
                            }
                        }
                    }
                    else
                    {   if(cx2==0)
                        {   if(cz2==cz1)
                            {   if(cz1>=0&&cz1<=nz-1)// 31
                                {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    rx=(cx2-xx1)/d;
                                    ix=cx2;iz=cz1;y[id]+=(1-rx)*l*x[iz*ny*nx+iy*nx+ix];
                                }
                            }
                            else
                            {   if(cz2>0&&cz2<nz)// 32
                                {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    rx=(cx2-xx1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                    if(rx>rz)
                                    {   ix=cx2;iz=cz2;y[id]+=(1-rx)*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                    else
                                    {   ix=cx2;iz=cz1;y[id]+=(rz-rx)*l*x[iz*ny*nx+iy*nx+ix];
                                        ix=cx2;iz=cz2;y[id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                }
                                else
                                {   if(cz2==0)// 33
                                    {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        rx=(cx2-xx1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                        if(rx>rz)
                                        {   ix=cx2;iz=cz2;y[id]+=(1-rx)*l*x[iz*ny*nx+iy*nx+ix];
                                        }
                                        else
                                        {   ix=cx2;iz=cz2;y[id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                        }
                                    }
                                    if(cz2==nz)// 34
                                    {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        rx=(cx2-xx1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                        if(rx>rz)
                                        {
                                        }
                                        else
                                        {   ix=cx2;iz=cz1;y[id]+=(rz-rx)*l*x[iz*ny*nx+iy*nx+ix];
                                        }
                                    }
                                }
                            }
                        }

                        if(cx2==nx)
                        {   if(cz2==cz1)
                            {   if(cz1>=0&&cz1<=nz-1)// 41
                                {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    rx=(cx2-xx1)/d;
                                    ix=cx1;iz=cz1;y[id]+=rx*l*x[iz*ny*nx+iy*nx+ix];
                                }
                            }
                            else
                            {   if(cz2>0&&cz2<nz)// 42
                                {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                    rx=(cx2-xx1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                    if(rx>rz)
                                    {   ix=cx1;iz=cz1;y[id]+=rz*l*x[iz*ny*nx+iy*nx+ix];
                                        ix=cx1;iz=cz2;y[id]+=(rx-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                    else
                                    {   ix=cx1;iz=cz1;y[id]+=rx*l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                }
                                else
                                {   if(cz2==0)// 43
                                    {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        rx=(cx2-xx1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                        if(rx>rz)
                                        {   ix=cx1;iz=cz2;y[id]+=(rx-rz)*l*x[iz*ny*nx+iy*nx+ix];
                                        }
                                        else
                                        {
                                        }
                                    }
                                    if(cz2==nz)// 44
                                    {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        rx=(cx2-xx1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                        if(rx>rz)
                                        {   ix=cx1;iz=cz1;y[id]+=rz*l*x[iz*ny*nx+iy*nx+ix];
                                        }
                                        else
                                        {   ix=cx1;iz=cz1;y[id]+=rx*l*x[iz*ny*nx+iy*nx+ix];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        y[id]*=scale;
    }
}

void Ax_cone_mf_gpu_new(float *X,float *y,float *sd_phi,float *sd_z,float *y_det,float *z_det,int *id_Y,int *Nv,
float SO,float OD,float scale,float dz,int nx,int ny,int nz,int nt,int na,int nb,int nv,int tmp_size,int nv_block)
//assuming dx=dy; during the computation, dx=dy=1; after the computation, the "dx=dy=scale" is taken into account.
//the input can be multiple 3D images, e.g., temporally-resolved images from 4DCT.
//the user needs to supply
//      "id_X" -- the image index for each projection
//      "id_Y" -- the projection index for each image
//      "nv_block" -- to deal with the GPU time limit
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
			Ax_cone_mf_gpu_kernel_new<<<dimGrid_t, dimBlock>>>(x_d,y_d,sd_phi_d,sd_z_d,y_det_d,z_det_d,&id_Y_d[it*tmp_size+v0],
			SO,OD,scale,dz,nx,ny,nz,na,nb,nv2);
			cudaThreadSynchronize();
			cudaMemcpy(&y[id*nd],y_d,nv2*nd*sizeof(float),cudaMemcpyDeviceToHost);
			v0+=nv2;id+=nv2;
		}
	}

    cudaFree(y_d);cudaFree(x_d);cudaFree(sd_phi_d);cudaFree(sd_z_d);cudaFree(y_det_d);cudaFree(z_det_d);cudaFree(id_Y_d);
}
