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
#include "sort_alpha.h"

void Ax_fan_mf_cpu_siddon(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X)
// Siddon's method to compute the X-ray transform
{   int Ji,n,nx2,ny2,iv,id,i,j,n_alpha,tx,ty;
    float *x,*alpha_x,*alpha_y,*alpha;
    float Si,cos_phi,sin_phi,x1,y1,x2,y2,dx,dy,totallength,alpha_c;

    alpha_x=(float*)malloc((nx+1)*sizeof(float));
    alpha_y=(float*)malloc((ny+1)*sizeof(float));
    alpha=(float*)malloc((nx+ny+2)*sizeof(float));

    nx2=nx/2;ny2=ny/2;
    n=nx*ny;

    for(iv=0;iv<nv;iv++)
    {   x=&X[id_X[iv]*n];// id_X[iv] refers to the image frame index for View "iv"
        cos_phi=(float)cos(sd_phi[iv]);sin_phi=(float)sin(sd_phi[iv]);
        x1=cos_phi*(-SO);// x coordinate of source
        y1=sin_phi*(-SO);// y coordinate of source
        for(id=0;id<nd;id++)
        {   x2=cos_phi*OD-sin_phi*y_det[id];// x coordinate of detector
            y2=sin_phi*OD+cos_phi*y_det[id];// y coordinate of detector
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
            totallength=(float)sqrt(dx*dx+dy*dy);
            y[iv*nd+id]=0;
            for(i=0;i<n_alpha-1;i++)
            {   alpha_c=(float)0.5*(alpha[i]+alpha[i+1]);
                tx=(int)floor(x1+alpha_c*dx+nx2);
                if(tx>=0&&tx<nx)
                {   ty=(int)floor(y1+alpha_c*dy+ny2);
                    if(ty>=0&&ty<ny)
                    {   Ji=ty*nx+tx;// the image pixel index
                        Si=(alpha[i+1]-alpha[i])*totallength;// the intersection length
                        y[iv*nd+id]+=Si*x[Ji];// sum of the length weighted by the attenuation coefficient
                    }
                }
            }
            y[iv*nd+id]*=scale;// take the length scale into account
        }
    }
    free(alpha);free(alpha_x);free(alpha_y);
}
