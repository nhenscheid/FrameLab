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

void Ax_cone_mf_cpu_siddon(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nz,int nv,float *sd_phi,float *sd_z,int na,int nb,float *y_det,float *z_det,int *id_X,float dz0)
// Siddon's method to compute the X-ray transform
{   int Ji,nd,n,nx2,ny2,nz2,iv,id,ia,ib,i,j,nxy,n_alpha,tx,ty,tz;
    float *x,*alpha_x,*alpha_y,*alpha_z,*alpha_xy,*alpha;
    float Si,cos_phi,sin_phi,x1,y1,z1,x2,y2,z2,dx,dy,dz,totallength,alpha_c;

    alpha_x=(float*)malloc((nx+1)*sizeof(float));
    alpha_y=(float*)malloc((ny+1)*sizeof(float));
    alpha_z=(float*)malloc((nz+1)*sizeof(float));
    alpha_xy=(float*)malloc((nx+ny+2)*sizeof(float));
    alpha=(float*)malloc((nx+ny+nz+3)*sizeof(float));

    nx2=nx/2;ny2=ny/2;nz2=nz/2;
    n=nx*ny*nz;
    nd=na*nb;

    for(iv=0;iv<nv;iv++)
    {   x=&X[id_X[iv]*n];// id_X[iv] refers to the image frame index for View "iv"
        cos_phi=(float)cos(sd_phi[iv]);sin_phi=(float)sin(sd_phi[iv]);
        x1=cos_phi*(-SO);// x coordinate of source
        y1=sin_phi*(-SO);// y coordinate of source
        z1=sd_z[iv];// z coordinate of source
        for(ib=0;ib<nb;ib++)
        {   for(ia=0;ia<na;ia++)
            {   x2=cos_phi*OD-sin_phi*y_det[ia];// x coordinate of detector
                y2=sin_phi*OD+cos_phi*y_det[ia];// y coordinate of detector
                z2=z_det[ib]+sd_z[iv];// z coordinate of detector
                id=ib*na+ia;
                // step 1: compute alpha_x and alpha_y
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
                    {alpha_z[j]=(i*dz0-z1)/dz;}
                }
                else if(z1>z2)
                {   for(i=nz2,j=0;i>=-nz2;i--,j++)
                    {alpha_z[j]=(i*dz0-z1)/dz;}
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
                totallength=(float)sqrt(dx*dx+dy*dy+dz*dz);
                y[iv*nd+id]=0;
                for(i=0;i<n_alpha-1;i++)
                {   alpha_c=(alpha[i]+alpha[i+1])/2;
                    tx=(int)floor(x1+alpha_c*dx+nx2);
                    if(tx>=0&&tx<nx)
                    {   ty=(int)floor(y1+alpha_c*dy+ny2);
                        if(ty>=0&&ty<ny)
                        {   tz=(int)floor((z1+alpha_c*dz)/dz0+nz2);
                            if(tz>=0&&tz<nz)
                            {   Ji=tz*ny*nx+ty*nx+tx;
                                Si=(alpha[i+1]-alpha[i])*totallength;
                                y[iv*nd+id]+=Si*x[Ji];
                            }
                        }
                    }
                }
                y[iv*nd+id]*=scale;
            }
        }
    }
    free(alpha);free(alpha_x);free(alpha_y);free(alpha_z);free(alpha_xy);
}
