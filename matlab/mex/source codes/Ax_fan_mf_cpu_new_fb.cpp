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
#include "find_area.h"
#define ABS(a) (a>0?a:-(a))

void Ax_fan_mf_cpu_new_fb(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X,float dy_det)
// A new method for computing the X-ray transform (finite-size beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
{   int n,nx2,ny2,ix,iy,iv,id,ymin,ymax,xmin,xmax,i,xx[4],yy[4];
    float *x,cos_phi,sin_phi,x1,y1,x2[2],y2[2],xx1,yy1,xx2,yy2,xc,yc,slope[2],SD,Px[2],Py[2],Pd[2],angle_det,l;

    nx2=nx/2;ny2=ny/2;
    n=nx*ny;
    SD=SO+OD;

    for(iv=0;iv<nv;iv++)
    {   x=&X[id_X[iv]*n];
        cos_phi=(float)cos(sd_phi[iv]);sin_phi=(float)sin(sd_phi[iv]);
        x1=cos_phi*(-SO);
        y1=sin_phi*(-SO);
        for(id=0;id<nd;id++)
        {   x2[0]=cos_phi*OD-sin_phi*(y_det[id]-dy_det/2);
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
}
