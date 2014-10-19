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
#define ABS(a) (a>0?a:-(a))

void Ax_fan_mf_cpu_new(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X)
// A new method for computing the X-ray transform (infinitely-narrow beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
{   int n,nx2,ny2,ix,iy,iv,id,c1,c2;
    float *x,cos_phi,sin_phi,x1,y1,x2,y2,xx1,yy1,xx2,yy2,slope,l,d;

    nx2=nx/2;ny2=ny/2;
    n=nx*ny;

    for(iv=0;iv<nv;iv++)
    {   x=&X[id_X[iv]*n];
        cos_phi=(float)cos(sd_phi[iv]);sin_phi=(float)sin(sd_phi[iv]);
        x1=cos_phi*(-SO);
        y1=sin_phi*(-SO);
        for(id=0;id<nd;id++)
        {   x2=cos_phi*OD-sin_phi*y_det[id];
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

                    if(c2==c1)
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

                    if(c2==c1)
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
}
