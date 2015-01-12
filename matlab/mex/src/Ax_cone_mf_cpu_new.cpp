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
#include "mex.h"
#define ABS(a) (a>0?a:-(a))

void Ax_cone_mf_cpu_new(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nz,int nv,float *sd_phi,float *sd_z,int na,int nb,float *y_det,float *z_det,int *id_X,float dz)
// A new method for computing the X-ray transform (infinitely-narrow beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
{   int n,nd,nx2,ny2,nz2,ix,iy,iz,iv,id,ia,ib,cx1,cx2,cy1,cy2,cz1,cz2;
    float *x,cos_phi,sin_phi,x1,y1,x2,y2,z1,z2,xx1,yy1,zz1,xx2,yy2,zz2,slope1,slope2,l,d,tmp,rx,ry,rz;

    nx2=nx/2;ny2=ny/2;nz2=nz/2;
    n=nx*ny*nz;
    nd=na*nb;

    for(iv=0;iv<nv;iv++)
    {   
        mexPrintf("%s%i\n","iv = ", iv);
        mexPrintf("%s%f\n","x1 = ",x1);
        mexPrintf("%s%f\n","y1 = ",y1);
        mexPrintf("%s%f\n","z1 = ",z1);
        x=&X[id_X[iv]*n];
        cos_phi=(float)cos(sd_phi[iv]);sin_phi=(float)sin(sd_phi[iv]);
        x1=cos_phi*(-SO);
        y1=sin_phi*(-SO);
        z1=sd_z[iv];
        for(ib=0;ib<nb;ib++)
        {   for(ia=0;ia<na;ia++)
            {   x2=cos_phi*OD-sin_phi*y_det[ia];
                y2=sin_phi*OD+cos_phi*y_det[ia];
                z2=z_det[ib]+sd_z[iv];
                //mexPrintf("%s%f\n","x1 = ",x1);
                //mexPrintf("%s%f\n","y1 = ",y1);
                //mexPrintf("%s%f\n","z1 = ",z1);
                //mexPrintf("%s%f\n","x2 = ",x2);
                //mexPrintf("%s%f\n","y2 = ",y2);
                //mexPrintf("%s%f\n","z2 = ",z2);
                id=ib*na+ia;
                y[iv*nd+id]=0;
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
                                        iy=cy1;iz=cz1;y[iv*nd+id]+=l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                }
                                else
                                {   if(cz1>=-1&&cz1<=nz-1)// 12
                                    {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        rz=(cz2-zz1)/(zz2-zz1);
                                        if(cz1>-1)
                                        {iy=cy1;iz=cz1;y[iv*nd+id]+=rz*l*x[iz*ny*nx+iy*nx+ix];}
                                        if(cz2<nz)
                                        {iy=cy1;iz=cz2;y[iv*nd+id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];}
                                    }
                                }
                            }
                        }
                        else
                        {   if(cy1>=-1&&cy1<=ny-1)
                            {   if(cz2==cz1)
                                {   if(cz1>=0&&cz1<=nz-1)// 21
                                    {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        ry=(cy2-yy1)/d;
                                        if(cy1>-1)
                                        {iy=cy1;iz=cz1;y[iv*nd+id]+=ry*l*x[iz*ny*nx+iy*nx+ix];}
                                        if(cy2<ny)
                                        {iy=cy2;iz=cz1;y[iv*nd+id]+=(1-ry)*l*x[iz*ny*nx+iy*nx+ix];}
                                    }
                                }
                                else
                                {   if(cz1>=-1&&cz1<=nz-1)// 22
                                    {   d=yy2-yy1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        ry=(cy2-yy1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                        if(ry>rz)
                                        {   if(cy1>-1&&cz1>-1)
                                            {iy=cy1;iz=cz1;y[iv*nd+id]+=rz*l*x[iz*ny*nx+iy*nx+ix];}
                                            if(cy1>-1&&cz2<nz)
                                            {iy=cy1;iz=cz2;y[iv*nd+id]+=(ry-rz)*l*x[iz*ny*nx+iy*nx+ix];}
                                            if(cy2<ny&&cz2<nz)
                                            {iy=cy2;iz=cz2;y[iv*nd+id]+=(1-ry)*l*x[iz*ny*nx+iy*nx+ix];}
                                        }
                                        else
                                        {   if(cy1>-1&&cz1>-1)
                                            {iy=cy1;iz=cz1;y[iv*nd+id]+=ry*l*x[iz*ny*nx+iy*nx+ix];}
                                            if(cy2<ny&&cz1>-1)
                                            {iy=cy2;iz=cz1;y[iv*nd+id]+=(rz-ry)*l*x[iz*ny*nx+iy*nx+ix];}
                                            if(cy2<ny&&cz2<nz)
                                            {iy=cy2;iz=cz2;y[iv*nd+id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];}
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
                                        ix=cx1;iz=cz1;y[iv*nd+id]+=l*x[iz*ny*nx+iy*nx+ix];
                                    }
                                }
                                else
                                {   if(cz1>=-1&&cz1<=nz-1)// 12
                                    {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        rz=(cz2-zz1)/(zz2-zz1);
                                        if(cz1>-1)
                                        {ix=cx1;iz=cz1;y[iv*nd+id]+=rz*l*x[iz*ny*nx+iy*nx+ix];}
                                        if(cz2<nz)
                                        {ix=cx1;iz=cz2;y[iv*nd+id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];}
                                    }
                                }
                            }
                        }
                        else
                        {   if(cx1>=-1&&cx1<=nx-1)
                            {   if(cz2==cz1)
                                {   if(cz1>=0&&cz1<=nz-1)// 21
                                    {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        rx=(cx2-xx1)/d;
                                        if(cx1>-1)
                                        {ix=cx1;iz=cz1;y[iv*nd+id]+=rx*l*x[iz*ny*nx+iy*nx+ix];}
                                        if(cx2<nx)
                                        {ix=cx2;iz=cz1;y[iv*nd+id]+=(1-rx)*l*x[iz*ny*nx+iy*nx+ix];}
                                    }
                                }
                                else
                                {   if(cz1>=-1&&cz1<=nz-1)// 22
                                    {   d=xx2-xx1;tmp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);l=(float)sqrt((d*d+1)*(tmp+(z1-z2)*(z1-z2))/tmp);
                                        rx=(cx2-xx1)/d;rz=(cz2-zz1)/(zz2-zz1);
                                        if(rx>rz)
                                        {   if(cx1>-1&&cz1>-1)
                                            {ix=cx1;iz=cz1;y[iv*nd+id]+=rz*l*x[iz*ny*nx+iy*nx+ix];}
                                            if(cx1>-1&&cz2<nz)
                                            {ix=cx1;iz=cz2;y[iv*nd+id]+=(rx-rz)*l*x[iz*ny*nx+iy*nx+ix];}
                                            if(cx2<nx&&cz2<nz)
                                            {ix=cx2;iz=cz2;y[iv*nd+id]+=(1-rx)*l*x[iz*ny*nx+iy*nx+ix];}
                                        }
                                        else
                                        {   if(cx1>-1&&cz1>-1)
                                            {ix=cx1;iz=cz1;y[iv*nd+id]+=rx*l*x[iz*ny*nx+iy*nx+ix];}
                                            if(cx2<nx&&cz1>-1)
                                            {ix=cx2;iz=cz1;y[iv*nd+id]+=(rz-rx)*l*x[iz*ny*nx+iy*nx+ix];}
                                            if(cx2<nx&&cz2<nz)
                                            {ix=cx2;iz=cz2;y[iv*nd+id]+=(1-rz)*l*x[iz*ny*nx+iy*nx+ix];}
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                y[iv*nd+id]*=scale;
            }
        }
    }
}
