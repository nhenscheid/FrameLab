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
#define MAX(a,b) (a>b?a:b)
#define MIN(a,b) (a<b?a:b)

void Atx_fan_mf_cpu_new_fb(float *x,float *Y,float SO,float OD,float scale,int nx,int ny,float *sd_phi,int nd,float *y_det,float dy_det,int nt,
int *id_Y,int *Nv,int tmp_size,int nv,float *cos_phi,float *sin_phi,float *cos_det,float *sin_det)
// A new method for computing the adjoint X-ray transform (finite-size beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
{
    int nx2,ny2,nd2,n,it,iv,id,ix,iy,nd_min,nd_max,j;
    float *y,xc,yc,xr,yr,SD,Px[2],Py[2],Pd[2],l,tmp;

    nd2=nd/2;
    SD=SO+OD;
    nx2=nx/2;ny2=ny/2;
    n=nx*ny;

    for(iy=0;iy<ny;iy++)
    {   yc=(float)(iy+0.5-ny2);
        for(ix=0;ix<nx;ix++)
        {   xc=(float)(ix+0.5-nx2);
            for(it=0;it<nt;it++)
            {   x[it*n+iy*nx+ix]=0;
                for(j=0;j<Nv[it];j++)
                {   iv=id_Y[it*tmp_size+j];
                    y=&Y[iv*nd];

                    xr=cos_phi[iv]*xc+sin_phi[iv]*yc;
                    yr=-sin_phi[iv]*xc+cos_phi[iv]*yc;
                    tmp=SD/((xr+SO)*dy_det);
                    nd_max=(int)floor((yr+1)*tmp+nd2);
                    nd_min=(int)floor((yr-1)*tmp+nd2);

                    for(id=MAX(0,nd_min);id<=MIN(nd_max,nd-1);id++)
                    {   Px[0]=-(sin_det[id]*cos_phi[iv]+cos_det[id]*sin_phi[iv]);
                        Py[0]=cos_det[id]*cos_phi[iv]-sin_det[id]*sin_phi[iv];
                        Pd[0]=SO*sin_det[id];
                        Px[1]=-(sin_det[id+1]*cos_phi[iv]+cos_det[id+1]*sin_phi[iv]);
                        Py[1]=cos_det[id+1]*cos_phi[iv]-sin_det[id+1]*sin_phi[iv];
                        Pd[1]=SO*sin_det[id+1];
                        l=find_area(Px[0],Py[0],Pd[0],Px[1],Py[1],Pd[1],xc,yc);
                        if(l>0)
                        {x[it*n+iy*nx+ix]+=l*y[id];}
                    }
                }
                x[it*n+iy*nx+ix]*=scale;
            }
        }
    }
}

