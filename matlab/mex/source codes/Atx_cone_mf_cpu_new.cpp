
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

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "find_l_3d.h"
#define MAX(a,b) (a>b?a:b)
#define MIN(a,b) (a<b?a:b)

void Atx_cone_mf_cpu_new(float *x,float *Y,float SO,float OD,float scale,int nx,int ny,int nz,float *sd_z,float *y_det,float *z_det,
float dy_det,float y_os,float dz_det,float dz,int nt,int *id_Y,int *Nv,int tmp_size,int nv,float *cos_phi,float *sin_phi,int na,int nb)
// A new method for computing the adjoint X-ray transform (infinitely-narrow beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
{
    int nd,nx2,ny2,nz2,na2,nb2,n,it,iv,ia,ib,id,ix,iy,iz,na_min,na_max,nb_min,nb_max,j,idx;
    float *y,xc,yc,zc,xr,yr,SD,l,tmp,x1,y1,z1,x2,y2,z2,d;

    SD=SO+OD;
    na2=na/2;nb2=nb/2;
    nx2=nx/2;ny2=ny/2;nz2=nz/2;
    n=nx*ny*nz;
    nd=na*nb;
    d=(float)sqrt((1+dz*dz)/2);

    for(iz=0;iz<nz;iz++)
    {   zc=(float)(iz+0.5-nz2)*dz;
        for(iy=0;iy<ny;iy++)
        {   yc=(float)(iy+0.5-ny2);
            for(ix=0;ix<nx;ix++)
            {   xc=(float)(ix+0.5-nx2);

                idx=iz*ny*nx+iy*nx+ix;
                for(it=0;it<nt;it++)
                {   x[it*n+idx]=0;
                    for(j=0;j<Nv[it];j++)
                    {   iv=id_Y[it*tmp_size+j];
                        y=&Y[iv*nd];

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
                            {   id=ib*na+ia;
                                x1=cos_phi[iv]*(-SO);
                                y1=sin_phi[iv]*(-SO);
                                z1=sd_z[iv];
                                x2=cos_phi[iv]*OD-sin_phi[iv]*y_det[ia];
                                y2=sin_phi[iv]*OD+cos_phi[iv]*y_det[ia];
                                z2=z_det[ib]+sd_z[iv];

                                l=find_l_3d(x1,y1,z1,x2,y2,z2,1.0,1.0,dz,xc,yc,zc);
                                x[it*n+idx]+=l*y[id];
                            }
                        }
                    }
                    x[it*n+idx]*=scale;
                }
            }
        }
    }
}

