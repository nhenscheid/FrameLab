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

inline __device__ float find_l(float x1_0,float y1_0,float x2_0,float y2_0,float dx,float dy,float x,float y)
{   float l=0,dx2,dy2,a,b,slope,tmp,tmp2,xi[2],yi[2],x1,y1;
    int i;

    a=x2_0-x1_0;b=y2_0-y1_0;
    dx2=dx/2;dy2=dy/2;

    if(a==0)
    {   tmp=ABS(x1_0-x);
        if(tmp<=dx2)
        {l=dy;}
    }
    else
    {   if(b==0)
        {   tmp=ABS(y1_0-y);
            if(tmp<=dy2)
            {l=dx;}
        }
        else
        {   x1=x1_0-x;y1=y1_0-y;
//            x2=x2_0-x;y2=y2_0-y;

            i=0;
            if(ABS(a)>ABS(b))
            {   slope=b/a;
                tmp=slope*(-x1)+y1;
                tmp2=slope*dx2;
                if(ABS(tmp-tmp2)<=dy2)
                {xi[i]=-dx2;yi[i]=tmp-tmp2;i++;}
                if(ABS(tmp+tmp2)<=dy2)
                {xi[i]=dx2;yi[i]=tmp+tmp2;i++;}

                if(i<2)
                {   slope=a/b;
                    tmp=slope*(-y1)+x1;
                    tmp2=slope*dy2;
                    if(ABS(tmp-tmp2)<=dx2)
                    {yi[i]=-dy2;xi[i]=tmp-tmp2;i++;}
                    if(i<2)
                    {   if(ABS(tmp+tmp2)<=dx2)
                        {yi[i]=dy2;xi[i]=tmp+tmp2;i++;}
                    }
                }
            }
            else
            {   slope=a/b;
                tmp=slope*(-y1)+x1;
                tmp2=slope*dy2;
                if(ABS(tmp-tmp2)<=dx2)
                {yi[i]=-dy2;xi[i]=tmp-tmp2;i++;}
                if(ABS(tmp+tmp2)<=dx2)
                {yi[i]=dy2;xi[i]=tmp+tmp2;i++;}

                if(i<2)
                {   slope=b/a;
                    tmp=slope*(-x1)+y1;
                    tmp2=slope*dx2;
                    if(ABS(tmp-tmp2)<=dy2)
                    {xi[i]=-dx2;yi[i]=tmp-tmp2;i++;}
                    if(i<2)
                    {   if(ABS(tmp+tmp2)<=dy2)
                        {xi[i]=dx2;yi[i]=tmp+tmp2;i++;}
                    }
                }
            }

            if(i==2)
            {   tmp=xi[1]-xi[0];tmp2=yi[1]-yi[0];
                l=(float)sqrt(tmp*tmp+tmp2*tmp2);
            }
        }
    }
    return l;
}

inline __device__ float find_l_3d(float x1_0,float y1_0,float z1_0,float x2_0,float y2_0,float z2_0,float dx,float dy,float dz,float x,float y,float z)
// assuming c~=0
// A method for computing the intersecting length of a voxel with a infinitely-narrow beam
// A better formula will be supplied to improve the speed.
{   float l=0,dx2,dy2,dz2,a,b,c,slope,tmp[2],tmp2[2],tmpx,tmpy,tmpz,xi[2],yi[2],zi[2],x1,y1,z1;
    int i;

    a=x2_0-x1_0;b=y2_0-y1_0;c=z2_0-z1_0;
    dx2=dx/2;dy2=dy/2;dz2=dz/2;

    if(a==0)
    {l=find_l(y1_0,z1_0,y2_0,z2_0,dy,dz,y,z);}
    else
    {   if(b==0)
        {l=find_l(x1_0,z1_0,x2_0,z2_0,dx,dz,x,z);}
        else
        {   x1=x1_0-x;y1=y1_0-y;z1=z1_0-z;
//            x2=x2_0-x;y2=y2_0-y;z2=z2_0-z;

            i=0;
            if(ABS(a)>ABS(b))
            {   slope=b/a;tmp[0]=slope*(-x1)+y1;tmp2[0]=slope*dx2;
                slope=c/a;tmp[1]=slope*(-x1)+z1;tmp2[1]=slope*dx2;
                if(ABS(tmp[0]-tmp2[0])<=dy2&&ABS(tmp[1]-tmp2[1])<=dz2)
                {xi[i]=-dx2;yi[i]=tmp[0]-tmp2[0];zi[i]=tmp[1]-tmp2[1];i++;}
                if(ABS(tmp[0]+tmp2[0])<=dy2&&ABS(tmp[1]+tmp2[1])<=dz2)
                {xi[i]=dx2;yi[i]=tmp[0]+tmp2[0];zi[i]=tmp[1]+tmp2[1];i++;}

                if(i<2)
                {   slope=a/b;tmp[0]=slope*(-y1)+x1;tmp2[0]=slope*dy2;
                    slope=c/b;tmp[1]=slope*(-y1)+z1;tmp2[1]=slope*dy2;
                    if(ABS(tmp[0]-tmp2[0])<=dx2&&ABS(tmp[1]-tmp2[1])<=dz2)
                    {xi[i]=tmp[0]-tmp2[0];yi[i]=-dy2;zi[i]=tmp[1]-tmp2[1];i++;}
                    if(i<2)
                    {   if(ABS(tmp[0]+tmp2[0])<=dx2&&ABS(tmp[1]+tmp2[1])<=dz2)
                        {xi[i]=tmp[0]+tmp2[0];yi[i]=dy2;zi[i]=tmp[1]+tmp2[1];i++;}
                    }
                }

                if(i<2)
                {   slope=a/c;tmp[0]=slope*(-z1)+x1;tmp2[0]=slope*dz2;
                    slope=b/c;tmp[1]=slope*(-z1)+y1;tmp2[1]=slope*dz2;
                    if(ABS(tmp[0]-tmp2[0])<=dx2&&ABS(tmp[1]-tmp2[1])<=dy2)
                    {xi[i]=tmp[0]-tmp2[0];yi[i]=tmp[1]-tmp2[1];zi[i]=-dz2;i++;}
                    if(i<2)
                    {   if(ABS(tmp[0]+tmp2[0])<=dx2&&ABS(tmp[1]+tmp2[1])<=dy2)
                        {xi[i]=tmp[0]+tmp2[0];yi[i]=tmp[1]+tmp2[1];zi[i]=dz2;i++;}
                    }
                }
            }
            else
            {   slope=a/b;tmp[0]=slope*(-y1)+x1;tmp2[0]=slope*dy2;
                slope=c/b;tmp[1]=slope*(-y1)+z1;tmp2[1]=slope*dy2;
                if(ABS(tmp[0]-tmp2[0])<=dx2&&ABS(tmp[1]-tmp2[1])<=dz2)
                {xi[i]=tmp[0]-tmp2[0];yi[i]=-dy2;zi[i]=tmp[1]-tmp2[1];i++;}
                if(ABS(tmp[0]+tmp2[0])<=dx2&&ABS(tmp[1]+tmp2[1])<=dz2)
                {xi[i]=tmp[0]+tmp2[0];yi[i]=dy2;zi[i]=tmp[1]+tmp2[1];i++;}

                if(i<2)
                {   slope=b/a;tmp[0]=slope*(-x1)+y1;tmp2[0]=slope*dx2;
                    slope=c/a;tmp[1]=slope*(-x1)+z1;tmp2[1]=slope*dx2;
                    if(ABS(tmp[0]-tmp2[0])<=dy2&&ABS(tmp[1]-tmp2[1])<=dz2)
                    {xi[i]=-dx2;yi[i]=tmp[0]-tmp2[0];zi[i]=tmp[1]-tmp2[1];i++;}
                    if(i<2)
                    {   if(ABS(tmp[0]+tmp2[0])<=dy2&&ABS(tmp[1]+tmp2[1])<=dz2)
                        {xi[i]=dx2;yi[i]=tmp[0]+tmp2[0];zi[i]=tmp[1]+tmp2[1];i++;}
                    }
                }

                if(i<2)
                {   slope=a/c;tmp[0]=slope*(-z1)+x1;tmp2[0]=slope*dz2;
                    slope=b/c;tmp[1]=slope*(-z1)+y1;tmp2[1]=slope*dz2;
                    if(ABS(tmp[0]-tmp2[0])<=dx2&&ABS(tmp[1]-tmp2[1])<=dy2)
                    {xi[i]=tmp[0]-tmp2[0];yi[i]=tmp[1]-tmp2[1];zi[i]=-dz2;i++;}
                    if(i<2)
                    {   if(ABS(tmp[0]+tmp2[0])<=dx2&&ABS(tmp[1]+tmp2[1])<=dy2)
                        {xi[i]=tmp[0]+tmp2[0];yi[i]=tmp[1]+tmp2[1];zi[i]=dz2;i++;}
                    }
                }
            }

            if(i==2)
            {   tmpx=xi[1]-xi[0];tmpy=yi[1]-yi[0];tmpz=zi[1]-zi[0];
                l=(float)sqrt(tmpx*tmpx+tmpy*tmpy+tmpz*tmpz);
            }
        }
    }
    return l;
}
