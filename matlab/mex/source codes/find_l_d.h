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

#define ABS(a) (a>0?a:-(a))

inline __device__ float find_l(float a,float b,float c,float x,float y)
// A method for computing the intersecting length of a pixel with a infinitely-narrow beam
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
{   float l=0,tmp,d,d0,dmax,a2,b2,tmpx,tmpy;

    a2=ABS(a);b2=ABS(b);
    tmpx=a2/2;tmpy=b2/2;
    dmax=tmpx+tmpy;
    tmp=c-a*x-b*y;
    d=ABS(tmp);

    if(d<dmax)
    {   tmp=tmpx-tmpy;
        d0=ABS(tmp);
        if(tmpx<tmpy)
        {tmp=(float)1.0/b2;}
        else
        {tmp=(float)1.0/a2;}
        if(d<=d0)
        {l=tmp;}
        else
        {l=(dmax-d)/(a2*b2);}
    }

    return l;
}
