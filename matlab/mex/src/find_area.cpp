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
#define SGN(a) (a>0?1:-1)

float find_area(float a1,float b1,float c1,float a2,float b2,float c2,float x,float y)
// A method for computing the intersecting area of a pixel with a finite-size beam
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
{   float tmp,d,d0,dmax,aa,bb,tmpx,tmpy,A1,A2,A,sgn1,sgn2;
    // step 1: A1 and sgn1
    aa=ABS(a1);bb=ABS(b1);
    tmpx=aa/2;tmpy=bb/2;
    dmax=tmpx+tmpy;
    tmp=c1-a1*x-b1*y;
    d=ABS(tmp);
    sgn1=SGN(tmp);
    if(d<dmax)
    {   tmp=tmpx-tmpy;
        d0=ABS(tmp);
        if(tmpx<tmpy)
        {tmp=(float)1.0/bb;}
        else
        {tmp=(float)1.0/aa;}
        if(d<=d0)
        {A1=d*tmp;}
        else
        {A1=(float)0.5-(float)0.5*(dmax-d)*(dmax-d)/(aa*bb);}
    }
    else
    {A1=(float)0.5;}
    // step 2: A2 and sgn2
    aa=ABS(a2);bb=ABS(b2);
    tmpx=aa/2;tmpy=bb/2;
    dmax=tmpx+tmpy;
    tmp=c2-a2*x-b2*y;
    d=ABS(tmp);
    sgn2=SGN(tmp);
    if(d<dmax)
    {   tmp=tmpx-tmpy;
        d0=ABS(tmp);
        if(tmpx<tmpy)
        {tmp=(float)1.0/bb;}
        else
        {tmp=(float)1.0/aa;}
        if(d<=d0)
        {A2=d*tmp;}
        else
        {A2=(float)0.5-(float)0.5*(dmax-d)*(dmax-d)/(aa*bb);}
    }
    else
    {A2=(float)0.5;}
    // step 3: A
    tmp=A1-sgn1*sgn2*A2;
    A=ABS(tmp);

    return A;
}
