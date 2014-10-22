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

int sort_alpha(int Nx,int Ny,float *alpha_x,float *alpha_y,float *alpha)
// Part of Siddon's algorithm: sort and bin alpha_x, alpha_y into alpha
// Here it is assumed that alpha_x (alpha_y) has the ascending order.
{   int n_alpha,ix,iy;
    float alpha_max;
    bool search_alpha_x;

    n_alpha=0;ix=0;iy=0;
    if(alpha_x[0]<alpha_y[0])
    {search_alpha_x=true;alpha_max=alpha_y[0];}
    else
    {search_alpha_x=false;alpha_max=alpha_x[0];}
    while(ix<Nx||iy<Ny)
    {   if(search_alpha_x)
        {   if(iy>=Ny)
            {   while(ix<Nx)
                {alpha[n_alpha]=alpha_x[ix];n_alpha++;ix++;}
            }
            else
            {   while(alpha_x[ix]<=alpha_max&&ix<Nx)
                {alpha[n_alpha]=alpha_x[ix];n_alpha++;ix++;}
                if(alpha_x[ix-1]==alpha_max)
                {iy++;}
                if(ix<Nx)
                {alpha_max=alpha_x[ix];}
                search_alpha_x=false;
            }
        }
        else
        {   if(ix>=Nx)
            {   while(iy<Ny)
                {alpha[n_alpha]=alpha_y[iy];n_alpha++;iy++;}
            }
            else
            {   while(alpha_y[iy]<=alpha_max&&iy<Ny)
                {alpha[n_alpha]=alpha_y[iy];n_alpha++;iy++;}
                if(alpha_y[iy-1]==alpha_max)
                {ix++;}
                if(iy<Ny)
                {alpha_max=alpha_y[iy];}
                search_alpha_x=true;
            }
        }
    }
    return n_alpha;
}
