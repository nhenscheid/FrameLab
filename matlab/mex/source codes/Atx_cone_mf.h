
void Atx_cone_mf_cpu_new(float *x,float *Y,float SO,float OD,float scale,int nx,int ny,int nz,float *sd_z,float *y_det,float *z_det,
float dy_det,float y_os,float dz_det,float dz,int nt,int *id_Y,int *Nv,int tmp_size,int nv,float *cos_phi,float *sin_phi,int na,int nb);
extern "C" void Atx_cone_mf_gpu_new(float *X,float *y,float *cos_phi,float *sin_phi,float *sd_z,float *y_det,float *z_det,int *id_Y,int *Nv,
float SO,float OD,float scale,float dy_det,float y_os,float dz_det,float dz,int nx,int ny,int nz,int nt,int na,int nb,int nv,int tmp_size,int nv_block);
// A new method for computing the adjoint X-ray transform (infinitely-narrow beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
