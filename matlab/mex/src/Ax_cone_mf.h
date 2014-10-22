
void Ax_cone_mf_cpu_new(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nz,int nv,float *sd_phi,float *sd_z,int na,int nb,float *y_det,float *z_det,int *id_X,float dz);
extern "C" void Ax_cone_mf_gpu_new(float *X,float *y,float *sd_phi,float *sd_z,float *y_det,float *z_det,int *id_Y,int *Nv,
float SO,float OD,float scale,float dz,int nx,int ny,int nz,int nt,int na,int nb,int nv,int tmp_size,int nv_block);
// A new method for computing the X-ray transform (infinitely-narrow beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).

void Ax_cone_mf_cpu_siddon(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nz,int nv,float *sd_phi,float *sd_z,int na,int nb,float *y_det,float *z_det,int *id_X,float dz0);
extern "C" void Ax_cone_mf_gpu_siddon(float *X,float *y,float *sd_phi,float *sd_z,float *y_det,float *z_det,int *id_Y,int *Nv,
float SO,float OD,float scale,float dz,int nx,int ny,int nz,int nt,int na,int nb,int nv,int tmp_size,int nv_block);
// Siddon's method to compute the X-ray transform
