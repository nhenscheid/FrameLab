
void Ax_fan_mf_cpu_siddon(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X);
extern "C" void Ax_fan_mf_gpu_siddon(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X,int nt);
// Siddon's method to compute the X-ray transform

void Ax_fan_mf_cpu_new(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X);
extern "C" void Ax_fan_mf_gpu_new(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X,int nt);
// A new method for computing the X-ray transform (infinitely-narrow beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).

void Ax_fan_mf_cpu_new_fb(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X,float dy_det);
extern "C" void Ax_fan_mf_gpu_new_fb(float *X,float *y,float SO,float OD,float scale,int nx,int ny,int nv,float *sd_phi,int nd,float *y_det,int *id_X,int nt,float dy_det);
// A new method for computing the X-ray transform (finite-size beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
