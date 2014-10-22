
void Atx_fan_mf_cpu_new(float *x,float *Y,float SO,float OD,float scale,int nx,int ny,float *sd_phi,int nd,float *y_det,float dy_det,float y_os,
int nt,int *id_Y,int *Nv,int tmp_size,int nv,float *cos_phi,float *sin_phi,float *cos_det,float *sin_det);

/*
extern "C" void Atx_fan_mf_gpu_new(float *x,float *Y,float SO,float OD,float scale,int nx,int ny,float *sd_phi,int nd,float *y_det,float dy_det,float y_os,
int nt,int *id_Y,int *Nv,int tmp_size,int nv,float *cos_phi,float *sin_phi,float *cos_det,float *sin_det);
// A new method for computing the adjoint X-ray transform (infinitely-narrow beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).


void Atx_fan_mf_cpu_new_fb(float *x,float *Y,float SO,float OD,float scale,int nx,int ny,float *sd_phi,int nd,float *y_det,float dy_det,int nt,
int *id_Y,int *Nv,int tmp_size,int nv,float *cos_phi,float *sin_phi,float *cos_det,float *sin_det);
extern "C" void Atx_fan_mf_gpu_new_fb(float *x,float *Y,float SO,float OD,float scale,int nx,int ny,float *sd_phi,int nd,float *y_det,float dy_det,int nt,
int *id_Y,int *Nv,int tmp_size,int nv,float *cos_phi,float *sin_phi,float *cos_det,float *sin_det);
// A new method for computing the adjoint X-ray transform (finite-size beam)
// The algorithm details are available in
// H. Gao. "Fast parallel algorithms for the X-ray transform and its adjoint", Medical Physics (2012).
*/