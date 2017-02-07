//==========================================================
//==========================================================
struct cmplx {
  double r;
  double i;
};
//==========================================================
//==========================================================
struct cmplx cmplx_init(double r,double i) {
  struct cmplx rtn;
  rtn.r=r;
  rtn.i=i;
  return rtn;    
}
//==========================================================
//==========================================================
struct cmplx cmplx_add(struct cmplx z1, struct cmplx z2) {
  struct cmplx rtn;
  rtn.r=z1.r+z2.r;
  rtn.i=z1.i+z2.i;
  return rtn;    
}
//==========================================================
//==========================================================
struct cmplx cmplx_prod(struct cmplx z1, struct cmplx z2) {
  struct cmplx rtn;
  rtn.r=z1.r*z2.r-z1.i*z2.i;
  rtn.i=z1.r*z2.i+z1.i*z2.r;
  return rtn;    
}
//==========================================================
//==========================================================
double cmplx_mag(struct cmplx z) {
  double rtn;
  rtn=sqrt(z.r*z.r+z.i*z.i);
  return rtn;    
}
//==========================================================
//==========================================================
struct cmplx cmplx_conjug(struct cmplx z) {
  struct cmplx rtn;
  rtn.r=z.r;
  rtn.i=-z.i;
  return rtn;    
}
//==========================================================
//==========================================================
struct cmplx cmplx_phasor(double mag,double phase) {
  struct cmplx rtn;
  rtn.r=mag*cos(phase);
  rtn.i=mag*sin(phase);
  return rtn;
}
//==========================================================
//==========================================================
