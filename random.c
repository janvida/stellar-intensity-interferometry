#define RNDMTWOPI 6.283185307179586
double rndm(void);
double rndmexp(void);
int rndmpoiss(double mu);
double rndmgauss(void);
//==========================================================
//==========================================================
double rndm(void) {
  //Returns random number with uniform deviate in [0,1[
  return ((float)random())/RAND_MAX;
}
//==========================================================
//==========================================================
double rndmexp(void) {
  //Return a random number with exponential deviate of mean 1
  double u;
  do{u=rndm();} while(u==0.0||u==1); 
  return -log(u);
}
//==========================================================
//==========================================================
int rndmpoiss(double mu){
  //Returns an random integer with Poisson distribution of mean mu
  int rtn=0;
  double x=0;
  /*Count the number of random number (with exponential distribution 
    of mean 1/mu) that must be added to reach 1.0 */
  while ((x=x+rndmexp()/mu)<1.0) rtn++;
  return rtn;
}
//==========================================================
//==========================================================
double rndmgauss(void) {
  //Returns a random number of mean 0 and standar deviation 1. 
  return sqrt(-2*log(rndm()))*cos(RNDMTWOPI*rndm());
}
//==========================================================
//==========================================================
