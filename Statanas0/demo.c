#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "./statanas.c"
#include "./random.c"

int main(void) {
  struct h1d hist;
  hist=hcreate(5,0,1);
  for(int i=1;i<=10000;i++){
    double val=rndm_();
    h1dfill(&hist,val,1);
  }

  h1dprint(hist);


  struct h1d gauss_hist;
  gauss_hist=hcreate(20,-6,6);
  for(int i=1;i<=10000;i++){
    double val=rndmgauss();
    h1dfill(&gauss_hist,val,1);
  }
  h1dprint(gauss_hist);
  return 0;
}
