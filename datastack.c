#define STACKOK 1
#define STACKNOTOK 0
struct runavg {
  int width; //Number of samples in the stack
  int current; //Index of current sample 
  int safecount; //Count entries until the stack has been filled. 
  double avg;
  double *smpl;  //Table containing the averaged numbers. 
};
void runavg_init(struct runavg *a,int l);
int fill_running_avg(struct runavg *a,double s);
int example(void);

//============================================
//============================================
int example(void) {
  struct runavg avg;
  runavg_init(&avg,50);

  for(int i=0;i<400;i++) {
    double s=(float)random()/RAND_MAX;
    if(fill_running_avg(&avg,s)) fprintf(stdout,"%d %f\n",i,avg.avg);
  }
  return 0;
}
//============================================
//============================================
void runavg_init(struct runavg *a,int l) {
  a->width=l;
  a->current=0;
  a->safecount=0;
  a->avg=0;
  a->smpl=new double [a->width];
  for(int i=0;i<a->width;i++) a->smpl[i]=0.0;
  return;
}
//============================================
//============================================
int fill_running_avg(struct runavg *a,double s) {
  s=s/a->width;
  a->avg=a->avg+s-a->smpl[a->current];
  a->smpl[a->current]=s;  
  a->current=(a->current+1)%a->width;
  if(a->safecount<a->width) {
    a->safecount++;
    return STACKNOTOK;
  }
  else
    return STACKOK;
}
//============================================
//============================================
