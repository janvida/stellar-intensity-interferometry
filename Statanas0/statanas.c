#define NSIM 1000  //This must desappear
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr 

struct h1d {
  float min;
  float max;
  int nbin;
  float underflow;
  float overflow;
  float *c;
  float cmin;
  float cmax;
  float sum;
  int nentries;
  //Could add average and standard deviation
};
float LiandMa(int on,int off,float beta);
struct h1d hcreate(int n,float min,float max);
void h1dreset(struct h1d *h);
void h1dfill(struct h1d *h,float x,float w);
void h1dprint(struct h1d h);
void h1dadd(struct h1d *h1,struct h1d h2,float c1,float c2);
void autocor(struct h1d histoplnt,struct h1d histosim);
void fourier(struct h1d histoplnt,struct h1d histosim);
float h1autocor(struct h1d hi,int tau);
float lagrange(float x, struct h1d hi);
float rndm_(void);
float deviate(struct h1d hi);
void fourn(float data[], unsigned long nn[], int ndim, int isign);
void h1four(struct h1d h,struct h1d *rtn);
//==============================================================
//==============================================================
float LiandMa(int on,int off,float beta){
  // return significance calculated with the Li and Ma formula
  double betasq,oneplusbeta,oneplusbetaoverbeta,Ntot;
  betasq=beta*beta;
  oneplusbeta=1.0+beta;
  oneplusbetaoverbeta=oneplusbeta/beta;
  //if(on*off<1) return 0.;
  if(beta==0.) return 0;
  double Non=double(on);
  double Noff=double(off);
  double Nsig=0;
  Ntot=Non+Noff;
  Nsig   = Non - beta*Noff;
  double sig;
  if(Ntot == 0.) sig=0.;
  else if( Non == 0 && Noff != 0. )
    {
      sig = sqrt(2.*( Noff*log(oneplusbeta*(Noff/Ntot)) ) );
    }
  else if( Non != 0 && Noff == 0. )
    {
      sig = sqrt(2.*( Non *log(oneplusbetaoverbeta*(Non/Ntot)) ) );
    }
  else
    {
      sig = sqrt(2.*( Non *log(oneplusbetaoverbeta*(Non/Ntot)) + 
		      Noff*log(oneplusbeta*(Noff/Ntot)) ));
    }
  if( Nsig < 0 ) sig=-sig;
  return sig;
}
//==============================================================
//==============================================================
struct h1d hcreate(int n,float min,float max) {
  struct h1d rtn;
  if(n<0) {
    fprintf(stderr,"In hcreate the number of bins must be strictly positive\n");
    exit(0);
  }
  if(max<=min) {
    fprintf(stderr,"In hcreate, max must be strictly greater than min\n");
    exit(0);
  }
  rtn.max=max;
  rtn.min=min;
  rtn.nbin=n;
  rtn.underflow=0;
  rtn.overflow=0;
  rtn.nentries=0;
  rtn.sum=0;
  rtn.c=new float[n];
  for(int i=0;i<n;i++) rtn.c[i]=0;
  rtn.cmin=0;
  rtn.cmax=0;
  return rtn;
}
//==============================================================
//==============================================================
void h1dadd(struct h1d *h1,struct h1d h2,float c1,float c2){
  for(int i=0;i<h1->nbin && i<h2.nbin; i++) 
    h1->c[i]=c1*h1->c[i]+c2*h2.c[i];
  return;
}
//==============================================================
//==============================================================
void h1dreset(struct h1d *h){
    for(int i=0;i<h->nbin;i++) h->c[i]=0;
    h->cmin=0;
    h->cmax=0;
  return;
}
//==============================================================
//==============================================================
void h1dfill(struct h1d *h,float x,float w){
  h->nentries++;
  h->sum+=w;
  if(x<h->min) {h->underflow=h->underflow+w;return;}
  if(x>h->max) {h->overflow=h->overflow+w;return;}
  int k=int(floor(h->nbin*(x-h->min)/(h->max-h->min)));
  h->c[k]=h->c[k]+w;
  if(h->c[k]>h->cmax) h->cmax=h->c[k];
  if(h->c[k]<h->cmin) h->cmin=h->c[k];

  return;
}
//==============================================================
//==============================================================
void h1dprint(struct h1d h) {
  fprintf(stdout,"MIN=%f\nMAX=%f\n",h.min,h.max);
  fprintf(stdout,"NBIN=%d\n",h.nbin);
  for(int i=0;i<h.nbin;i++)
    fprintf(stdout, "%d %f\n",i+1,h.c[i]);
  fprintf(stdout,"UNDERFLOW=%f\nOVERFLOW=%f\n",h.underflow,h.overflow);
  fprintf(stdout,"ENTRIES:%d\nSUM=%f\n",h.nentries,h.sum);
  return;
}
//==============================================================
//==============================================================
void autocor(struct h1d histoplnt,struct h1d histosim){
  
  //Autocorrelation reference
  float acorzero=h1autocor(histoplnt,0);
  /*  for(int i=0;i<histoplnt.nbin;i++)
    fprintf(stdout,"%d %f\n",i,h1autocor(histoplnt,i)/acorzero);
  */


  //storage for the average correlation
  float *mbin;
  mbin= new float[histoplnt.nbin];
  for(int i=0;i<histoplnt.nbin;i++) mbin[i]=0;
  //storage for the average square correlation
  float *mbin2;
  mbin2= new float[histoplnt.nbin];
  for(int i=0;i<histoplnt.nbin;i++) mbin2[i]=0;
  //Start the simulations
  struct h1d himc;
  himc=hcreate(histoplnt.nbin,histoplnt.min,histoplnt.max);
  // struct h1d avg;
  // avg=hcreate(histoplnt.nbin,histoplnt.min,histoplnt.max);
  for(int k=0;k<NSIM;k++) {
    h1dreset(&himc);
    for(int j=0;j<histoplnt.nentries;j++)
      h1dfill(&himc,deviate(histosim),1.0);
    // h1dadd(&avg,himc,1.0,1.0/NSIM);
    float rescale=h1autocor(himc,0);
    for(int j=0;j<histoplnt.nbin;j++) {
	float dummy=h1autocor(himc,j)/rescale;
	mbin[j]=mbin[j]+dummy;
	mbin2[j]=mbin2[j]+dummy*dummy;
    }
  }
  for(int j=0;j<histoplnt.nbin;j++) {
    mbin[j]=mbin[j]/NSIM;
    mbin2[j]=mbin2[j]/NSIM;
    fprintf(stdout,"%d 0.5 %f %f\n",j,h1autocor(histoplnt,j)/acorzero,sqrt(mbin2[j]-mbin[j]*mbin[j]));
  }
  //  h1dprint(avg);
  //  h1dprint(histoplnt);
  return;
}
//==============================================================
//==============================================================
float h1autocor(struct h1d hi,int tau) {
  float rtn=0;
  for(int i=0;i<hi.nbin-tau;i++) 
    rtn=rtn+hi.c[i]*hi.c[i+tau];
  return rtn;  
}
//==============================================================
//==============================================================
float rndm_(void) {
  return float(random())/RAND_MAX;
}
//==============================================================
//==============================================================
float lagrange(float x, struct h1d hi) {
  float rtn=0;
  for(int i=0;i<hi.nbin;i++) {
    float t=hi.c[i];
    float xi=hi.min+(i+0.5)*(hi.max-hi.min)/hi.nbin; //i-th bin center
    for(int k=0;k<hi.nbin;k++){
      if(k==i) continue;
      float xk=hi.min+(k+0.5)*(hi.max-hi.min)/hi.nbin; //k-th bin center
      if (hi.c[i]-hi.c[k]!=0)
	t=t*(x-hi.c[k])/(hi.c[i]-hi.c[k]);
    }
    rtn=rtn+t;
  }
  return rtn;
}
//==============================================================
//==============================================================
float deviate(struct h1d hi) {
  float rtn;
  int dummy;
  do{
    rtn=hi.min+rndm_()*(hi.max-hi.min);
    dummy=int(floor(hi.nbin*(rtn-hi.min)/(hi.max-hi.min)));
  } while(hi.c[dummy]/hi.cmax<rndm_());

    //  } while(lagrange(rtn,hi)/hi.cmax<rndm());
  return rtn;
}
//==============================================================
//==============================================================
void fourier(struct h1d histoplnt,struct h1d histosim){

  struct h1d plntfour;
  plntfour.nbin=histoplnt.nbin/2+1;
  plntfour.c=new float[histoplnt.nbin/2+1];
  plntfour.cmin=0;
  plntfour.cmax=0.5*histoplnt.nbin/(histoplnt.max-histoplnt.min);
  struct h1d plntsimfour;
  plntsimfour.nbin=histoplnt.nbin/2+1;
  plntsimfour.c=new float[histoplnt.nbin/2+1];
  plntsimfour.cmin=0;
  plntsimfour.cmax=0.5*histoplnt.nbin/(histoplnt.max-histoplnt.min);

  h1four(histoplnt,&plntfour);

  //Start the simulations
  struct h1d simfouravg;
  simfouravg.nbin=histoplnt.nbin/2+1;
  simfouravg.c=new float[histoplnt.nbin/2+1];
  simfouravg.cmin=0;
  simfouravg.cmax=0.5*histoplnt.nbin/(histoplnt.max-histoplnt.min);
  struct h1d simfour2avg;
  simfour2avg.nbin=histoplnt.nbin/2+1;
  simfour2avg.c=new float[histoplnt.nbin/2+1];
  simfour2avg.cmin=0;
  simfour2avg.cmax=0.5*histoplnt.nbin/(histoplnt.max-histoplnt.min);
  for(int i=0;i<plntfour.nbin;i++) {
    simfouravg.c[i]=0.0;    
    simfour2avg.c[i]=0.0;
  }
  struct h1d himc;
  himc=hcreate(histoplnt.nbin,histoplnt.min,histoplnt.max);
  //  for(int k=0;k<NSIM;k++) {
  for(int k=0;k<NSIM;k++) {
    h1dreset(&himc);
    for(int j=0;j<histoplnt.nentries;j++)
      h1dfill(&himc,deviate(histosim),1.0);
    //Take the fourier transform of the simulation
    h1four(himc,&plntsimfour);
    //Calculate the variance
    for(int i=0;i<plntfour.nbin;i++) {
      simfouravg.c[i]=simfouravg.c[i]+plntsimfour.c[i];    
      simfour2avg.c[i]=simfour2avg.c[i]+plntsimfour.c[i]*plntsimfour.c[i];
    }
  }
  
  for(int i=0;i<plntfour.nbin;i++) {
    simfouravg.c[i]=simfouravg.c[i]/NSIM;
    simfour2avg.c[i]=simfour2avg.c[i]/NSIM;
    fprintf(stderr,"%d %f %g %g\n",i,i*plntfour.cmax/(plntfour.nbin-1),
	    plntfour.c[i],
	    sqrt(simfour2avg.c[i]-simfouravg.c[i]*simfouravg.c[i]));
  }
  
  return;
}
//==============================================================
void h1four(struct h1d h,struct h1d *rtn) {
  //Returns an array of length h.nbn/2+1 containing the fourier power from
  //frequency 0 to Nyquist frequency 1/(2*dx). Point j of the fourier power 
  //table forresponds to frequency f[j]=j/(nbin*dx) with nbin the number of 
  //points in the original data. The fourier power is normalized to the 
  //fourier power of frequency 0.
  //  struct h1d rtn;

  //Create an array to be used to compute the fourier transform
  float *d;
  d= new float[2*h.nbin];
  for(int i=0;i<h.nbin;i++) {
    d[2*i]=h.c[i]; //Real part
    d[2*i+1]=0.0;          //Imaginary part
  }

  //Call the FFT algorithm
  unsigned long ndim[1];
  ndim[0]=h.nbin;
  fourn(d,ndim,1,1);

  for(int i=0;i<h.nbin/2+1;i++) {
    rtn->c[i]=d[2*i]*d[2*i]+d[2*i+1]*d[2*i+1];  
  }
  for(int i=h.nbin/2;i>=0;i--) {
    rtn->c[i]=rtn->c[i]/rtn->c[0];  
  }

  return;
}
//==============================================================
//======================================================================
void fourn(float data[], unsigned long nn[], int ndim, int isign) {
  /*
    Replaces data by its ndim-dimensional discrete Fourier transform, 
    if isign is input as 1. 
    nn[1..ndim] is an integer array containing the lengths of each dimension 
    (number of complex values), which MUST all be powers of 2. 
    data is a real array of length twice the product of these lengths, 
    in which the data are stored as in a multidimensional complex array: 
    real and imaginary parts of each element are in consecutive locations, 
    and the rightmost index of the array increases most rapidly as one 
    proceeds along data. 
    For a two-dimensional array, this is equivalent to storing the array by 
    rows. 
    If isign is input as −1, data is replaced by its inverse transform times 
    the product of the lengths of all dimensions.
*/


  int idim;
  unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
  float tempi,tempr;
  //Double precision for trigonometric references
  double theta,wi,wpi,wpr,wr,wtemp; 
  //Compute the total number of complexes values
  for(ntot=1,idim=0;idim<ndim;idim++) 
    ntot*=nn[idim];

    
  nprev=1;
  //Main loop over the dimensions
  for (idim=ndim-1;idim>=0;idim--) { 
    n=nn[idim];
    nrem=ntot/(n*nprev);
    ip1=nprev << 1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev=0;
    //This is the bit-reversal section of the routine.
    for (i2=0;i2<ip2;i2+=ip1) { 
      if (i2 < i2rev) { 
	for (i1=i2;i1<i2+ip1-1;i1+=2) {
	  for (i3=i1;i3<ip3;i3+=ip2) {
	    i3rev=i2rev+i3-i2;
	    SWAP(data[i3],data[i3rev]);
	    SWAP(data[i3+1],data[i3rev+1]);
	  }
	}
      }
      ibit=ip2 >> 1;
      while (ibit >= ip1 && i2rev+1 > ibit) {
	i2rev -= ibit;
	ibit >>= 1;
      }
      i2rev += ibit;
    }
    //Here begins the Danielson-Lanczos section of the routine.
    ifp1=ip1; 
    while(ifp1 < ip2) { 
      ifp2=ifp1 << 1;
      //Initialize for the trig. recurrence.
      theta=isign*6.28318530717959/(ifp2/ip1); 
      wtemp=sin(0.5*theta); 
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (i3=0;i3<ifp1;i3+=ip1) {
	for (i1=i3;i1<i3+ip1-1;i1+=2) {
	  for (i2=i1;i2<ip3;i2+=ifp2) {
	    //Danielson-Lanczos formula:
	    k1=i2; 
	    k2=k1+ifp1;
	    tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
	    tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
	    data[k2]=data[k1]-tempr;
	    data[k2+1]=data[k1+1]-tempi;
	    data[k1] += tempr;
	    data[k1+1] += tempi;
	  }
	}
	//Trigonometric recurrence.
	wr=(wtemp=wr)*wpr-wi*wpi+wr; 
	wi=wi*wpr+wtemp*wpi+wi;
      }
      ifp1=ifp2;
    }
    nprev *= n;
  }
}
//======================================================================
//======================================================================
