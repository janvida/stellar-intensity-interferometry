#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "./cmplx.c"
#include "./random.c"
#include "./linkedlist.c"
#include "./datastack.c"
#include "./misc.c"

#define TWOPI 6.283185307179586
#define SPEED_OF_LIGHT 3E8 //Units of m/s

#define NBPOINTS 100 //Number of point sources which make up the star
#define NBTELS 2 //Number of telescopes in array

double *psignal; 

//-------Define Structures--------------------
struct star_indiv {
  double flux; //Flux of each individual point source
  double amp; //Amplitude of electric field, GAUSS
  double freq; //Frequency of light in Hz
  double phi; //Phase, chosen between 0,2Pi rad
  double thetax; //Defines x position of source on star, chosen between 0 and the size of the star in rad
  double thetay; //Defines y position of source on star, chosen between 0 and the size of the star in rad
  double rad; //Maximum angular radius of the star in radians
};
struct source {
  int npts; //Number of point sources which make up star
  double flux; //Total flux, photon*s^-1*m^-2
  double meanfreq; //Mean frequency of star light
  double rad;
  struct star_indiv *star;
};
struct tel_indiv {
  double xpos; //Position of telescope on the ground in m
  double ypos; //Position of telescope on the ground in m
  double r; //Radius of telescope in m
  double dr; //Step size to loop over area in m^2
  double a; //Area of telescope in m^2
  double mean;
  int nphot;
  double pmtsgnl; //PMT signal
  int u; //Number of step sizes taken when integrating over telescope
};
struct array {
  int ntel; //Number of telescopes in array
  double tr; //Radius of telescopes (m)
  double stepsz; //Step size for integrating over telescope area
  struct tel_indiv *tel;  
};
struct parameters {
  double dt; //Time steps in units of s
  double T; //Total time to run program in s
  double pdur; //Duration of main pulse (ns)
  double pulsedur; //Duration of entire pulse (ns)
  double pampl; //Amplitude of main pulse
  double bglight; //Amount of stray background light as a multiple of source flulx
  double exnoise; //Amount of excess noise present in system
  double maxbl; //Maximum baseline (m)
  double blstep; //Step for obtaining points along baseline (m)
  double cohertime; //Coherence time of light in given bandwidth
};
struct statistics {
  int n;
  double sum;
  double sum2;
  double mean;
  double stdev;
};
struct run {
  int exnonoff; //Turns loop on excess noise on or off
  int bslonoff; //Turns loop on baseline on or off
  char pilotfile[MAXFILENAMELENGTH];
};
 


int rndmlight(struct source *src,struct parameters p);
double pulse(double t,struct parameters p);
double test_pulse(double t,struct parameters p);
int photons(int dti,struct array *ta,struct source src,struct parameters par,int clr);
int pulseinit(struct parameters p,float a, float b, float rct, float tstep);
int source_init(struct source *src,struct parameters p,int np);
int array_init(struct array *ta,int nt);
int par_init(struct parameters *p);
int run_init(struct run *r);
int stat_init(struct statistics *stat);
int statistics(struct statistics *stat,double x);
double gauss(double x, double mean, double sdev);
int readpilot(struct parameters *p,struct run *r,struct source *src,struct array *ta);




int main(int argc,char **argv) {
 

  //Initialize loops
  struct run loop;
  run_init(&loop);

  //Declare structures
  struct parameters par;
  struct source src;
  struct array ta;

  for(int i=1;i<argc;i++){
    if(strcmp(argv[i],"-h")==0) {
      showerror("-h gives this help page");
      showerror("-p to specify a pilot file other than default ./default.pilot");
      exit(0);
    }
    if(strcmp(argv[i],"-p")==0){
      if(i+1<argc){
	i=i+1;
	strcpy(loop.pilotfile,argv[i]);
      } else {
	showxerror("Argument -p must be followed by a pilot file name");
      }
    }
  }
  
  if(readpilot(&par,&loop,&src,&ta)) showxerror("Error reading pilot file");
  

  //Initialize source
  source_init(&src,par,NBPOINTS);


  //Initialize telescopes
  array_init(&ta,NBTELS);


  //Mean (star) photon flux on each telescope
  //This is the mean number of photons from the star only, 
  //does not include background light
  for(int tel=0;tel<ta.ntel;tel++) {
    ta.tel[tel].mean=src.flux*par.dt*ta.tel[tel].a;
  }
  
  
  //Print out parameters for each file
  //fprintf(stdout,"#Number Star Points=%d\n#Integration time=%.1e\n#Timestep=%.1e\n#Coherence Time=%.1e\n#Frequency=%.2e\n#Flux=%.1e\n#Star Radius=%.1e\n#Telescope Radius=%g\n#Excess Noise=%g\n#Pulse Duration=%d\n",NBPOINTS,INTIME,TIMESTEP,COHERTIME,FREQUENCY,src.flux,SRAD,TELRAD,EXNOISE,PDUR);
  
  
  
  /*int nexpsl;
  //double *bglight_array;
  //bglight_array=new double[loop.numbg];
  double bglight_array[loop.numbg]={0.0,0.3,0.5,1.0,2.0,5.0,10.0};
  if(loop.bglonoff!=0) nexpsl=loop.numbg-1;
  else nexpsl=0;
  for(int m=0;m<=nexpsl;m++) {
  //double bglight=BGLIGHT;
  if(loop.bglonoff!=0){
  par.bglight=bglight_array[m];
      }*/
  
  // double exnoise;    
  int nexnoise;
  if (loop.exnonoff==0) {
    nexnoise=0; 
  }
  else nexnoise=10;
  for(int a=0;a<=nexnoise;a++) {
    if (nexnoise!=0) par.exnoise=a*0.1;
    
    int nbl;
    if(loop.bslonoff==0){
      for(int tel=0;tel<ta.ntel;tel++){
	ta.tel[tel].ypos=0;
      }
      nbl=0;
    }
    else nbl=floor(par.maxbl/par.blstep);
    for(int bl=0;bl<=nbl;bl++){
      ta.tel[1].ypos=bl*par.blstep;
      
      
      
      //Initialize statistics
      struct statistics photstat;
      stat_init(&photstat);
      struct statistics sgnlstat;
      stat_init(&sgnlstat);
      
      
      //Initialize pulse
      int pulsedur=par.pulsedur;
      psignal=new double[pulsedur];
      //Define pulse shape
      pulseinit(par,1.3,2.5,85.0,1.0); //PDUR=20.0 PULSEDUR=200
      
      
      //Get the mean pulse amplitude (truncation of amplitude fluctuation effect)
      double meanamp=par.pampl;
      if(par.exnoise>0.0) {
	meanamp=0;
	double sum1=0;
	double epsilon=par.exnoise/20.0;
	for(double s1=0;s1<=1+7*par.exnoise;s1=s1+epsilon) {
	  double w1=epsilon*gauss(s1,1.0,par.exnoise);
	  sum1=sum1+w1;
	  meanamp=meanamp+s1*w1;
	}
	
	meanamp=meanamp/sum1;
      }
      //fprintf(stderr,"meanamp=%f\n",meanamp);
      
      //Integral of pulses
      double intpulse=0; //Integral of pulse
      double intpulse2=0; //Integral of pulse squared
      for(double dti=0;dti<par.pulsedur;dti++) {
	double pls=pulse(dti,par);
	intpulse=intpulse+pls;
	intpulse2=intpulse2+pls*pls;
      }  
      intpulse=intpulse*meanamp;
      intpulse2=intpulse2*meanamp*meanamp;
      //Just to check that AC coupling works, the integral should be zero and the square of the integral should not be zero
      //fprintf(stderr,"intpulse=%f intpulse2=%f\n",intpulse,intpulse2);
      
      
      
      
      
      //Loop to calculate statistics on |g|^2
      for(int u=0;u<=50;u++) {
	//Reset photon list at each experiment
	photons(0,&ta,src,par,1);
	
	//Construct a past to signals (duration of one pulselength)
	for(int dti=-par.pulsedur;dti<0;dti++) {
	  rndmlight(&src,par);
	  photons(dti,&ta,src,par,0);
	}
	
	//Store both photon correlation and signal correlation
	double alpha=sqrt(par.cohertime/par.dt);
	double crphot=0;
	double crsgnl=0;
	//For normalization, the mean number of photons incident on telescope 
	//(i.e. from the star and from the background contamination)
	double norm=src.flux*src.flux*ta.tel[0].a*ta.tel[1].a*par.dt*par.T*(1+par.bglight)*(1+par.bglight);
	
	//Loop over time
	int nsteps=par.T/par.dt;
	for(int t=0;t<nsteps;t++) {
	  
	  //Reset star at each timestep
	  rndmlight(&src,par);
	  
	  //Calculate number of photons at each timestep
	  photons(t,&ta,src,par,0);
	  
	  
	  //Calculate photon correlation 
	  crphot=crphot+(ta.tel[0].nphot-ta.tel[0].mean*(1+par.bglight))*(ta.tel[1].nphot-ta.tel[1].mean*(1+par.bglight))/norm;
	  //Signal correlation calculated from pulse and subtracting average
	  //crsgnl=crsgnl+(ta.tel[0].pmtsgnl-intpulse*ta.tel[0].mean*(1+bglight))*(ta.tel[1].pmtsgnl-intpulse*ta.tel[1].mean*(1+bglight))/(intpulse2*norm);
	  
	  //Signal correlation calculated from AC coupled signal
	  crsgnl=crsgnl+(ta.tel[0].pmtsgnl*ta.tel[1].pmtsgnl)/(intpulse2*norm);
	  
	  
	}//End of loop on time
	
	crphot=crphot*(((1+par.bglight)*(1+par.bglight))/(alpha*alpha));
	crsgnl=crsgnl*(((1+par.bglight)*(1+par.bglight))/(alpha*alpha));
	
	//Calculate statistics for all experiments
	statistics(&photstat,crphot);
	statistics(&sgnlstat,crsgnl);
	
	
      }//End of u loop
      
      double param=0;
      //if(SLONOFF!=0) param=par.bglight;
      if(loop.exnonoff!=0) param=par.exnoise;
      if(loop.bslonoff!=0) param=ta.tel[1].ypos;
      fprintf(stderr,"%g %f %f %f %f\n",param,photstat.mean,photstat.stdev/sqrt(photstat.n+1),sgnlstat.mean,sgnlstat.stdev/sqrt(sgnlstat.n+1));
      fprintf(stdout,"%g %f %f %f %f\n",param,photstat.mean,photstat.stdev/sqrt(photstat.n+1),sgnlstat.mean,sgnlstat.stdev/sqrt(sgnlstat.n+1));
      double simSNR=1/photstat.stdev;
      double hbtSNR=ta.tel[0].a*src.flux*par.cohertime*(1/(1+par.bglight))*sqrt((1/par.dt)*par.T);
      fprintf(stderr,"SNR Sim %f +/- %f HBT %f\n",simSNR,simSNR*(1/sqrt(2*49)),hbtSNR);
      
    } //End of loop on BASELINE
    
  } //End of loop on EXNOISE
  
  
    //}//End of loop on STRAY LIGHT CONTAMINATION
  
  
  
  return 0;
} 

//==========================================================
//==========================================================
int rndmlight(struct source *src,struct parameters p) {
  for(int s=0;s<src->npts;s++) {
    double in=rndmexp()*src->star[s].flux;
    src->star[s].amp=sqrt(in);
    src->star[s].phi=rndm()*TWOPI;
    //Randomize the rest of the star as well
    src->star[s].freq=src->meanfreq-(rndm()*(1/p.cohertime)-0.5*(1/p.cohertime));
    double sr_ang=rndm()*TWOPI;
    double dr=src->star[s].rad*sqrt(rndm());
    src->star[s].thetax=dr*cos(sr_ang); //Set x position on the source randomly
    src->star[s].thetay=dr*sin(sr_ang); //Set y position on the source randomly

  }
  return 0;
}
//==========================================================
//==========================================================
int pulseinit(struct parameters p,float a, float b, float rct, float tstep) {
  a=a/tstep; //Related to rise time
  b=b/tstep; //Fall time
  rct=rct/tstep;
  //Define an array for tail
  int tlen=p.pulsedur-p.pdur;
  double *tail; tail=new double [tlen];
  double beta=1/(rct*(1-exp(-tlen/rct)));
  double tailint=0;
  for(int i=0;i<tlen-1;i++) {
    tail[i]=beta*exp(-i/rct);
    //printf("%d %f\n",i,tail[i]);
    tailint=tailint+tail[i];
  }
  tail[tlen-1]=1.0-tailint;
  //printf("tailint=%f\n",tailint+tail[tlen-1]);

  double signalint=0;
  
  //int pulsedur=PULSEDUR;
  //int pdur=PDUR;

  for(int i=0;i<p.pulsedur;i++) psignal[i]=0;
  
  for(int i=0;i<p.pulsedur;i++) {
    double t=0;
    if(i<p.pdur) {
      for(int j=0;j<=i;j++) t=t-tail[i-j]*(1-exp(-j/a))*exp(-j/b);
      psignal[i]=t+(1-exp(-i/a))*exp(-i/b);
    }
    else {
      for(int j=0;j<p.pdur;j++) {
	double s=-tail[i-j]*(1-exp(-j/a))*exp(-j/b);
	psignal[i]=psignal[i]+s;
      }
    }
    
    //Use this line to check what the signal looks like
    //printf("%d %f 0\n",i,psignal[i]);

    //signalint=signalint+psignal[i];
  }
  //fprintf(stderr,"signalint=%f\n",signalint);
  //exit(0);
    

  return 0;
}
//==========================================================
//==========================================================
double pulse(double t,struct parameters p) {
  //return psignal[(int)t];
  return test_pulse(t,p);
}
//==========================================================
//==========================================================
double test_pulse(double t,struct parameters p) {
  double rtn=0;
  if(t>p.pulsedur) return 0;
  if(t>=0 && t<p.pdur) return p.pampl; //Just for testing
  rtn=-p.pdur/(p.pulsedur-p.pdur);
  return rtn;
}
//==========================================================
//==========================================================
int photons(int dti,struct array *ta,struct source src,struct parameters par,int clr) {
  //Photon streams
  static struct list pstream[NBTELS];
  //Clear the stream if requested
  if (clr==1) {
    for(int tel=0;tel<NBTELS;tel++) {
      clearlist(&(pstream[tel]));
    }
    return 0;
  }

  //Calculate number of photons incident on telescope
  for(int tel=0;tel<ta->ntel;tel++) {
    double power=0;

    //Loop of area of detector


    for(double ix=-ta->tel[tel].r;ix<=ta->tel[tel].r;ix=ix+ta->tel[tel].dr){
      double x=ta->tel[tel].xpos+ix;
      for(double iy=-ta->tel[tel].r;iy<=ta->tel[tel].r;iy=iy+ta->tel[tel].dr){
	double y=ta->tel[tel].ypos+iy;
	
	double r2=ix*ix+iy*iy;
	if(r2>(ta->tel[tel].r*ta->tel[tel].r)) continue;


	//Area calculation works correctly now


	//Sum total electric field for each detector
	struct cmplx wv;
	//Initialize wave to zero for each area element of detector
	wv=cmplx_init(0.0,0.0);
	for(int s=0;s<src.npts;s++) {
	  double ps=src.star[s].phi+TWOPI*src.star[s].freq
	    *(x*src.star[s].thetax+y*src.star[s].thetay)/SPEED_OF_LIGHT;
	  //Sum of electric fields from each point on star
	  wv=cmplx_add(wv,cmplx_phasor(src.star[s].amp,ps));
	}
	//After summing the electric fields then get intensity
	double intens=cmplx_mag(wv)*cmplx_mag(wv);
	//Sum of intensities for entire telescope
	power=power+intens*ta->tel[tel].dr*ta->tel[tel].dr;
	
      }//End of loop on iy
    }//End of loop on ix

    //Pick number of photons from Poisson
    if(ta->tel[tel].r==0) ta->tel[tel].a=1;
    double alpha=sqrt(par.cohertime/par.dt);
    double mu=power*par.dt;
    double poissmu=ta->tel[tel].mean*(1-alpha)+mu*alpha;
    ta->tel[tel].nphot=rndmpoiss(poissmu+par.bglight*ta->tel[tel].mean);


    //Generate photons in linked list
    /*Create one node per photon depending on the number of photons that 
      are created. All nodes (photons) occur at the same timestep*/
    for(int phot=0;phot<ta->tel[tel].nphot;phot++) {
      makenode(&(pstream[tel]));
      pstream[tel].current->data.time=dti;
      pstream[tel].current->data.ampl=-1.0;
      while(pstream[tel].current->data.ampl<0.0)
	pstream[tel].current->data.ampl=1.0+par.exnoise*rndmgauss();
    }

    
    //Calculate signals
    ta->tel[tel].pmtsgnl=0.0;
    if(pstream[tel].nbnodes>0) {
      reset_current(&(pstream[tel]));
      do{
	double tmpt=dti-pstream[tel].current->data.time;
	if (tmpt>par.pulsedur) delete_current(&(pstream[tel]));
	else if (tmpt>=0) ta->tel[tel].pmtsgnl=ta->tel[tel].pmtsgnl
			    +pstream[tel].current->data.ampl*pulse(tmpt,par);
      } while(next_one(&(pstream[tel]))==0);
    }//End if loop
    
  }//End of loop on telescopes
  

  return 0;
}
//==========================================================
//==========================================================
int source_init(struct source *src,struct parameters p,int np) {
  src->npts=np;
  src->star=new struct star_indiv [src->npts];
  for(int s=0;s<src->npts;s++) {
    src->star[s].rad=src->rad;
    src->star[s].flux=(src->flux/src->npts);
    double in=rndmexp()*src->star[s].flux;
    src->star[s].amp=sqrt(in);
    src->star[s].freq=src->meanfreq-(rndm()*(1/p.cohertime)-0.5*(1/p.cohertime));//Pick frequency randomly
    src->star[s].phi=rndm()*TWOPI;
    double sr_ang=rndm()*TWOPI;
    double dr=src->star[s].rad*sqrt(rndm());
    src->star[s].thetax=dr*cos(sr_ang); //Set x position on the source randomly
    src->star[s].thetay=dr*sin(sr_ang); //Set y position on the source randomly
  }

  return 0;
}
//==========================================================
//==========================================================
int array_init(struct array *ta,int nt) {
  ta->ntel=nt;
  ta->tel=new struct tel_indiv [ta->ntel];
  for(int tel=0;tel<ta->ntel;tel++) {
    ta->tel[tel].pmtsgnl=0;
    ta->tel[tel].xpos=0; //x position of telescope on ground
    ta->tel[tel].ypos=0; //y position of telescope on ground
    ta->tel[tel].r=ta->tr; //Radius of telescope
    ta->tel[tel].a=0;
    ta->tel[tel].u=0;

    if(ta->tel[tel].r==0) ta->tel[tel].dr=1.0;
    else ta->tel[tel].dr=ta->tel[tel].r/ta->stepsz; //Step size 

    //Calculate total area of detectors
    //Loop of area of detector


    for(double ix=-ta->tel[tel].r;ix<=ta->tel[tel].r;ix=ix+ta->tel[tel].dr){
      double x=ta->tel[tel].xpos+ix;
      for(double iy=-ta->tel[tel].r;iy<=ta->tel[tel].r;iy=iy+ta->tel[tel].dr){
	double y=ta->tel[tel].ypos+iy;
	
	double r2=ix*ix+iy*iy;
	if(r2>(ta->tel[tel].r*ta->tel[tel].r)) continue;
	//Count total area of telescope
	ta->tel[tel].a=ta->tel[tel].a+ta->tel[tel].dr*ta->tel[tel].dr;
	ta->tel[tel].u++;
      }//End loop iy
    }//End loop ix


    //Area calculation works correctly now


    //fprintf(stderr,"TEL %d NSTEPS %d Area %f\n",tel,ta->tel[tel].u,ta->tel[tel].a);
  }//End of loop on number of detectors
}
//==========================================================
//==========================================================
int par_init(struct parameters *p) {
  p->dt=0; //Time step
  p->T=0; //Total integration time
  p->pdur=0;
  return 0;
}
//==========================================================
//==========================================================
int run_init(struct run *r) {
  r->exnonoff=0; 
  r->bslonoff=0;
  strcpy(r->pilotfile,"./default.pilot");
  return 0;
}
//==========================================================
//==========================================================
int stat_init(struct statistics *stat) {
  stat->n=0;
  stat->sum=0;
  stat->sum2=0;
  return 0;
}
//==========================================================
//==========================================================
int statistics(struct statistics *stat,double x) {
  stat->n=stat->n++;
  stat->sum=stat->sum+x;
  stat->sum2=stat->sum2+x*x;
  stat->mean=stat->sum/(stat->n);
  stat->stdev=sqrt(stat->sum2/(stat->n)-stat->sum*stat->sum/((stat->n)*(stat->n)));  
  return 0;
}
//==========================================================
//==========================================================
double gauss(double x, double mean, double sdev) {
  double rtn=exp(-(mean-x)*(mean-x)/(2.0*sdev*sdev));
  //sqrt(2PI)=2.506628274631
  rtn=rtn/(2.506628274631*sdev);
  return rtn;
}
//==========================================================
//==========================================================
int readpilot(struct parameters *p,struct run *r,struct source *src,struct array *ta) {
  char buff[MAXMSGLENGTH];
  FILE *farray;
  farray=fopen(r->pilotfile,"r");
  int ieof=0;
  char flag[FLAGLENGTH];
  while(!ieof) {
    look_for_flag(&ieof,flag,farray);
    if(ieof) break;
    
    if(strstr(flag,"TSTEP")!=NULL) {
      fscanf(farray,"%lf",&(p->dt));
      p->dt=p->dt*1e-9;
    }

    if(strstr(flag,"INTTM")!=NULL) {
      fscanf(farray,"%lf",&(p->T));
      p->T=p->T*1e-9;
    }

    if(strstr(flag,"PDURM")!=NULL) {
      fscanf(farray,"%lf",&(p->pdur));
    }

    if(strstr(flag,"PDURT")!=NULL) {
      fscanf(farray,"%lf",&(p->pulsedur));
    }

    if(strstr(flag,"PAMPL")!=NULL) {
      fscanf(farray,"%lf",&(p->pampl));
    }

    if(strstr(flag,"BGLIT")!=NULL) {
      fscanf(farray,"%lf",&(p->bglight));
    }

    if(strstr(flag,"EXNOS")!=NULL) {
      fscanf(farray,"%lf",&(p->exnoise));
    }

    if(strstr(flag,"EXNNF")!=NULL) {
      fscanf(farray,"%d",&(r->exnonoff));
    }

    if(strstr(flag,"BSLNF")!=NULL) {
      fscanf(farray,"%d",&(r->bslonoff));
    }

    if(strstr(flag,"MAXBL")!=NULL) {
      fscanf(farray,"%lf",&(p->maxbl));
    }

    if(strstr(flag,"BLSTP")!=NULL) {
      fscanf(farray,"%lf",&(p->blstep));
    }

    if(strstr(flag,"CTIME")!=NULL) {
      fscanf(farray,"%lf",&(p->cohertime));
      p->cohertime=p->cohertime*1e-9;
    }

    if(strstr(flag,"SFLUX")!=NULL) {
      fscanf(farray,"%lf",&(src->flux));
      src->flux=src->flux*1e9;
    }

    if(strstr(flag,"FREQN")!=NULL) {
      fscanf(farray,"%lf",&(src->meanfreq));
      src->meanfreq=src->meanfreq*1e14;
    }

    if(strstr(flag,"SRADI")!=NULL) {
      fscanf(farray,"%lf",&(src->rad));
      src->rad=src->rad*1e-9;
    }

    if(strstr(flag,"TERAD")!=NULL) {
      fscanf(farray,"%lf",&(ta->tr));
    }

    if(strstr(flag,"STEPS")!=NULL) {
      fscanf(farray,"%lf",&(ta->stepsz));
    }


  }//End of while loop
  fclose(farray);
  return 0;
}
//==========================================================
//==========================================================
