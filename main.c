#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include<time.h>
#define KB 1
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

typedef enum{berendsen,andersen,langevin,bussi} t_type;
typedef enum{FALSE,TRUE} boolean;
void read_freq_file(char filename[],double freq_vett[],int *n_particles,FILE *file_pointer);
void initialize_v(double v[],double m,int num_particles,double T,long *seed);
void initialize_q(double q[],int num_particles,double m,double T,double freq[],long *seed);
void update_qp(double q[],double v[],double ***m,int num_particles,double step,double freq[]);
void thermostat(double velocities[],int num_particles,double m,double targetT,double delta_t,int dof,double tau,t_type type,boolean update_p,long *seed);
void sort_increasing(double v[],int length);
double oscill_energy(double q,double v,double w,double mass);
double max(double v[],int length);
void print_traj(double q[],double v[],int particle_index,FILE* filept);
void comp_rot_matr(double ***m,double delta_t,double w[],int num_particles);
double resamplekin_sumnoises(int nn,long *seed);
double resamplekin(double kk,double sigma, int ndeg, double taut,long *seed);
float ran1(long *idum);
float gasdev(long *idum);
float gamdev(int ia, long *idum);

int main(int argc, char *argv[])
{

FILE *f_freq,*energy_before,*energy_after,*trajectory1,*trajectory2,*pcorr,*qcorr,*ph,*energy,*f_seed;
double *pos;     //phonon positions vector
double *vel;     //phonon velocities vector     //to change with dynamic allocation
double *freq_vett,t_tot,dt,tau,*energy1,*energy2,***matr_array,**pcorr_matrix,**qcorr_matrix,*sumsquare,*mean,*stdev;         //phonon frequency vector, energy without thermostat, energy with thermostat
int n_particles=0,n_steps,i,j,k,nruns=1,nr,intseed,narg,PLOT_TRAJ=0,np,COMP_EN=0;
double targetT=20.0,Tstart=10.0,m=1.0,a,tau1,tau2;
char fname[15]="phase",prtc_numb[6],seeds_fn[20],thermo[15],*ptr,data[17],trjname[21]="trj",ename[21]="en";
t_type thermo_type;
long SEED,*seeds;
time_t rawtime;
struct tm* ltime;
t_tot=atoi(argv[2]);
if(strcmp(argv[3],"bussi")==0){
    thermo_type=3;
    strcpy(thermo,"bussi");}
if(strcmp(argv[3],"langevin")==0){
    thermo_type=2;
    strcpy(thermo,"langevin");}
if(strcmp(argv[3],"andersen")==0){
    thermo_type=1;
    strcpy(thermo,"andersen");}
if(strcmp(argv[3],"berendsen")==0){
    thermo_type=0;
    strcpy(thermo,"berendsen");}
if(argc>4)
    tau2=strtod(argv[4],&ptr);
if(argc>5)
    PLOT_TRAJ=atoi(argv[5]);
if(argc>6)
    COMP_EN=atoi(argv[6]);
if(argc>7)
{
    nruns=atoi(argv[7]);
    strcpy(seeds_fn,argv[8]);
}
time (&rawtime);
ltime = localtime(&rawtime);
sprintf(data,"%d%d%d%d%d%s",ltime->tm_mday,ltime->tm_mon+1,ltime->tm_year+1900,ltime->tm_hour,ltime->tm_min,".dat");

//Frequency file reading
   if ((f_freq= fopen(argv[1],"r")) != NULL)
   {
       fscanf(f_freq,"%d",&n_particles );
       freq_vett=(double *)malloc(sizeof(double) * n_particles   );
       read_freq_file(argv[1],freq_vett,&n_particles,f_freq);
   }
   else
    {
        printf("\nError opening %s",argv[1]);
        exit(1);
    }
//Memory allocation for positions and velocities

pos= (double *)malloc(sizeof(double) * n_particles);
vel= (double *)malloc(sizeof(double) * n_particles);
energy1= (double *)malloc(sizeof(double) * n_particles);
energy2= (double *)malloc(sizeof(double) * n_particles);
sumsquare=(double *)malloc(sizeof(double) * n_particles);
stdev=(double *)malloc(sizeof(double) * n_particles);
mean=(double *)malloc(sizeof(double) * n_particles);
matr_array=(double ***)malloc(sizeof(double **) * n_particles);
for(i=0;i<n_particles;i++)
    matr_array[i]=malloc(2*sizeof(double *));
for(i=0;i<n_particles;i++)
    for(j=0;j<2;j++)
        matr_array[i][j]=malloc(2*sizeof(double));
/*pcorr_matrix=(double **)malloc(sizeof(double *) * n_particles);
qcorr_matrix=(double **)malloc(sizeof(double *) * n_particles);
for(i=0;i<n_particles;i++)
        {
        pcorr_matrix[i]=malloc(n_particles*sizeof(double));
        qcorr_matrix[i]=malloc(n_particles*sizeof(double));
        }*/

sort_increasing(freq_vett,n_particles);
//t_tot=100000;
tau1=2*M_PI/max(freq_vett,n_particles);
tau=tau2;
if(tau2>tau1)
	tau=tau1;
dt=tau/10;
n_steps=t_tot/dt+1;
comp_rot_matr(matr_array,dt,freq_vett,n_particles);
if(nruns>1)
{
    if((f_seed=fopen(seeds_fn,"r"))!=NULL)
        {
        fscanf(f_seed,"%d",&nruns);
        seeds=(long *)malloc(sizeof(long) * nruns);
        for(nr=0;nr<nruns;nr++)
            {
            fscanf(f_seed,"%d",&intseed);
            seeds[nr]=(long)intseed;
            }

        }
    for(i=0;i<n_particles;i++)
    {
        sumsquare[j]=0.0;
        mean[j]=0.0;
    }
fclose(f_seed);
}

nr=0;
if(PLOT_TRAJ==1)
{
    
    strcat(trjname,data);
    trajectory2= fopen(trjname,"w");
    fprintf(trajectory2,"#%d %d\n#Thermostat used: %s\n#n_steps=%d,tau=%f,nruns=%d\n#Freq.file used:\n",n_steps,2*n_particles,thermo,n_steps,tau2,nruns);
    for(np=0;np<n_particles;np++)
        fprintf(trajectory2,"#%f\n",freq_vett[np]);
}
while(nr<nruns)
{
	SEED=(long)-1;
	a=ran1(&SEED);
	if(nruns>1)
        SEED=seeds[nr];
    else
        SEED=(long)0;
	initialize_q(pos,n_particles,m,Tstart,freq_vett,&SEED);
    initialize_v(vel,m,n_particles,Tstart,&SEED);
    if(COMP_EN==1)
    {
    for(j=0;j<n_particles;j++)
        {
            //energy1[j]=0.0;
            energy2[j]=0.0;
           /* for(k=0;k<n_particles;k++)
                {
                pcorr_matrix[j][k]=0.0;
                qcorr_matrix[j][k]=0.0;
                }*/
        }
    }

        //thermo_type=bussi;


        int counter=0;
        boolean UPDATE_P=FALSE;
        for(i=0;i<n_steps;i++)
        {
            if(PLOT_TRAJ==1 && nruns==1)
            {
                for(np=0;np<n_particles;np++)
                fprintf(trajectory2,"%f %f ",pos[np],vel[np]);
                fprintf(trajectory2,"\n");
            }
            if(counter>=(tau2/dt))
            {
            UPDATE_P=TRUE;
            counter=0;
            }
            thermostat(vel,n_particles,m,targetT,dt/2,2,tau2,thermo_type,UPDATE_P,&SEED);
            update_qp(pos,vel,matr_array,n_particles,dt,freq_vett);
            thermostat(vel,n_particles,m,targetT,dt/2,2,tau2,thermo_type,UPDATE_P,&SEED);
            counter++;
            UPDATE_P=FALSE;
	    if(COMP_EN==1)
	    {
	      for(j=0;j<n_particles;j++)
	      {
		  energy2[j]+=oscill_energy(pos[j],vel[j],freq_vett[j],m)/n_steps;
		/* for(k=0;k<n_particles;k++)
		      {
			      pcorr_matrix[j][k]+=pow(m,2)*vel[j]*vel[k]/n_steps;
			      qcorr_matrix[j][k]+=vel[j]*vel[k]/n_steps;*/
			    /* sprintf(prtc_numb,"%d%c%d",j+1,'_',k+1);
			      strcat(fname,prtc_numb);
			      ph=fopen(fname,"a");
			      fprintf(ph,"%d %f\n",i,atan(vel[j]/pos[j])-atan(vel[k]/pos[k]));
			      fclose(ph);
			      strcpy(fname,"phase");*/
		      //}
	      }
	    }
        }


        //energy_before=fopen("energy_before.dat","w");
        //energy_after=fopen("energy_after.dat","w");
        //pcorr=fopen("pcorr_matr.dat","w");
        //qcorr=fopen("qcorr_matr.dat","w");
        //if( energy_before!=NULL && energy_after!=NULL && qcorr!=NULL && pcorr!=NULL)
            if(nruns>1)
            {
                for(i=0;i<n_particles;i++)
                {
                    //fprintf(energy_before,"%f %f\n",freq_vett[i],energy1[i]);
                    //fprintf(energy_after,"%f %f\n",freq_vett[i],energy2[i]);
                    sumsquare[i]+=pow(energy2[i],2)/nruns;
                    mean[i]+=energy2[i]/nruns;
                    /*for(j=0;j<n_particles;j++)
                        {
                        fprintf(pcorr,"%d %d %f\n",i,j,pcorr_matrix[i][j]);
                        fprintf(qcorr,"%d %d %f\n",i,j,qcorr_matrix[i][j]);
                        }*/
                }
            }

        //Test to plot trajectory without thermostat (only Hamiltonian steps)
	/*
        initialize_q(pos,n_particles,m,Tstart,freq_vett,&SEED);
        initialize_v(vel,m,n_particles,Tstart,&SEED);
        trajectory1= fopen("trajectory1.dat","w");
        for(i=0;i<n_steps;i++)
            {

                print_traj(pos,vel,0,trajectory1);
                for(j=0;j<n_particles;j++)
                        energy1[j]+=oscill_energy(pos[j],vel[j],freq_vett[j],m)/n_steps;
                update_qp(pos,vel,matr_array,n_particles,dt,freq_vett);

            }
        fclose(trajectory1);*/

        //fclose(energy_before);
       // fclose(energy_after);
        //fclose(pcorr);
        //fclose(qcorr);

nr++;
}
if(PLOT_TRAJ==1)
    fclose(trajectory2);
if(COMP_EN==1)
{
strcat(ename,data);
energy=fopen(ename,"w");
fprintf(energy,"#Thermostat used: %s\n#n_steps=%d,tau=%f,nruns=%d\n#Freq.file used:\n",thermo,n_steps,tau2,nruns);
for(np=0;np<n_particles;np++)
    fprintf(energy,"#%f\n",freq_vett[np]);
}
for(i=0;i<n_particles;i++)
{
    if(nruns>1)
    {
      stdev[i]=(sqrt(sumsquare[i]-pow(mean[i],2)))/sqrt(nruns);
      fprintf(energy,"%f %f %f\n",freq_vett[i],mean[i],stdev[i]);
    }
    else
    { 
      if(COMP_EN==1)
	fprintf(energy,"%f %f\n",freq_vett[i],energy2[i]);
      
    }
}
if(COMP_EN==1)
  fclose(energy);
return 0;

}









void read_freq_file(char filename[],double freq_vett[],int *n_particles,FILE *file_pointer)
{


   float freq;
   int np,part_check=0,i,j;

    //lettura file frequenze
          i=0;
       while(i<*n_particles)
       {
           if((fscanf(file_pointer,"%f%d", &freq,&np))==2)  //np = # phonons with frequency freq
            {
               part_check+=np;
               for(j=0;j<np;j++)
               {
                  freq_vett[i]=freq;
                  i++;
               }
            }
            else
            {
                freq_vett[i]=freq;
                i++;
            }
       }
fclose(file_pointer);
}

void update_qp(double q[],double v[],double ***m,int num_particles,double step,double freq[])
{
//q: positions vector
//p: velocities vector
//freq: Frequencies vector
//step: step of Hamiltonian dynamic
//m:array of rotation matrices
    int i;
    float qt,vt;
    for(i=0;i<num_particles;i++)
    {
        qt=q[i];
        vt=v[i];
        q[i]=qt*m[i][0][0]+vt*m[i][0][1];
        v[i]=qt*m[i][1][0]+vt*m[i][1][1];

    }
}
void thermostat(double velocities[],int num_particles,double m,double targetT,double delta_t,int dof,double tau,t_type type,boolean update_p,long *seed)
{
    int i;

    if (type==berendsen)
    {
    double kt=0.0,lambda,kt2,target=0.5*KB*num_particles*targetT;
    for(i=0;i<num_particles;i++)
            kt+=0.5*m*velocities[i]*velocities[i];
    kt2=target+(kt-target)*exp(-delta_t/tau);
    lambda=sqrt(kt2/kt);
    for(i=0;i<num_particles;i++)
            velocities[i]=velocities[i]*lambda;
    }
    else
    {
        if(type==langevin)
        {
        double c1,c2,radq=sqrt(KB*targetT/m);
        c1=exp(delta_t/tau);
        c2=sqrt(1-c1*c1);
        for(i=0;i<num_particles;i++)
            velocities[i]=velocities[i]*c1+c2*gasdev(seed)*radq;
        }
        else
        {
            if(type==andersen && update_p==TRUE )
            {
            double radq=sqrt(KB*targetT/m);
            for(i=0;i<num_particles;i++)
                velocities[i]=gasdev(seed)*radq;
            }
            else
                {
                    if(type==bussi)
                        {
                        double kt=0.0,lambda,kt2,target=num_particles*KB*0.5*targetT;
                        for(i=0;i<num_particles;i++)
                            kt+=0.5*m*velocities[i]*velocities[i];
                        kt2=resamplekin(kt,target,num_particles,tau/delta_t,seed);
                        lambda=sqrt(kt2/kt);
                        for(i=0;i<num_particles;i++)
                            velocities[i]=velocities[i]*lambda;

                        }
                    else
                        return;
                }

        }

    }

}

void initialize_v(double v[],double m,int num_particles,double T,long *seed)
{
    int i;
    double radq=sqrt(KB*T/m);
    for(i=0;i<num_particles;i++)
        v[i]=gasdev(seed)*radq;
}
void initialize_q(double q[],int num_particles,double m,double T,double freq[],long *seed)  //TODO: include SEED
{
    int i;
    double k;
    for(i=0;i<num_particles;i++)
    {
        k=m*pow(freq[i],2);
        q[i]=gasdev(seed)*sqrt(KB*T/k);
    }
}
void sort_increasing(double v[],int length)
{
    int i,j;
    float temp;
    for(i=0;i<length-1;i++)
        for(j=i+1;j<length;j++)
        {
            if(v[j]<=v[i])
            {
               temp=v[i];
               v[i]=v[j];
               v[j]=temp;
            }
        }
}
double oscill_energy(double q,double v,double w,double mass)
{
    double e,k;
    k=mass*pow(w,2);
    e=0.5*(mass*pow(v,2)+k*pow(q,2)) ;
    return e;
}
double resamplekin(double kk,double sigma, int ndeg, double taut,long *seed){
/*
  kk:    present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
  sigma: target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
  ndeg:  number of degrees of freedom of the atoms to be thermalized
  taut:  relaxation time of the thermostat, in units of 'how often this routine is called'
*/
  double factor,rr;
  if(taut>0.1){
    factor=exp(-1.0/taut);
  } else{
    factor=0.0;
  }
  rr = gasdev(seed);
  return kk + (1.0-factor)* (sigma*(resamplekin_sumnoises(ndeg-1,seed)+rr*rr)/ndeg-kk)
            + 2.0*rr*sqrt(kk*sigma/ndeg*(1.0-factor)*factor);
}

double resamplekin_sumnoises(int nn,long *seed){
/*
  returns the sum of n independent gaussian noises squared
   (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
*/
  double rr;
  if(nn==0) {
    return 0.0;
  } else if(nn==1) {
    rr=gasdev(seed);
    return rr*rr;
  } else if(nn%2==0) {
    return 2.0*gamdev(nn/2,seed);
  } else {
    rr=gasdev(seed);
    return 2.0*gamdev((nn-1)/2,seed) + rr*rr;
  }
}





double max(double v[],int length)
{
    double m=v[0];
    int i;
    for(i=0;i<length;i++)
    {
        if(v[i]>=m)
            m=v[i];
    }
    return m;
}
void print_traj(double q[],double v[],int particle_index,FILE* filept)
{

if (filept != NULL)
        fprintf(filept,"%f %f\n",q[particle_index],v[particle_index]);
else
    printf("Error creating trajectory file\n");

}
void comp_rot_matr(double ***m,double delta_t,double w[],int num_particles)
{
    int i;
    double s,c;
    for(i=0;i<num_particles;i++)
    {
        c=cos(w[i]*delta_t);
        s=sin(w[i]*delta_t);
        m[i][0][0]=c;
        m[i][0][1]=s/w[i];
        m[i][1][0]=-w[i]*s;
        m[i][1][1]=c;
    }
}


float ran1(long *idum)
/*Minimal” random number generator of Park and Miller with Bays-Durham shuffle and added
safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
successive deviates in a sequence. RNMX should approximate the largest floating value that is
less than 1.*/
{
int j;
long k;
static long iy=0;
static long iv[NTAB];
float temp;

if (*idum <= 0 || !iy)
{
	if (-(*idum) < 1)
		*idum=1;
	else
		*idum = -(*idum);
	for (j=NTAB+7;j>=0;j--)
	{
		k=(*idum)/IQ;
		*idum=IA*(*idum-k*IQ)-IR*k;
		if (*idum < 0)
			*idum += IM;
		if (j < NTAB)
			iv[j] = *idum;
	}
	iy=iv[0];
}
k=(*idum)/IQ;
*idum=IA*(*idum-k*IQ)-IR*k;
if (*idum < 0)
	*idum += IM;
j=iy/NDIV;
iy=iv[j];
iv[j] = *idum;
if ((temp=AM*iy) > RNMX)
	return RNMX;
else
	return temp;
}
float gasdev(long *idum)
/*Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum) as the source of uniform deviates.*/
{
float ran1(long *idum);
static int iset=0;
static float gset;
float fac,rsq,v1,v2;
if (*idum < 0)
	iset=0;
if (iset == 0)
	{
	do
		{
		v1=2.0*ran1(idum)-1.0;
		v2=2.0*ran1(idum)-1.0;
		rsq=v1*v1+v2*v2;
		}
	while (rsq >= 1.0 || rsq == 0.0);
	fac=sqrt(-2.0*log(rsq)/rsq);
	gset=v1*fac;
	iset=1;
	return v2*fac;
	}
else
	{
	iset=0;
	return gset;
	}
}

float gamdev(int ia, long *idum)
/*Returns a deviate distributed as a gamma distribution of integer order ia, i.e., a waiting time
+to the iath event in a Poisson process of unit mean, using ran1(idum) as the source of
+uniform deviates.*/
{
float ran1(long *idum);
//void nrerror(char error_text[]);
int j;
float am,e,s,v1,v2,x,y;
if (ia < 1) {} //ERROR nrerror("Error in routine gamdev");
if (ia < 6) {
x=1.0;
for (j=1;j<=ia;j++) x *= ran1(idum);
x = -log(x);
} else {
do {
do {
do {
v1=ran1(idum);
v2=2.0*ran1(idum)-1.0;
} while (v1*v1+v2*v2 > 1.0);
y=v2/v1;
am=ia-1;
s=sqrt(2.0*am+1.0);
x=s*y+am;
} while (x <= 0.0);
e=(1.0+y*y)*exp(am*log(x/am)-s*y);
} while (ran1(idum) > e);
}
return x;

}