#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define KB 1

typedef enum{berendsen,andersen,langevin,bussi} t_type;
typedef enum{FALSE,TRUE} boolean;
void read_freq_file(char filename[],double freq_vett[],int *n_particles,FILE *file_pointer);
void initialize_v(double v[],double m,int num_particles,double T);
void initialize_q(double q[],int num_particles,double m,double T,double freq[]);
void update_qp(double q[],double v[],double ***m,int num_particles,double step,double freq[]);
void thermostat(double velocities[],int num_particles,double m,double targetT,double delta_t,int dof,double tau,t_type type,boolean update_p);
void sort_increasing(double v[],int length);
double oscill_energy(double q,double v,double w,double mass);
double max(double v[],int length);
void print_traj(double q[],double v[],int particle_index,FILE* filept);
void comp_rot_matr(double ***m,double delta_t,double w[],int num_particles);

double ran1();
double gasdev();
double gamdev(const int ia);
double resamplekin_sumnoises(int nn);
double resamplekin(double kk,double sigma, int ndeg, double taut);
void print_traj2(double v1[],double v2[],double v3[],int l,int m,int n,FILE* filept);
int main(int argc, char *argv[])
{

FILE *f_freq,*energy_before,*energy_after,*trajectory1,*trajectory2,*pcorr,*qcorr,*ph,*traj3,*mom3f,*mom4f;
double *pos;     //phonon positions vector
double *vel;     //phonon velocities vector     //to change with dynamic allocation
double *freq_vett,t_tot,dt,tau,tau1,tau2,*energy1,*energy2,***matr_array,**pcorr_matrix,**qcorr_matrix,***mom3,****mom4;//*avg,*z;         //phonon frequency vector, energy without thermostat, energy with thermostat
int n_particles=0,n_steps,i,j,k,l,n;
double targetT=20.0,Tstart=10.0,m=1.0;  //target energy
char fname[15]="phase",prtc_numb[6];
t_type thermo_type;
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
matr_array=(double ***)malloc(sizeof(double **) * n_particles);
mom3=(double ***)malloc(sizeof(double ** ) * n_particles);
mom4=(double ****)malloc(sizeof(double***)*n_particles);
for(i=0;i<n_particles;i++)
{
    matr_array[i]=malloc(2*sizeof(double *));
    mom3[i]=malloc(n_particles*sizeof(double *));
    mom4[i]=malloc(n_particles*sizeof(double **));
    for(j=0;j<2;j++)
        matr_array[i][j]=malloc(2*sizeof(double));
    for(j=0;j<n_particles;j++)
        {
        mom3[i][j]=malloc(n_particles*sizeof(double));
        mom4[i][j]=malloc(n_particles*sizeof(double *));
        for(k=0;k<n_particles;k++)
            mom4[i][j][k]=malloc(n_particles*sizeof(double));
        }

}

pcorr_matrix=(double **)malloc(sizeof(double *) * n_particles);
qcorr_matrix=(double **)malloc(sizeof(double *) * n_particles);
for(i=0;i<n_particles;i++)
        {
        pcorr_matrix[i]=malloc(n_particles*sizeof(double));
        qcorr_matrix[i]=malloc(n_particles*sizeof(double));
        }

energy1= (double *)malloc(sizeof(double) * n_particles);
energy2= (double *)malloc(sizeof(double) * n_particles);

initialize_q(pos,n_particles,m,Tstart,freq_vett);
initialize_v(vel,m,n_particles,Tstart);
sort_increasing(freq_vett,n_particles);

t_tot=200000;
tau1=2*M_PI/max(freq_vett,n_particles);
tau2=10; /*is the relaxation time of the thermostat to the target temperature*/
tau=tau2;
if(tau2>tau1)
	tau=tau1;
dt=tau/10;
energy_after=fopen("ber_t_03dt.dat","w");
n_steps=t_tot/dt+1;

for(j=0;j<n_particles;j++)
{
    energy1[j]=0.0;
    energy2[j]=0.0;
    for(k=0;k<n_particles;k++)
        {
        pcorr_matrix[j][k]=0.0;
        qcorr_matrix[j][k]=0.0;
        for(l=0;l<n_particles;l++)
            {
            mom3[j][k][l]=0.0;
            for(n=0;n<n_particles;n++)
                mom4[j][k][l][n]=0.0;
            }
        }
}
comp_rot_matr(matr_array,dt,freq_vett,n_particles);

thermo_type=berendsen;

trajectory2= fopen("trajectory2.dat","w");
traj3=fopen("traj3d_ber.dat","w");
int counter=0;
boolean UPDATE_P=FALSE;
for(i=0;i<n_steps;i++)
{
    print_traj(pos,vel,0,trajectory2);
    print_traj2(vel,vel,vel,1,2,3,traj3);
    if(counter>=(tau/dt))
    {
    UPDATE_P=TRUE;
    counter=0;
    }
    thermostat(vel,n_particles,m,targetT,dt/2,2,tau2,thermo_type,UPDATE_P);
    update_qp(pos,vel,matr_array,n_particles,dt,freq_vett);
    thermostat(vel,n_particles,m,targetT,dt/2,2,tau2,thermo_type,UPDATE_P);
    counter++;
    UPDATE_P=FALSE;
    for(j=0;j<n_particles;j++)
    {
        energy2[j]+=oscill_energy(pos[j],vel[j],freq_vett[j],m)/n_steps;
        for(k=0;k<n_particles;k++)
            {
                    pcorr_matrix[j][k]+=pow(m,2)*vel[j]*vel[k]/n_steps;
                    qcorr_matrix[j][k]+=vel[j]*vel[k]/n_steps;
                    for(l=0;l<n_particles;l++)
                                {
                                   mom3[j][k][l]+=vel[j]*vel[k]*vel[l]/n_steps;
                                    for(n=0;n<n_particles;n++)
                                        mom4[j][k][l][n]+=vel[j]*vel[k]*vel[l]*vel[n]/n_steps;
                                }

                    /*sprintf(prtc_numb,"%d%c%d",j+1,'_',k+1);
                    strcat(fname,prtc_numb);
                    ph=fopen(fname,"a");
                    fprintf(ph,"%d %f\n",i,atan(vel[j]/pos[j])-atan(vel[k]/pos[k]));
                    fclose(ph);
                    strcpy(fname,"phase");*/
            }
    }
}

fclose(trajectory2);
energy_before=fopen("energy_before.dat","w");

pcorr=fopen("pcorr_matr.dat","w");
qcorr=fopen("qcorr_matr.dat","w");
mom3f=fopen("mom3.dat","w");
mom4f=fopen("mom4.dat","w");
if( energy_before!=NULL && energy_after!=NULL && qcorr!=NULL && pcorr!=NULL && mom3f !=NULL )
{
    for(i=0;i<n_particles;i++)
    {
        fprintf(energy_before,"%f %f\n",freq_vett[i],energy1[i]);
        fprintf(energy_after,"%f %f\n",freq_vett[i],energy2[i]);
        for(j=0;j<n_particles;j++)
            {
            fprintf(pcorr,"%d %d %f\n",i,j,pcorr_matrix[i][j]);
            fprintf(qcorr,"%d %d %f\n",i,j,qcorr_matrix[i][j]);
            for(l=0;l<n_particles;l++)
                {
                fprintf(mom3f,"%d%d%d: %f\n",i,j,l,mom3[i][j][l]);
                for(n=0;n<n_particles;n++)
                    fprintf(mom4f,"%d%d%d%d: %f\n",i,j,l,n,mom4[i][j][l][n]);
                }

            }
    }
}
/*
//Test to plot trajectory without thermostat (only Hamiltonian steps)

initialize_q(pos,n_particles,m,Tstart,freq_vett);
initialize_v(vel,m,n_particles,Tstart);
trajectory1= fopen("trajectory1.dat","w");
for(i=0;i<n_steps;i++)
    {

        print_traj(pos,vel,0,trajectory1);
        for(j=0;j<n_particles;j++)
                energy1[j]+=oscill_energy(pos[j],vel[j],freq_vett[j],m)/n_steps;
        update_qp(pos,vel,matr_array,n_particles,dt,freq_vett);

    }
fclose(trajectory1);
*/



fclose(energy_before);
fclose(energy_after);
fclose(pcorr);
fclose(qcorr);


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
void thermostat(double velocities[],int num_particles,double m,double targetT,double delta_t,int dof,double tau,t_type type,boolean update_p)
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
            velocities[i]=velocities[i]*c1+c2*gasdev()*radq; //because sqrt(T)=v[i]
        }
        else
        {
            if(type==andersen && update_p==TRUE )
            {
            double radq=sqrt(KB*targetT/m);
            for(i=0;i<num_particles;i++)
                velocities[i]=gasdev()*radq;
            }
            else
                {
                    if(type==bussi)
                        {
                        double kt=0.0,lambda,kt2,target=num_particles*KB*0.5*targetT;  //?
                        for(i=0;i<num_particles;i++)
                            kt+=0.5*m*velocities[i]*velocities[i];
                        kt2=resamplekin(kt,target,num_particles,tau/delta_t);
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

void initialize_v(double v[],double m,int num_particles,double T)
{
    int i;
    double radq=sqrt(KB*T/m);
    for(i=0;i<num_particles;i++)
        v[i]=gasdev()*radq;
}
void initialize_q(double q[],int num_particles,double m,double T,double freq[])  //TODO: include SEED
{
    int i;
    double k;
    for(i=0;i<num_particles;i++)
    {
        k=m*pow(freq[i],2);
        q[i]=gasdev()*sqrt(KB*T/k);
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
double resamplekin(double kk,double sigma, int ndeg, double taut){
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
  rr = gasdev();
  return kk + (1.0-factor)* (sigma*(resamplekin_sumnoises(ndeg-1)+rr*rr)/ndeg-kk)
            + 2.0*rr*sqrt(kk*sigma/ndeg*(1.0-factor)*factor);
}

double resamplekin_sumnoises(int nn){
/*
  returns the sum of n independent gaussian noises squared
   (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
*/
  double rr;
  if(nn==0) {
    return 0.0;
  } else if(nn==1) {
    rr=gasdev();
    return rr*rr;
  } else if(nn%2==0) {
    return 2.0*gamdev(nn/2);
  } else {
    rr=gasdev();
    return 2.0*gamdev((nn-1)/2) + rr*rr;
  }
}



double gamdev(const int ia)
{
	int j;
	double am,e,s,v1,v2,x,y;

	if (ia < 1) {}; // FATAL ERROR
	if (ia < 6) {
		x=1.0;
		for (j=1;j<=ia;j++) x *= ran1();
		x = -log(x);
	} else {
		do {
			do {
				do {
					v1=ran1();
					v2=2.0*ran1()-1.0;
				} while (v1*v1+v2*v2 > 1.0);
				y=v2/v1;
				am=ia-1;
				s=sqrt(2.0*am+1.0);
				x=s*y+am;
			} while (x <= 0.0);
			e=(1.0+y*y)*exp(am*log(x/am)-s*y);
		} while (ran1() > e);
	}
	return x;
}



double gasdev()
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (iset == 0) {
		do {
			v1=2.0*ran1()-1.0;
			v2=2.0*ran1()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

#define NTAB 32
double ran1()
{
	const int IA=16807,IM=2147483647,IQ=127773,IR=2836;
	const int NDIV=(1+(IM-1)/NTAB);
	const double EPS=3.0e-16,AM=1.0/IM,RNMX=(1.0-EPS);
	static int iy=0;
	//int iv[NTAB];
	static int iv[NTAB];
	int j,k;
	double temp;
        static int idum=0; /* ATTENTION: THE SEED IS HARDCODED */

	if (idum <= 0 || !iy) {
		if (-idum < 1) idum=1;
		else idum = -idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ;
			idum=IA*(idum-k*IQ)-IR*k;
			if (idum < 0) idum += IM;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
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
void print_traj2(double v1[],double v2[],double v3[],int l,int m,int n,FILE* filept)
{

if (filept != NULL)
        fprintf(filept,"%f %f %f\n",v1[l],v2[m],v3[n]);
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
