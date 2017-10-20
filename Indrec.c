// Individually based forward in time program of evolution of recombination rate between, and of, plasticity and target locus under the genomic storage effect.  

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.141592654
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS (1.2E-07)
#define RNMX (1.0-EPS)
#define E 2.71828182845904523536
#define NLs 3                                         
#define Np0 25000                           ///////////////////////////////////    define population size here    ////////////////////////////////////////////////////////

long double ran1(long *idum); long double poidev(long double xm, long *idum); long double gammln(long double xx);

int main(void){

long REPS, TIME, BURNIN, seed;     int PHASE, sscd;    double REC, REC2;           
FILE * par, * inp , * out2;
par = fopen("pars1", "r"); inp = fopen("trajectory", "r"); out2 = fopen("outfr", "w"); fscanf(par, "%ld", &seed);

//////////////////////////////////// define parameters : ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double mu = 0.000004;         ///////////////////////////////////////define  mutation rate                                                                            
REPS =1;                      /// replicates, recomended to use 1 but many parallel simulations
TIME = 200*Np0;               /// time of simulation, the statinary distribution is recorded over the 100N last generations
BURNIN = 100*Np0;             /// burning period
seed = -1000 - seed;          /// to generate negative seed 
PHASE = 10;                   ////////////////////////////////////////define PHASE
REC = 0.5;                    //////////////////////////////////////// define recombination rate between the rec modifier and pla modifier locus
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////all of the parameters are now defined///////////////

long double mutants, tra[3*PHASE], hetrec = 0.0, Het = 0.0, ps, Phet = 0.0, Plas, muts = mu*Np0;
static double parents[NLs+2][Np0], offspring[NLs][Np0], sumrec[51];
double fits, maxf, cond;
long ispoly, k, c, i, j, y, recvec[51], mutated[NLs], nmutated, pick1, pick2, plafreq, freq, notfinished=0, nopoly =0, isthere, pladone, fixed =0, starttime = 0, losstime = 0, fixtime =0, numfix =0, numloss =0, N2 = (Np0/2);
int z, m, ti, maxti, notgood=0, second;

for (k = 0; k < (3*PHASE); k++){ fscanf(inp, "%Lf", &tra[k]);}   // read in trajectory of choice

for(c =0 ; c < REPS; c++){     //////////////////////////////////////////////////////////////////////////////////////////// REPLICATES START

m = 0;

REC2 = ran1(&seed)*0.5;
	for(j=0; j<Np0; j++){
                        parents[0][j]=REC2;
 }

	for(i=1; i<NLs; i++){
	for(j=0; j<Np0; j++){
			parents[i][j]=-1.0;
	}}

for(k=0; k<TIME; k++){   						//////////////////// TIME STARTS

freq = 0; plafreq = 0; ispoly = 0;

//////////////////////recurrent mutation:            

mutants = poidev(muts ,&seed);
for(y = 0; y < mutants; y++){
j = (long)(ran1(&seed)*Np0);
parents[0][j]= ran1(&seed)*0.5;        // mutation at a recombination (Pla with Target) locus
}

for(i=1; i<NLs; i++){       /// rest of the loci
mutants = poidev(muts ,&seed);
for(y = 0; y < mutants; y++){
j = (long)(ran1(&seed)*Np0);
parents[i][j]=-parents[i][j];
}
}                                          
 

/////////////////selection

if(m == PHASE) m = 0;
fits = tra[m];

for(j=0; j<Np0; j++){
parents[NLs][j]=0.0;
}

for(j=0; j<Np0; j++){
if(parents[1][j]!=1.0){ 
parents[NLs][j]+=log(1.0+parents[2][j]*fits);
}
}

maxf = parents[NLs][0];
for(j=1; j<Np0; j++){
if(parents[NLs][j] > maxf)maxf = parents[NLs][j];
}


for(j=0; j<Np0; j++){
parents[NLs+1][j]=pow(E,parents[NLs][j])/pow(E,maxf);
}

m += 1;

/////////////////////////////// recombination and repro

for(j=0; j<N2; j++){

do{
pick1 = (long)(ran1(&seed)*Np0);
}while(ran1(&seed)> parents[NLs+1][pick1]);

do{
pick2 = (long)(ran1(&seed)*Np0);
}while(ran1(&seed)> parents[NLs+1][pick2]);

offspring[0][j]= parents[0][pick1];   
offspring[0][j+N2]=parents[0][pick2];

if(ran1(&seed)<REC) {isthere = pick1; pick1 = pick2; pick2 = isthere;}
offspring[1][j]= parents[1][pick1];
offspring[1][j+N2] = parents[1][pick2];


if(ran1(&seed)<(0.5*parents[0][pick1]+0.5*parents[0][pick2])){isthere = pick1; pick1 = pick2; pick2 = isthere;}
offspring[2][j]= parents[2][pick1];
offspring[2][j+N2] = parents[2][pick2];

}   // end for all individuals

/// offspring become parents in the next generation

for(i=0; i<NLs; i++){
for(j=0; j<Np0; j++){
parents[i][j]=offspring[i][j];
}}

		/// freq and hets at target
			for(j=0; j < Np0; j++) if(parents[2][j]==1.0)freq+=1;
			ps = (long double)freq/(long double)Np0;
			Het += 2.0*ps*(1.0-ps);
		/// freq and hets at pla modifier
			for(j=0; j < Np0; j++) if(parents[1][j]==1.0)plafreq+=1;
			Plas = (long double)plafreq/(long double)Np0;
			Phet += 2.0*Plas*(1.0-Plas);
		/// freqs and hets at recombination
			for(ti = 0; ti < 51; ti ++) recvec[ti]=0;
			for(j=0; j < Np0; j++){
			ti =(long)(parents[0][j]*100.00);
			recvec[ti]+=1;
			}
			for(ti = 0; ti <51; ti++){
			hetrec += ((double)recvec[ti]/(double)Np0)*(1.0-(double)recvec[ti]/(double)Np0);
			}


		/// record rec freqs right after burnin 
			if(k >= BURNIN){
			for(ti = 0; ti < 51; ti ++) recvec[ti]=0;
			        for(j=0; j < Np0; j++){
        			ti =(long)(parents[0][j]*100.00);
        			recvec[ti]+=1;
        		}

			for(ti = 0; ti < 51; ti++){sumrec[ti]+=(double)recvec[ti]/(double)Np0;}
			}

			if(k == (BURNIN + 3*PHASE -1)){
			for(ti = 0; ti < 51; ti++){fprintf(out2, "%d %f \n", ti, sumrec[ti]/(double)(3*PHASE));}}
	
}// end time

}// end reps 
	/// final freqs at the rec modifier
	for(ti = 0; ti < 51; ti++){
	fprintf(out2, "%d %f \n", 51+ti, sumrec[ti]/(double)(TIME-BURNIN)) ;
	}

fprintf(out2, "hetsR-P-T %Lf %Lf %Lf\n", hetrec/(long double)REPS, Phet/(long double)REPS, Het/(long double)REPS);

fclose(par);
fclose(out2);
fclose(inp);
return 0;
}  // end main
	

long double poidev(long double xm, long *idum)
{
   long double gammln(long double xx);
   long double ran1(long *idum);
   static long double sq,alxm,g,oldm=(-1.0);
   long double em,t,y;

   if (xm < 12.0) {
     if (xm != oldm) {
       oldm=xm;
       g=exp(-xm);
     }
     em = -1;
     t=1.0;
     do {
       ++em;
       t *= ran1(idum);
     } while (t > g);
   } else {
     if (xm != oldm) {
       oldm=xm;
      sq=sqrt(2.0*xm);
       alxm=log(xm);
       g=xm*alxm-gammln(xm+1.0);
     }
     do {
       do {
        y=tan(PI*ran1(idum));
         em=sq*y+xm;
       } while (em < 0.0);
       em=floor(em);
       t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
     } while (ran1(idum) > t);
   }
return em;
}



long double gammln(long double xx)
{
        long double x,y,tmp,ser;
        static long double cof[6]={76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
                0.1208650973866179e-2,-0.5395239384953e-5};
        int j;

        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
}


long double ran1(long *idum)
{
	long long j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	long double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

