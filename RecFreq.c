///// freq based evolution of rec2new recombianation rate introduced to a monomorphic population of rec2 rate, when a recombination modifier 
///// is recombining with the selected sequence with the rate rec1. Please see (and edit) "parameters" below.  The "output" file
///// reports the frequency of ancestral rate at  and at equilibrium. 

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

double bnldev(double pp, long n, long *idum);
double gammln(long double xx);
double ran1(long *idum);
double poidev(double xm, long *idum);

int main(void){

long REPS2, Np0, TIME, seed;    int PHASE;        double REC, C, swc, mu, rec1, rec2, rec2new, plas; 

FILE * inp, * out2;        out2 = fopen("output", "w");      inp = fopen("trajectory", "r");    seed = -1022;

////////////////////////// paramethers ///////////////////////////////////////////////////

PHASE = 20;                	// lenght oscillating cycle
rec2 = 0.0;                	// ancesteral rec rate
rec2new = 0.125;          	// invater rec rate
Np0 = 100000;             	// population size
TIME =  110*Np0;                // simulation time 
REPS2 = 100000;                	// number replicate simulation runs
mu =  0.1/(double)Np0;       	// mutation rate = numerator mutants per generation
rec1 = 0.5;                  	// recombiantion rate between rec modifier and sequence
plas = 1.0;               	// the effect of plasticity

/////////////////////////////////////////////////////////////////////////////////////////////

long mutants, N1, cc, k;  int y, z, m;
double  muts = mu*Np0, RF1, RF2, r1, r2, r3, r4, avefit, tra[3*PHASE], f1a, s2a, t3a, cond, aaa, ada, daa, dda, aad, add, dad, ddd, faaa, fada, fdaa, fdda, faad, fadd, fdad, fddd, paaa, pada, pdaa, pdda, paad, padd, pdad, pddd, fd, fa, avef1a = 0.0, avecf1a = 0.0;

for (k = 0; k < (3*PHASE); k++){
fscanf(inp, "%lf", &tra[k]);
}
swc = tra[PHASE/4];

r1 = rec1; r2 = rec2; r3= 0.5*rec2+0.5*rec2new; r4 = rec2new; RF1 = 1.0-r1-r3+2*r1*r3; RF2 = r1*r3-r1-r3;

for(cc =0 ; cc < REPS2; cc++){                // //////////////REPLICATOINS START///////////

aaa = 1.00; daa = 0.0; ada = 0.0; dda = 0.0; aad = 0.0; dad = 0.0; add = 0.0; ddd = 0.0; m = (int)(ran1(&seed)*PHASE);

        for( k = 0 ; k< TIME ; k++){
if (m == PHASE) m = 0;
////////////////////////////////////MUTATION 

	///////////////////////////////////////////// Rec locus

	mutants = poidev((muts),&seed);
	for(y = 0; y < mutants; y++){
	cond = ran1(&seed);
	if(cond<aaa){aaa-=1.00/(double)Np0; daa +=1.00/(double)Np0;}
	else if(cond < aaa+ada){ada -=1.00/(double)Np0; dda+=1.00/(double)Np0;}
	else if(cond < aaa+ada+aad){aad -=1.00/(double)Np0; dad+=1.00/(double)Np0;}
	else if(cond < aaa+ada+aad+add){add -=1.00/(double)Np0; ddd+=1.00/(double)Np0;}
	else if(cond < aaa+ada+aad+add+daa){daa -=1.00/(double)Np0; aaa+=1.00/(double)Np0; }
	else if(cond < aaa+ada+aad+add+daa+dad){dad -=1.00/(double)Np0; aad+=1.00/(double)Np0; }
	else if(cond < aaa+ada+aad+add+daa+dad+dda){dda -=1.00/(double)Np0; ada+=1.00/(double)Np0; }
	else {ddd -=1.00/(double)Np0; add+=1.00/(double)Np0; }
	}

	////////////////////////////////////////////target locus

	mutants = poidev((muts),&seed);
	for(y = 0; y < mutants; y++){
	cond = ran1(&seed);
	if(cond<aaa){aaa-=1.00/(double)Np0; aad+=1.00/(double)Np0; }
	else if(cond < aaa+ada){ada -=1.00/(double)Np0;add+=1.00/(double)Np0; }
	else if(cond < aaa+ada+aad){aad -=1.00/(double)Np0;aaa+=1.00/(double)Np0; }
	else if(cond < aaa+ada+aad+add){add -=1.00/(double)Np0;ada+=1.00/(double)Np0; }
	else if(cond < aaa+ada+aad+add+daa){daa -=1.00/(double)Np0; dad+=1.00/(double)Np0;}
	else if(cond < aaa+ada+aad+add+daa+dda){dda -=1.00/(double)Np0; ddd+=1.00/(double)Np0;}
	else if(cond < aaa+ada+aad+add+daa+dda+dad){dad -=1.00/(double)Np0; daa+=1.00/(double)Np0;}
	else {ddd-=1.00/(double)Np0; dda+=1.00/(double)Np0;}
	}

	///////////////////////////////////////////plasticity locus
	mutants = poidev((muts),&seed);
	for(y = 0; y < mutants; y++){
	cond = ran1(&seed);
	if(cond<aaa){aaa-=1.00/(double)Np0; ada+=1.00/(double)Np0;}
	else if(cond < aaa+ada){ada -=1.00/(double)Np0;aaa+=1.00/(double)Np0; }
	else if(cond < aaa+ada+aad){aad -=1.00/(double)Np0;add+=1.00/(double)Np0; }
	else if(cond < aaa+ada+aad+add){add -=1.00/(double)Np0;aad+=1.00/(double)Np0; }
	else if(cond < aaa+ada+aad+add+daa){daa -=1.00/(double)Np0;dda+=1.00/(double)Np0; }
	else if(cond < aaa+ada+aad+add+daa+dda){dda -=1.00/(double)Np0;daa+=1.00/(double)Np0; }
	else if(cond < aaa+ada+aad+add+daa+dda+dad){dad -=1.00/(double)Np0;ddd+=1.00/(double)Np0; }
	else {ddd-=1.00/(double)Np0; dad+=1.00/(double)Np0;}
	}

f1a = aaa+ada+aad+add;
s2a = aaa+daa+aad+dad;
t3a = aaa+ada+daa+dda;

if(f1a*(1.0-f1a)+s2a*(1.0-s2a)+t3a*(1.0-t3a)==0)goto WCONT;

	//////////////////////////////////SELECTION   


	fa = 1.0 - tra[m];
	fd = 1.0 + tra[m];

	faaa = fa; fada = 1.0-(1.0-plas)*tra[m]; fdaa = fa; fdda =1.0-(1.0-plas)*tra[m]; faad = fd; fadd = 1.0+(1.0-plas)*tra[m]; fdad = fd; fddd = 1.0+(1.0-plas)*tra[m];   

	avefit = aaa*faaa + ada*fada + daa *fdaa + dda*fdda + aad*faad   + add*fadd  + dad*fdad  + ddd *fddd;
	
	aaa = aaa*faaa/avefit;
	ada = ada*fada/avefit;
	daa = daa*fdaa/avefit;
	dda = dda*fdda/avefit;
	aad = aad*faad/avefit;
	add = add*fadd/avefit;
	dad = dad*fdad/avefit;
        ddd = ddd*fddd/avefit;

	////////////////////////////RECOMBINATION

	paaa = aaa*(1.0-dad-ddd-r2*add-r1*dda+RF1*dad+(1.0-r1)*(1.0-r3)*ddd)+ada*(r2*aad+r1*daa+r1*r3*dad)+aad*((1.0-RF1)*daa+(1.0-r1)*r3*dda)+r1*(1.0-r3)*add*daa;

	pada = ada*(1.0-dad-ddd-r2*aad-r1*daa+(1.0-r1)*(1.0-r3)*dad+RF1*ddd)+ aaa*(r2*add+r1*dda+r1*r3*ddd)+r1*(1.0-r3)*aad*dda+add*((1.0-r1)*r3*daa+(1.0-RF1)*dda);

	paad = aad*(1.0-r1*ddd-daa+RF2*dda+RF1*daa-r2*ada)+aaa*(r2*add+(1.0-RF1)*dad+(1.0-r1)*r3*ddd)+ada*dad*r1*(1.0-r3)+add*(daa*r1*r3+r1*dad);

	padd = add*(1.0-r2*aaa+RF2*daa+(RF1-1.0)*dda-r1*dad)+aaa*ddd*r1*(1.0-r3)+ada*(aad*r2+dad*(1.0-r1)*r3+ddd*(1.0-RF1))+aad*(dda*r1*r3+ddd*r1);

	pdaa = daa*(1.0-r1*ada+RF2*add-r4*ddd+(RF1-1.0)*aad)+aaa*(r1*dda+(1.0-RF1)*dad+r1*(1.0-r3)*ddd) +ada*dad*(1.0-r1)*r3+aad*dda*r1*r3+dda*dad*r4; 

	pdda = dda*(1.0-r1*aaa+RF2*aad+(RF1-1.0)*add-r4*dad)+(1.0-r1)*r3*aaa*ddd+ada*(daa*r1+dad*r1*(1.0-r3)+ddd*(1.0-RF1))+daa*(add*r1*r3+ddd*r4);

	pdad = dad*(1.0+(RF1-1.0)*aaa+RF2*ada-r1*add-r4*dda)+aaa*ddd*r1*r3+aad*(daa*(1.0-RF1)+dda*r1*(1.0-r3)+ddd*r1)+daa*(add*(1.0-r1)*r3+ddd*r4);

	pddd = ddd*(1.0+RF2*aaa+(RF1-1.0)*ada-r1*aad-r4*daa)+add*(daa*r1*(1.0-r3)+dda*(1.0-RF1)+dad*r1)+dad*(dda*r4+ada*r1*r3)+aad*dda*(1.0-r1)*r3;


        ///////////// REPRODUCTION
	aaa = paaa;
	if(paaa<1.0)ada = pada/(1.0-paaa);
	if((paaa+pada)<1.0)aad = paad/(1.0-paaa-pada);
	if((paaa+pada+paad)<1.0)add = padd/(1.0-paaa-pada-paad);
	if((paaa+pada+paad+padd)<1.0)daa = pdaa/(1.0-paaa-pada-paad-padd);
	if((paaa+pada+paad+padd+pdaa)<1.0)dda = pdda/(1.0-paaa-pada-paad-padd-pdaa);
	if((paaa+pada+paad+padd+pdaa+pdda)<1.0)dad = pdad/(1.0-paaa-pada-paad-padd-pdaa-pdda);

	if(aaa>0.0)aaa = bnldev(aaa, Np0, &seed);
	N1 = Np0 - (long)aaa;

	if(N1>0)ada = bnldev(ada, N1, &seed);
	else{ada = 0.0;}
	N1 = N1 - (long)ada;

	if(N1>0) aad  = bnldev( aad, N1, &seed);
	else{aad=0.0;}
	N1 = N1 -(long)aad ;

	if(N1>0)add  = bnldev(add , N1, &seed);
	else{add = 0.0;}
	N1 = N1 - (long)add;

	if(N1>0) daa = bnldev( daa, N1, &seed);
	else{daa =0.0;}
	N1 = N1 - (long)daa;

	if(N1>0) dda = bnldev( dda, N1, &seed);
	else{dda = 0.0;}
	N1 = N1 - (long)dda;

	if(N1>0) dad = bnldev( dad, N1, &seed);
	else{dad = 0.0;}
	N1 = N1 - (long)dad;

	aaa = aaa/(double)Np0; ada = ada/(double)Np0; aad = aad/(double)Np0; add = add/(double)Np0; daa = daa/(double)Np0; dda = dda/(double)Np0; dad = dad/(double)Np0;

	ddd = 1.0 - aaa - ada - aad - add - daa - dda - dad;

/////////////// summaries

f1a = aaa+ada+aad+add;
if(k>100*Np0)avef1a += f1a;
if(k==100*Np0)avecf1a += f1a;

WCONT:;
m+=1;
	} // end k
}// end reps 

fprintf(out2, "P %d s %.3f N %ld mu %.10lf pla %lf  reps %ld  RecR %lf  RecPT %lf RecPTnew %lf freqAncRate %lf freqAncRateatBurnin %lf \n", PHASE, swc, Np0,mu, plas, REPS2, rec1, rec2, rec2new, avef1a/(double)(REPS2*(10*Np0-1)), avecf1a/(double)REPS2);

fclose(inp);
fclose(out2);
return 0;
}  // end main
	



double bnldev( double pp, long n, long *idum)
{

long long j;

static long long nold =(-1);

double am, em, g, angle, p, bnl, sq, t, y;
static double pold=(-1.0), pc, plog, pclog, en, oldg;

p=(pp <= 0.5 ? pp : 1.0-pp);

am = n*p;

if(n<25){
	bnl = 0.0;

	for(j=1; j<=n; j++)
		if(ran1(idum) < p) ++bnl;
}else if (am < 1.0) {
  	g=exp(-am);
  	t = 1.0;
	for (j=0; j<=n; j++){
		t *= ran1(idum);
		if (t < g) break;
}
bnl=(j<=n ? j:n);
}else {


if (n != nold) {
	en = n;
	oldg = gammln(en+1.0);
	nold = n;
}

if (p != pold) {
	pc = 1.0-p;
	plog = log(p);
	pclog = log(pc);
	pold = p;
}

sq = sqrt(2.0*am*pc);
do{
	do{
		angle = PI * ran1(idum);
		y = tan(angle);
		em = sq*y+am;
	} while (em < 0.0 || em >= (en + 1.0));
	
	em = floor(em);
	t = 1.2 * sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
}while (ran1(idum) > t);


bnl = em;

}

if (p != pp) bnl = n - bnl;
return bnl;

}


double gammln(long double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	long long j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}




double ran1(long *idum)
{
	long long j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

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

double poidev(double xm, long *idum)
{
   double gammln(long double xx);
   double ran1(long *idum);
   static double sq,alxm,g,oldm=(-1.0);
   double em,t,y;

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


