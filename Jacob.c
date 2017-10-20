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

int main(void){

long TIME; int PHASE; long double rec1, rec2, rec2new, plas, START, END, rp, rnewp;   

FILE * inp, * out2, *par;
out2 = fopen("output", "w"); inp = fopen("trajectory", "r"); par = fopen("pars1", "r");

/////////////////////////////////////////////////////paramethers 

fscanf(par, "%d", &PHASE);
fscanf(par, "%Lf", &START);
fscanf(par, "%Lf", &END);

TIME = 100000000000;   //// if no equilibrium is reached by TIME, "output" will report it 
rec1 = 0.5;            //// recombination between the rec. modifier and the pla-tg sequence
plas = 1.00;           //// relative effect of plasticity
rp = 0.99;             //// ancestral rec. rate
rnewp = 1.00 - rp;     //// invader rec. rate

/////////////////////////////////////////////////////

long k; 
int z, m, dojacob, jacobdone; 
long double  swc, a, p, RF1, RF2, r1, r2, r3, r4, avefit, tra[PHASE], f1a, s2a, t3a, cond1, aaa, ada, daa, dda, aad, add, dad, ddd, faaa, fada, eqaaa , eqada ,  eqdaa, eqdda, eqaad, eqadd, eqdad, eqddd, olaaa, olada, oldaa, oldda, olaad, oladd, oldad, olddd, fdaa, fdda, faad, fadd, fdad, fddd, paaa, pada, pdaa, pdda, paad, padd, pdad, pddd, fd, fa;
long double epsilon, e11, e12, e13, e14, e15, e16, e17, e21, e22, e23, e24, e25, e26, e27, e31, e32, e33, e34, e35, e36, e37, e41, e42, e43, e44, e45, e46, e47, e51, e52, e53, e54, e55, e56, e57, e61, e62, e63, e64, e65, e66, e67, e71, e72, e73, e74, e75, e76, e77, dx01aaa, dx01daa, dx01ada, dx01dda, dx01aad, dx01dad, dx01add, dx01ddd, dx02aaa, dx02daa, dx02ada, dx02dda, dx02aad, dx02dad, dx02add, dx02ddd, dx03aaa, dx03daa, dx03ada, dx03dda, dx03aad, dx03dad, dx03add, dx03ddd, dx04aaa, dx04daa, dx04ada, dx04dda, dx04aad, dx04dad, dx04add, dx04ddd, dx05aaa, dx05daa, dx05ada, dx05dda, dx05aad, dx05dad, dx05add, dx05ddd, dx06aaa, dx06daa, dx06ada, dx06dda, dx06aad, dx06dad, dx06add, dx06ddd, dx07aaa, dx07daa, dx07ada, dx07dda, dx07aad, dx07dad, dx07add, dx07ddd, dx08aaa, dx08daa, dx08ada, dx08dda, dx08aad, dx08dad, dx08add, dx08ddd, dx09aaa, dx09daa, dx09ada, dx09dda, dx09aad, dx09dad, dx09add, dx09ddd, dx10aaa, dx10daa, dx10ada, dx10dda, dx10aad, dx10dad, dx10add, dx10ddd, dx11aaa, dx11daa, dx11ada, dx11dda, dx11aad, dx11dad, dx11add, dx11ddd, dx12aaa, dx12daa, dx12ada, dx12dda, dx12aad, dx12dad, dx12add, dx12ddd, dx13aaa, dx13daa, dx13ada, dx13dda, dx13aad, dx13dad, dx13add, dx13ddd, dx14aaa, dx14daa, dx14ada, dx14dda, dx14aad, dx14dad, dx14add, dx14ddd;  

for (k = 0; k < PHASE; k++){
fscanf(inp, "%Lf", &tra[k]);
}

swc = tra[PHASE/4];

for(rec2 =START; rec2 <END; rec2 += 0.01){
	for(rec2new = 0.00; rec2new < 0.51; rec2new +=0.01){

if((rec2new - rec2)*(rec2new - rec2) < 0.00001) goto WCONT;

r1 = rec1; r2 = rec2; r3= 0.5*rec2+0.5*rec2new; r4 = rec2new; RF1 = 1.0-r1-r3+2*r1*r3; RF2 = r1*r3-r1-r3;

eqaaa = 0.0; eqada = 0.0;  eqdaa=0.0; eqdda=0.0; eqaad=0.0; eqadd=0.0; eqdad=0.0; eqddd=0.0;
 
		for(p=0.05;p<1.0;p+=0.1){
			for(a=0.05;a<1.0;a+=0.1){
//////starting haplotype freqs.

aaa = rp*p*a; 
daa = rnewp*p*a; 
ada = rp*(1.00-p)*a; 
dda = rnewp*(1.00-p)*a; 
aad = rp*p*(1.00-a); 
dad = rnewp*p*(1.00-a); 
add = rp*(1.00-p)*(1.00-a); 

m =0; dojacob =0; jacobdone = 0; olaaa=1.0; olada=1.0; oldaa=1.0; oldda=1.0; olaad=1.0; oladd=1.0; oldad=1.0; olddd=1.0;

for(k=0; k<TIME; k++){  // evolve and identify new equilibrium
// selection
	if (m == PHASE) m = 0;
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
        m+=1;
// recombination
        paaa = aaa*(1.0-dad-ddd-r2*add-r1*dda+RF1*dad+(1.0-r1)*(1.0-r3)*ddd)+ada*(r2*aad+r1*daa+r1*r3*dad)+aad*((1.0-RF1)*daa+(1.0-r1)*r3*dda)+r1*(1.0-r3)*add*daa;
	pada = ada*(1.0-dad-ddd-r2*aad-r1*daa+(1.0-r1)*(1.0-r3)*dad+RF1*ddd)+ aaa*(r2*add+r1*dda+r1*r3*ddd)+r1*(1.0-r3)*aad*dda+add*((1.0-r1)*r3*daa+(1.0-RF1)*dda);
        paad = aad*(1.0-r1*ddd-daa+RF2*dda+RF1*daa-r2*ada)+aaa*(r2*add+(1.0-RF1)*dad+(1.0-r1)*r3*ddd)+ada*dad*r1*(1.0-r3)+add*(daa*r1*r3+r1*dad);
        padd = add*(1.0-r2*aaa+RF2*daa+(RF1-1.0)*dda-r1*dad)+aaa*ddd*r1*(1.0-r3)+ada*(aad*r2+dad*(1.0-r1)*r3+ddd*(1.0-RF1))+aad*(dda*r1*r3+ddd*r1);
        pdaa = daa*(1.0-r1*ada+RF2*add-r4*ddd+(RF1-1.0)*aad)+aaa*(r1*dda+(1.0-RF1)*dad+r1*(1.0-r3)*ddd) +ada*dad*(1.0-r1)*r3+aad*dda*r1*r3+dda*dad*r4;
        pdda = dda*(1.0-r1*aaa+RF2*aad+(RF1-1.0)*add-r4*dad)+(1.0-r1)*r3*aaa*ddd+ada*(daa*r1+dad*r1*(1.0-r3)+ddd*(1.0-RF1))+daa*(add*r1*r3+ddd*r4);
        pdad = dad*(1.0+(RF1-1.0)*aaa+RF2*ada-r1*add-r4*dda)+aaa*ddd*r1*r3+aad*(daa*(1.0-RF1)+dda*r1*(1.0-r3)+ddd*r1)+daa*(add*(1.0-r1)*r3+ddd*r4);
        pddd = ddd*(1.0+RF2*aaa+(RF1-1.0)*ada-r1*aad-r4*daa)+add*(daa*r1*(1.0-r3)+dda*(1.0-RF1)+dad*r1)+dad*(dda*r4+ada*r1*r3)+aad*dda*(1.0-r1)*r3;
	aaa = paaa; ada = pada; daa = pdaa; dda = pdda; aad = paad; add = padd; dad = pdad; ddd = pddd;

// is it equilibrium?

	if(m==PHASE && dojacob ==0){
	
	if((olaaa-aaa) < 0.000000001 && (olada-ada) < 0.000000001 && (oldaa-daa) < 0.000000001 && (oldda-dda) < 0.000000001 && (olaad-aad) < 0.000000001 && (oladd-add) < 0.000000001 && (oldad-dad) < 0.000000001 && (olddd-ddd)< 0.000000001){
		f1a = aaa+ada+aad+add; s2a = aaa+daa+aad+dad; t3a = aaa+ada+daa+dda;
	
		if((s2a*(1.0-s2a))< 0.01 || (t3a*(1.0-t3a)) < 0.01)jacobdone=1;
		
                else if((aaa-eqaaa) < 0.01 && (ada- eqada) < 0.01 && (daa-eqdaa) < 0.01 && (dda-eqdda) < 0.01 && (aad-eqaad) < 0.01 && (add-eqadd) < 0.01 && (dad-eqdad) < 0.01 && (ddd-eqddd)<0.01){
		//fprintf(out2, "#");
		jacobdone = 1;
		}

		else{
		eqaaa = aaa; eqada = ada;  eqdaa=daa; eqdda=dda; eqaad=aad; eqadd=add; eqdad=dad; eqddd=ddd;	
		fprintf(out2, "\nP %d Smax %Lf      r-ancestral %Lf r-invader %Lf  equilibrium freqs at PHASE: Rec_Anc %.2Lf Pla_Anc %.2Lf Tg_Anc %.2Lf, R code:  \n\n", PHASE, swc, rec2, rec2new, f1a, s2a, t3a);  
		dojacob =1;
	}}
	else{
	olaaa=aaa; olada=ada; oldaa=daa; oldda=dda; olaad=aad; oladd=add; oldad=dad; olddd=ddd;
	}}

// jacobian calc.

if(dojacob==1){
for(z=0; z<3; z++){
if(z==0)epsilon = 0.0000000001L; if(z==1)epsilon= 0.00000000001L; if(z==2) epsilon = 0.000000000001L;
for(m=0; m<PHASE; m++){

        fa = 1.0 - tra[m];
        fd = 1.0 + tra[m];
        faaa = fa; fada = 1.0-(1.0-plas)*tra[m]; fdaa = fa; fdda =1.0-(1.0-plas)*tra[m]; faad = fd; fadd = 1.0+(1.0-plas)*tra[m]; fdad = fd; fddd = 1.0+(1.0-plas)*tra[m];

/////////////////// first column
if(m==0){
dx01aaa = eqaaa + epsilon;
dx01daa = eqdaa;
dx01ada = eqada;
dx01dda = eqdda;
dx01aad = eqaad;
dx01dad = eqdad;
dx01add = eqadd;
dx01ddd = 1.00 - dx01aaa - dx01daa - dx01ada - dx01dda - dx01aad - dx01dad - dx01add;
}

        avefit = dx01aaa*faaa + dx01ada*fada + dx01daa *fdaa + dx01dda*fdda + dx01aad*faad   + dx01add*fadd  + dx01dad*fdad  + dx01ddd *fddd;

        dx01aaa = dx01aaa*faaa/avefit;
        dx01ada = dx01ada*fada/avefit;
        dx01daa = dx01daa*fdaa/avefit;
        dx01dda = dx01dda*fdda/avefit;
        dx01aad = dx01aad*faad/avefit;
        dx01add = dx01add*fadd/avefit;
        dx01dad = dx01dad*fdad/avefit;
        dx01ddd = dx01ddd*fddd/avefit;

        paaa = dx01aaa*(1.0-dx01dad-dx01ddd-r2*dx01add-r1*dx01dda+RF1*dx01dad+(1.0-r1)*(1.0-r3)*dx01ddd)+dx01ada*(r2*dx01aad+r1*dx01daa+r1*r3*dx01dad)+dx01aad*((1.0-RF1)*dx01daa+(1.0-r1)*r3*dx01dda)+r1*(1.0-r3)*dx01add*dx01daa;
        pada = dx01ada*(1.0-dx01dad-dx01ddd-r2*dx01aad-r1*dx01daa+(1.0-r1)*(1.0-r3)*dx01dad+RF1*dx01ddd)+ dx01aaa*(r2*dx01add+r1*dx01dda+r1*r3*dx01ddd)+r1*(1.0-r3)*dx01aad*dx01dda+dx01add*((1.0-r1)*r3*dx01daa+(1.0-RF1)*dx01dda);
        paad = dx01aad*(1.0-r1*dx01ddd-dx01daa+RF2*dx01dda+RF1*dx01daa-r2*dx01ada)+dx01aaa*(r2*dx01add+(1.0-RF1)*dx01dad+(1.0-r1)*r3*dx01ddd)+dx01ada*dx01dad*r1*(1.0-r3)+dx01add*(dx01daa*r1*r3+r1*dx01dad);
        padd = dx01add*(1.0-r2*dx01aaa+RF2*dx01daa+(RF1-1.0)*dx01dda-r1*dx01dad)+dx01aaa*dx01ddd*r1*(1.0-r3)+dx01ada*(dx01aad*r2+dx01dad*(1.0-r1)*r3+dx01ddd*(1.0-RF1))+dx01aad*(dx01dda*r1*r3+dx01ddd*r1);
        pdaa = dx01daa*(1.0-r1*dx01ada+RF2*dx01add-r4*dx01ddd+(RF1-1.0)*dx01aad)+dx01aaa*(r1*dx01dda+(1.0-RF1)*dx01dad+r1*(1.0-r3)*dx01ddd) +dx01ada*dx01dad*(1.0-r1)*r3+dx01aad*dx01dda*r1*r3+dx01dda*dx01dad*r4;
        pdda = dx01dda*(1.0-r1*dx01aaa+RF2*dx01aad+(RF1-1.0)*dx01add-r4*dx01dad)+(1.0-r1)*r3*dx01aaa*dx01ddd+dx01ada*(dx01daa*r1+dx01dad*r1*(1.0-r3)+dx01ddd*(1.0-RF1))+dx01daa*(dx01add*r1*r3+dx01ddd*r4);
        pdad = dx01dad*(1.0+(RF1-1.0)*dx01aaa+RF2*dx01ada-r1*dx01add-r4*dx01dda)+dx01aaa*dx01ddd*r1*r3+dx01aad*(dx01daa*(1.0-RF1)+dx01dda*r1*(1.0-r3)+dx01ddd*r1)+dx01daa*(dx01add*(1.0-r1)*r3+dx01ddd*r4);
        dx01aaa = paaa; dx01ada = pada; dx01daa = pdaa; dx01dda = pdda; dx01aad = paad; dx01add = padd; dx01dad = pdad; dx01ddd = 1.00 - dx01aaa - dx01daa - dx01ada - dx01dda - dx01aad - dx01dad - dx01add;

if(m==0){
dx02aaa = eqaaa - epsilon;
dx02daa = eqdaa;
dx02ada = eqada;
dx02dda = eqdda;
dx02aad = eqaad;
dx02dad = eqdad;
dx02add = eqadd;
dx02ddd = 1.00 - dx02aaa - dx02daa - dx02ada - dx02dda - dx02aad - dx02dad - dx02add;
}

        avefit = dx02aaa*faaa + dx02ada*fada + dx02daa *fdaa + dx02dda*fdda + dx02aad*faad   + dx02add*fadd  + dx02dad*fdad  + dx02ddd *fddd;

        dx02aaa = dx02aaa*faaa/avefit;
        dx02ada = dx02ada*fada/avefit;
        dx02daa = dx02daa*fdaa/avefit;
        dx02dda = dx02dda*fdda/avefit;
        dx02aad = dx02aad*faad/avefit;
        dx02add = dx02add*fadd/avefit;
        dx02dad = dx02dad*fdad/avefit;
        dx02ddd = dx02ddd*fddd/avefit;

        paaa = dx02aaa*(1.0-dx02dad-dx02ddd-r2*dx02add-r1*dx02dda+RF1*dx02dad+(1.0-r1)*(1.0-r3)*dx02ddd)+dx02ada*(r2*dx02aad+r1*dx02daa+r1*r3*dx02dad)+dx02aad*((1.0-RF1)*dx02daa+(1.0-r1)*r3*dx02dda)+r1*(1.0-r3)*dx02add*dx02daa;
        pada = dx02ada*(1.0-dx02dad-dx02ddd-r2*dx02aad-r1*dx02daa+(1.0-r1)*(1.0-r3)*dx02dad+RF1*dx02ddd)+ dx02aaa*(r2*dx02add+r1*dx02dda+r1*r3*dx02ddd)+r1*(1.0-r3)*dx02aad*dx02dda+dx02add*((1.0-r1)*r3*dx02daa+(1.0-RF1)*dx02dda);
        paad = dx02aad*(1.0-r1*dx02ddd-dx02daa+RF2*dx02dda+RF1*dx02daa-r2*dx02ada)+dx02aaa*(r2*dx02add+(1.0-RF1)*dx02dad+(1.0-r1)*r3*dx02ddd)+dx02ada*dx02dad*r1*(1.0-r3)+dx02add*(dx02daa*r1*r3+r1*dx02dad);
        padd = dx02add*(1.0-r2*dx02aaa+RF2*dx02daa+(RF1-1.0)*dx02dda-r1*dx02dad)+dx02aaa*dx02ddd*r1*(1.0-r3)+dx02ada*(dx02aad*r2+dx02dad*(1.0-r1)*r3+dx02ddd*(1.0-RF1))+dx02aad*(dx02dda*r1*r3+dx02ddd*r1);
        pdaa = dx02daa*(1.0-r1*dx02ada+RF2*dx02add-r4*dx02ddd+(RF1-1.0)*dx02aad)+dx02aaa*(r1*dx02dda+(1.0-RF1)*dx02dad+r1*(1.0-r3)*dx02ddd) +dx02ada*dx02dad*(1.0-r1)*r3+dx02aad*dx02dda*r1*r3+dx02dda*dx02dad*r4;
        pdda = dx02dda*(1.0-r1*dx02aaa+RF2*dx02aad+(RF1-1.0)*dx02add-r4*dx02dad)+(1.0-r1)*r3*dx02aaa*dx02ddd+dx02ada*(dx02daa*r1+dx02dad*r1*(1.0-r3)+dx02ddd*(1.0-RF1))+dx02daa*(dx02add*r1*r3+dx02ddd*r4);
        pdad = dx02dad*(1.0+(RF1-1.0)*dx02aaa+RF2*dx02ada-r1*dx02add-r4*dx02dda)+dx02aaa*dx02ddd*r1*r3+dx02aad*(dx02daa*(1.0-RF1)+dx02dda*r1*(1.0-r3)+dx02ddd*r1)+dx02daa*(dx02add*(1.0-r1)*r3+dx02ddd*r4);
        dx02aaa = paaa; dx02ada = pada; dx02daa = pdaa; dx02dda = pdda; dx02aad = paad; dx02add = padd; dx02dad = pdad; dx02ddd =1.00 - dx02aaa - dx02daa - dx02ada - dx02dda - dx02aad - dx02dad - dx02add;

////// second column
if(m==0){
dx03aaa = eqaaa;
dx03daa = eqdaa + epsilon;
dx03ada = eqada;
dx03dda = eqdda;
dx03aad = eqaad;
dx03dad = eqdad;
dx03add = eqadd;
dx03ddd = 1.00 - dx03aaa - dx03daa - dx03ada - dx03dda - dx03aad - dx03dad - dx03add;
}
        avefit = dx03aaa*faaa + dx03ada*fada + dx03daa *fdaa + dx03dda*fdda + dx03aad*faad   + dx03add*fadd  + dx03dad*fdad  + dx03ddd *fddd;

        dx03aaa = dx03aaa*faaa/avefit;
        dx03ada = dx03ada*fada/avefit;
        dx03daa = dx03daa*fdaa/avefit;
        dx03dda = dx03dda*fdda/avefit;
        dx03aad = dx03aad*faad/avefit;
        dx03add = dx03add*fadd/avefit;
        dx03dad = dx03dad*fdad/avefit;
        dx03ddd = dx03ddd*fddd/avefit;

        paaa = dx03aaa*(1.0-dx03dad-dx03ddd-r2*dx03add-r1*dx03dda+RF1*dx03dad+(1.0-r1)*(1.0-r3)*dx03ddd)+dx03ada*(r2*dx03aad+r1*dx03daa+r1*r3*dx03dad)+dx03aad*((1.0-RF1)*dx03daa+(1.0-r1)*r3*dx03dda)+r1*(1.0-r3)*dx03add*dx03daa;
        pada = dx03ada*(1.0-dx03dad-dx03ddd-r2*dx03aad-r1*dx03daa+(1.0-r1)*(1.0-r3)*dx03dad+RF1*dx03ddd)+ dx03aaa*(r2*dx03add+r1*dx03dda+r1*r3*dx03ddd)+r1*(1.0-r3)*dx03aad*dx03dda+dx03add*((1.0-r1)*r3*dx03daa+(1.0-RF1)*dx03dda);
        paad = dx03aad*(1.0-r1*dx03ddd-dx03daa+RF2*dx03dda+RF1*dx03daa-r2*dx03ada)+dx03aaa*(r2*dx03add+(1.0-RF1)*dx03dad+(1.0-r1)*r3*dx03ddd)+dx03ada*dx03dad*r1*(1.0-r3)+dx03add*(dx03daa*r1*r3+r1*dx03dad);
        padd = dx03add*(1.0-r2*dx03aaa+RF2*dx03daa+(RF1-1.0)*dx03dda-r1*dx03dad)+dx03aaa*dx03ddd*r1*(1.0-r3)+dx03ada*(dx03aad*r2+dx03dad*(1.0-r1)*r3+dx03ddd*(1.0-RF1))+dx03aad*(dx03dda*r1*r3+dx03ddd*r1);
        pdaa = dx03daa*(1.0-r1*dx03ada+RF2*dx03add-r4*dx03ddd+(RF1-1.0)*dx03aad)+dx03aaa*(r1*dx03dda+(1.0-RF1)*dx03dad+r1*(1.0-r3)*dx03ddd) +dx03ada*dx03dad*(1.0-r1)*r3+dx03aad*dx03dda*r1*r3+dx03dda*dx03dad*r4;
        pdda = dx03dda*(1.0-r1*dx03aaa+RF2*dx03aad+(RF1-1.0)*dx03add-r4*dx03dad)+(1.0-r1)*r3*dx03aaa*dx03ddd+dx03ada*(dx03daa*r1+dx03dad*r1*(1.0-r3)+dx03ddd*(1.0-RF1))+dx03daa*(dx03add*r1*r3+dx03ddd*r4);
        pdad = dx03dad*(1.0+(RF1-1.0)*dx03aaa+RF2*dx03ada-r1*dx03add-r4*dx03dda)+dx03aaa*dx03ddd*r1*r3+dx03aad*(dx03daa*(1.0-RF1)+dx03dda*r1*(1.0-r3)+dx03ddd*r1)+dx03daa*(dx03add*(1.0-r1)*r3+dx03ddd*r4);
        dx03aaa = paaa; dx03ada = pada; dx03daa = pdaa; dx03dda = pdda; dx03aad = paad; dx03add = padd; dx03dad = pdad; dx03ddd = 1.00 - dx03aaa - dx03daa - dx03ada - dx03dda - dx03aad - dx03dad - dx03add;

if(m==0){
dx04aaa = eqaaa;
dx04daa = eqdaa - epsilon;
dx04ada = eqada;
dx04dda = eqdda;
dx04aad = eqaad;
dx04dad = eqdad;
dx04add = eqadd;
dx04ddd = 1.00 - dx04aaa - dx04daa - dx04ada - dx04dda - dx04aad - dx04dad - dx04add;
}
        avefit = dx04aaa*faaa + dx04ada*fada + dx04daa *fdaa + dx04dda*fdda + dx04aad*faad   + dx04add*fadd  + dx04dad*fdad  + dx04ddd *fddd;

        dx04aaa = dx04aaa*faaa/avefit;
        dx04ada = dx04ada*fada/avefit;
        dx04daa = dx04daa*fdaa/avefit;
        dx04dda = dx04dda*fdda/avefit;
        dx04aad = dx04aad*faad/avefit;
        dx04add = dx04add*fadd/avefit;
        dx04dad = dx04dad*fdad/avefit;
        dx04ddd = dx04ddd*fddd/avefit;

        paaa = dx04aaa*(1.0-dx04dad-dx04ddd-r2*dx04add-r1*dx04dda+RF1*dx04dad+(1.0-r1)*(1.0-r3)*dx04ddd)+dx04ada*(r2*dx04aad+r1*dx04daa+r1*r3*dx04dad)+dx04aad*((1.0-RF1)*dx04daa+(1.0-r1)*r3*dx04dda)+r1*(1.0-r3)*dx04add*dx04daa;
        pada = dx04ada*(1.0-dx04dad-dx04ddd-r2*dx04aad-r1*dx04daa+(1.0-r1)*(1.0-r3)*dx04dad+RF1*dx04ddd)+ dx04aaa*(r2*dx04add+r1*dx04dda+r1*r3*dx04ddd)+r1*(1.0-r3)*dx04aad*dx04dda+dx04add*((1.0-r1)*r3*dx04daa+(1.0-RF1)*dx04dda);
        paad = dx04aad*(1.0-r1*dx04ddd-dx04daa+RF2*dx04dda+RF1*dx04daa-r2*dx04ada)+dx04aaa*(r2*dx04add+(1.0-RF1)*dx04dad+(1.0-r1)*r3*dx04ddd)+dx04ada*dx04dad*r1*(1.0-r3)+dx04add*(dx04daa*r1*r3+r1*dx04dad);
        padd = dx04add*(1.0-r2*dx04aaa+RF2*dx04daa+(RF1-1.0)*dx04dda-r1*dx04dad)+dx04aaa*dx04ddd*r1*(1.0-r3)+dx04ada*(dx04aad*r2+dx04dad*(1.0-r1)*r3+dx04ddd*(1.0-RF1))+dx04aad*(dx04dda*r1*r3+dx04ddd*r1);
        pdaa = dx04daa*(1.0-r1*dx04ada+RF2*dx04add-r4*dx04ddd+(RF1-1.0)*dx04aad)+dx04aaa*(r1*dx04dda+(1.0-RF1)*dx04dad+r1*(1.0-r3)*dx04ddd) +dx04ada*dx04dad*(1.0-r1)*r3+dx04aad*dx04dda*r1*r3+dx04dda*dx04dad*r4;
        pdda = dx04dda*(1.0-r1*dx04aaa+RF2*dx04aad+(RF1-1.0)*dx04add-r4*dx04dad)+(1.0-r1)*r3*dx04aaa*dx04ddd+dx04ada*(dx04daa*r1+dx04dad*r1*(1.0-r3)+dx04ddd*(1.0-RF1))+dx04daa*(dx04add*r1*r3+dx04ddd*r4);
        pdad = dx04dad*(1.0+(RF1-1.0)*dx04aaa+RF2*dx04ada-r1*dx04add-r4*dx04dda)+dx04aaa*dx04ddd*r1*r3+dx04aad*(dx04daa*(1.0-RF1)+dx04dda*r1*(1.0-r3)+dx04ddd*r1)+dx04daa*(dx04add*(1.0-r1)*r3+dx04ddd*r4);
        dx04aaa = paaa; dx04ada = pada; dx04daa = pdaa; dx04dda = pdda; dx04aad = paad; dx04add = padd; dx04dad = pdad; dx04ddd = 1.00 - dx04aaa - dx04daa - dx04ada - dx04dda - dx04aad - dx04dad - dx04add;

/////////////third column

if(m==0){
dx05aaa = eqaaa;
dx05daa = eqdaa;
dx05ada = eqada + epsilon;
dx05dda = eqdda;
dx05aad = eqaad;
dx05dad = eqdad;
dx05add = eqadd;
dx05ddd = 1.00 - dx05aaa - dx05daa - dx05ada - dx05dda - dx05aad - dx05dad - dx05add;
}

        avefit = dx05aaa*faaa + dx05ada*fada + dx05daa *fdaa + dx05dda*fdda + dx05aad*faad   + dx05add*fadd  + dx05dad*fdad  + dx05ddd *fddd;

        dx05aaa = dx05aaa*faaa/avefit;
        dx05ada = dx05ada*fada/avefit;
        dx05daa = dx05daa*fdaa/avefit;
        dx05dda = dx05dda*fdda/avefit;
        dx05aad = dx05aad*faad/avefit;
        dx05add = dx05add*fadd/avefit;
        dx05dad = dx05dad*fdad/avefit;
        dx05ddd = dx05ddd*fddd/avefit;

        paaa = dx05aaa*(1.0-dx05dad-dx05ddd-r2*dx05add-r1*dx05dda+RF1*dx05dad+(1.0-r1)*(1.0-r3)*dx05ddd)+dx05ada*(r2*dx05aad+r1*dx05daa+r1*r3*dx05dad)+dx05aad*((1.0-RF1)*dx05daa+(1.0-r1)*r3*dx05dda)+r1*(1.0-r3)*dx05add*dx05daa;
        pada = dx05ada*(1.0-dx05dad-dx05ddd-r2*dx05aad-r1*dx05daa+(1.0-r1)*(1.0-r3)*dx05dad+RF1*dx05ddd)+ dx05aaa*(r2*dx05add+r1*dx05dda+r1*r3*dx05ddd)+r1*(1.0-r3)*dx05aad*dx05dda+dx05add*((1.0-r1)*r3*dx05daa+(1.0-RF1)*dx05dda);
        paad = dx05aad*(1.0-r1*dx05ddd-dx05daa+RF2*dx05dda+RF1*dx05daa-r2*dx05ada)+dx05aaa*(r2*dx05add+(1.0-RF1)*dx05dad+(1.0-r1)*r3*dx05ddd)+dx05ada*dx05dad*r1*(1.0-r3)+dx05add*(dx05daa*r1*r3+r1*dx05dad);
        padd = dx05add*(1.0-r2*dx05aaa+RF2*dx05daa+(RF1-1.0)*dx05dda-r1*dx05dad)+dx05aaa*dx05ddd*r1*(1.0-r3)+dx05ada*(dx05aad*r2+dx05dad*(1.0-r1)*r3+dx05ddd*(1.0-RF1))+dx05aad*(dx05dda*r1*r3+dx05ddd*r1);
        pdaa = dx05daa*(1.0-r1*dx05ada+RF2*dx05add-r4*dx05ddd+(RF1-1.0)*dx05aad)+dx05aaa*(r1*dx05dda+(1.0-RF1)*dx05dad+r1*(1.0-r3)*dx05ddd) +dx05ada*dx05dad*(1.0-r1)*r3+dx05aad*dx05dda*r1*r3+dx05dda*dx05dad*r4;
        pdda = dx05dda*(1.0-r1*dx05aaa+RF2*dx05aad+(RF1-1.0)*dx05add-r4*dx05dad)+(1.0-r1)*r3*dx05aaa*dx05ddd+dx05ada*(dx05daa*r1+dx05dad*r1*(1.0-r3)+dx05ddd*(1.0-RF1))+dx05daa*(dx05add*r1*r3+dx05ddd*r4);
        pdad = dx05dad*(1.0+(RF1-1.0)*dx05aaa+RF2*dx05ada-r1*dx05add-r4*dx05dda)+dx05aaa*dx05ddd*r1*r3+dx05aad*(dx05daa*(1.0-RF1)+dx05dda*r1*(1.0-r3)+dx05ddd*r1)+dx05daa*(dx05add*(1.0-r1)*r3+dx05ddd*r4);
        dx05aaa = paaa; dx05ada = pada; dx05daa = pdaa; dx05dda = pdda; dx05aad = paad; dx05add = padd; dx05dad = pdad; dx05ddd = 1.00 - dx05aaa - dx05daa - dx05ada - dx05dda - dx05aad - dx05dad - dx05add;

if(m==0){
dx06aaa = eqaaa;
dx06daa = eqdaa;
dx06ada = eqada - epsilon;
dx06dda = eqdda;
dx06aad = eqaad;
dx06dad = eqdad;
dx06add = eqadd;
dx06ddd = 1.00 - dx06aaa - dx06daa - dx06ada - dx06dda - dx06aad - dx06dad - dx06add;
}
        avefit = dx06aaa*faaa + dx06ada*fada + dx06daa *fdaa + dx06dda*fdda + dx06aad*faad   + dx06add*fadd  + dx06dad*fdad  + dx06ddd *fddd;

        dx06aaa = dx06aaa*faaa/avefit;
        dx06ada = dx06ada*fada/avefit;
        dx06daa = dx06daa*fdaa/avefit;
        dx06dda = dx06dda*fdda/avefit;
        dx06aad = dx06aad*faad/avefit;
        dx06add = dx06add*fadd/avefit;
        dx06dad = dx06dad*fdad/avefit;
        dx06ddd = dx06ddd*fddd/avefit;

        paaa = dx06aaa*(1.0-dx06dad-dx06ddd-r2*dx06add-r1*dx06dda+RF1*dx06dad+(1.0-r1)*(1.0-r3)*dx06ddd)+dx06ada*(r2*dx06aad+r1*dx06daa+r1*r3*dx06dad)+dx06aad*((1.0-RF1)*dx06daa+(1.0-r1)*r3*dx06dda)+r1*(1.0-r3)*dx06add*dx06daa;
        pada = dx06ada*(1.0-dx06dad-dx06ddd-r2*dx06aad-r1*dx06daa+(1.0-r1)*(1.0-r3)*dx06dad+RF1*dx06ddd)+ dx06aaa*(r2*dx06add+r1*dx06dda+r1*r3*dx06ddd)+r1*(1.0-r3)*dx06aad*dx06dda+dx06add*((1.0-r1)*r3*dx06daa+(1.0-RF1)*dx06dda);
        paad = dx06aad*(1.0-r1*dx06ddd-dx06daa+RF2*dx06dda+RF1*dx06daa-r2*dx06ada)+dx06aaa*(r2*dx06add+(1.0-RF1)*dx06dad+(1.0-r1)*r3*dx06ddd)+dx06ada*dx06dad*r1*(1.0-r3)+dx06add*(dx06daa*r1*r3+r1*dx06dad);
        padd = dx06add*(1.0-r2*dx06aaa+RF2*dx06daa+(RF1-1.0)*dx06dda-r1*dx06dad)+dx06aaa*dx06ddd*r1*(1.0-r3)+dx06ada*(dx06aad*r2+dx06dad*(1.0-r1)*r3+dx06ddd*(1.0-RF1))+dx06aad*(dx06dda*r1*r3+dx06ddd*r1);
        pdaa = dx06daa*(1.0-r1*dx06ada+RF2*dx06add-r4*dx06ddd+(RF1-1.0)*dx06aad)+dx06aaa*(r1*dx06dda+(1.0-RF1)*dx06dad+r1*(1.0-r3)*dx06ddd) +dx06ada*dx06dad*(1.0-r1)*r3+dx06aad*dx06dda*r1*r3+dx06dda*dx06dad*r4;
        pdda = dx06dda*(1.0-r1*dx06aaa+RF2*dx06aad+(RF1-1.0)*dx06add-r4*dx06dad)+(1.0-r1)*r3*dx06aaa*dx06ddd+dx06ada*(dx06daa*r1+dx06dad*r1*(1.0-r3)+dx06ddd*(1.0-RF1))+dx06daa*(dx06add*r1*r3+dx06ddd*r4);
        pdad = dx06dad*(1.0+(RF1-1.0)*dx06aaa+RF2*dx06ada-r1*dx06add-r4*dx06dda)+dx06aaa*dx06ddd*r1*r3+dx06aad*(dx06daa*(1.0-RF1)+dx06dda*r1*(1.0-r3)+dx06ddd*r1)+dx06daa*(dx06add*(1.0-r1)*r3+dx06ddd*r4);
        dx06aaa = paaa; dx06ada = pada; dx06daa = pdaa; dx06dda = pdda; dx06aad = paad; dx06add = padd; dx06dad = pdad; dx06ddd = 1.00 - dx06aaa - dx06daa - dx06ada - dx06dda - dx06aad - dx06dad - dx06add;

//////////fourth column
if(m==0){
dx07aaa = eqaaa;
dx07daa = eqdaa;
dx07ada = eqada;
dx07dda = eqdda + epsilon;
dx07aad = eqaad;
dx07dad = eqdad;
dx07add = eqadd;
dx07ddd = 1.00 - dx07aaa - dx07daa - dx07ada - dx07dda - dx07aad - dx07dad - dx07add;
}

        avefit = dx07aaa*faaa + dx07ada*fada + dx07daa *fdaa + dx07dda*fdda + dx07aad*faad   + dx07add*fadd  + dx07dad*fdad  + dx07ddd *fddd;

        dx07aaa = dx07aaa*faaa/avefit;
        dx07ada = dx07ada*fada/avefit;
        dx07daa = dx07daa*fdaa/avefit;
        dx07dda = dx07dda*fdda/avefit;
        dx07aad = dx07aad*faad/avefit;
        dx07add = dx07add*fadd/avefit;
        dx07dad = dx07dad*fdad/avefit;
        dx07ddd = dx07ddd*fddd/avefit;

        paaa = dx07aaa*(1.0-dx07dad-dx07ddd-r2*dx07add-r1*dx07dda+RF1*dx07dad+(1.0-r1)*(1.0-r3)*dx07ddd)+dx07ada*(r2*dx07aad+r1*dx07daa+r1*r3*dx07dad)+dx07aad*((1.0-RF1)*dx07daa+(1.0-r1)*r3*dx07dda)+r1*(1.0-r3)*dx07add*dx07daa;
        pada = dx07ada*(1.0-dx07dad-dx07ddd-r2*dx07aad-r1*dx07daa+(1.0-r1)*(1.0-r3)*dx07dad+RF1*dx07ddd)+ dx07aaa*(r2*dx07add+r1*dx07dda+r1*r3*dx07ddd)+r1*(1.0-r3)*dx07aad*dx07dda+dx07add*((1.0-r1)*r3*dx07daa+(1.0-RF1)*dx07dda);
        paad = dx07aad*(1.0-r1*dx07ddd-dx07daa+RF2*dx07dda+RF1*dx07daa-r2*dx07ada)+dx07aaa*(r2*dx07add+(1.0-RF1)*dx07dad+(1.0-r1)*r3*dx07ddd)+dx07ada*dx07dad*r1*(1.0-r3)+dx07add*(dx07daa*r1*r3+r1*dx07dad);
        padd = dx07add*(1.0-r2*dx07aaa+RF2*dx07daa+(RF1-1.0)*dx07dda-r1*dx07dad)+dx07aaa*dx07ddd*r1*(1.0-r3)+dx07ada*(dx07aad*r2+dx07dad*(1.0-r1)*r3+dx07ddd*(1.0-RF1))+dx07aad*(dx07dda*r1*r3+dx07ddd*r1);
        pdaa = dx07daa*(1.0-r1*dx07ada+RF2*dx07add-r4*dx07ddd+(RF1-1.0)*dx07aad)+dx07aaa*(r1*dx07dda+(1.0-RF1)*dx07dad+r1*(1.0-r3)*dx07ddd) +dx07ada*dx07dad*(1.0-r1)*r3+dx07aad*dx07dda*r1*r3+dx07dda*dx07dad*r4;
        pdda = dx07dda*(1.0-r1*dx07aaa+RF2*dx07aad+(RF1-1.0)*dx07add-r4*dx07dad)+(1.0-r1)*r3*dx07aaa*dx07ddd+dx07ada*(dx07daa*r1+dx07dad*r1*(1.0-r3)+dx07ddd*(1.0-RF1))+dx07daa*(dx07add*r1*r3+dx07ddd*r4);
        pdad = dx07dad*(1.0+(RF1-1.0)*dx07aaa+RF2*dx07ada-r1*dx07add-r4*dx07dda)+dx07aaa*dx07ddd*r1*r3+dx07aad*(dx07daa*(1.0-RF1)+dx07dda*r1*(1.0-r3)+dx07ddd*r1)+dx07daa*(dx07add*(1.0-r1)*r3+dx07ddd*r4);
        dx07aaa = paaa; dx07ada = pada; dx07daa = pdaa; dx07dda = pdda; dx07aad = paad; dx07add = padd; dx07dad = pdad; dx07ddd = 1.00 - dx07aaa - dx07daa - dx07ada - dx07dda - dx07aad - dx07dad - dx07add;

if(m==0){
dx08aaa = eqaaa;
dx08daa = eqdaa;
dx08ada = eqada;
dx08dda = eqdda - epsilon;
dx08aad = eqaad;
dx08dad = eqdad;
dx08add = eqadd;
dx08ddd = 1.00 - dx08aaa - dx08daa - dx08ada - dx08dda - dx08aad - dx08dad - dx08add;
}


        avefit = dx08aaa*faaa + dx08ada*fada + dx08daa *fdaa + dx08dda*fdda + dx08aad*faad   + dx08add*fadd  + dx08dad*fdad  + dx08ddd *fddd;

        dx08aaa = dx08aaa*faaa/avefit;
        dx08ada = dx08ada*fada/avefit;
        dx08daa = dx08daa*fdaa/avefit;
        dx08dda = dx08dda*fdda/avefit;
        dx08aad = dx08aad*faad/avefit;
        dx08add = dx08add*fadd/avefit;
        dx08dad = dx08dad*fdad/avefit;
        dx08ddd = dx08ddd*fddd/avefit;

        paaa = dx08aaa*(1.0-dx08dad-dx08ddd-r2*dx08add-r1*dx08dda+RF1*dx08dad+(1.0-r1)*(1.0-r3)*dx08ddd)+dx08ada*(r2*dx08aad+r1*dx08daa+r1*r3*dx08dad)+dx08aad*((1.0-RF1)*dx08daa+(1.0-r1)*r3*dx08dda)+r1*(1.0-r3)*dx08add*dx08daa;
        pada = dx08ada*(1.0-dx08dad-dx08ddd-r2*dx08aad-r1*dx08daa+(1.0-r1)*(1.0-r3)*dx08dad+RF1*dx08ddd)+ dx08aaa*(r2*dx08add+r1*dx08dda+r1*r3*dx08ddd)+r1*(1.0-r3)*dx08aad*dx08dda+dx08add*((1.0-r1)*r3*dx08daa+(1.0-RF1)*dx08dda);
        paad = dx08aad*(1.0-r1*dx08ddd-dx08daa+RF2*dx08dda+RF1*dx08daa-r2*dx08ada)+dx08aaa*(r2*dx08add+(1.0-RF1)*dx08dad+(1.0-r1)*r3*dx08ddd)+dx08ada*dx08dad*r1*(1.0-r3)+dx08add*(dx08daa*r1*r3+r1*dx08dad);
        padd = dx08add*(1.0-r2*dx08aaa+RF2*dx08daa+(RF1-1.0)*dx08dda-r1*dx08dad)+dx08aaa*dx08ddd*r1*(1.0-r3)+dx08ada*(dx08aad*r2+dx08dad*(1.0-r1)*r3+dx08ddd*(1.0-RF1))+dx08aad*(dx08dda*r1*r3+dx08ddd*r1);
        pdaa = dx08daa*(1.0-r1*dx08ada+RF2*dx08add-r4*dx08ddd+(RF1-1.0)*dx08aad)+dx08aaa*(r1*dx08dda+(1.0-RF1)*dx08dad+r1*(1.0-r3)*dx08ddd) +dx08ada*dx08dad*(1.0-r1)*r3+dx08aad*dx08dda*r1*r3+dx08dda*dx08dad*r4;
        pdda = dx08dda*(1.0-r1*dx08aaa+RF2*dx08aad+(RF1-1.0)*dx08add-r4*dx08dad)+(1.0-r1)*r3*dx08aaa*dx08ddd+dx08ada*(dx08daa*r1+dx08dad*r1*(1.0-r3)+dx08ddd*(1.0-RF1))+dx08daa*(dx08add*r1*r3+dx08ddd*r4);
        pdad = dx08dad*(1.0+(RF1-1.0)*dx08aaa+RF2*dx08ada-r1*dx08add-r4*dx08dda)+dx08aaa*dx08ddd*r1*r3+dx08aad*(dx08daa*(1.0-RF1)+dx08dda*r1*(1.0-r3)+dx08ddd*r1)+dx08daa*(dx08add*(1.0-r1)*r3+dx08ddd*r4);
        dx08aaa = paaa; dx08ada = pada; dx08daa = pdaa; dx08dda = pdda; dx08aad = paad; dx08add = padd; dx08dad = pdad; dx08ddd = 1.00 - dx08aaa - dx08daa - dx08ada - dx08dda - dx08aad - dx08dad - dx08add;

///////// fifth column

if(m==0){
dx09aaa = eqaaa;
dx09daa = eqdaa;
dx09ada = eqada;
dx09dda = eqdda;
dx09aad = eqaad + epsilon;
dx09dad = eqdad;
dx09add = eqadd;
dx09ddd = 1.00 - dx09aaa - dx09daa - dx09ada - dx09dda - dx09aad - dx09dad - dx09add;
}

        avefit = dx09aaa*faaa + dx09ada*fada + dx09daa *fdaa + dx09dda*fdda + dx09aad*faad   + dx09add*fadd  + dx09dad*fdad  + dx09ddd *fddd;

        dx09aaa = dx09aaa*faaa/avefit;
        dx09ada = dx09ada*fada/avefit;
        dx09daa = dx09daa*fdaa/avefit;
        dx09dda = dx09dda*fdda/avefit;
        dx09aad = dx09aad*faad/avefit;
        dx09add = dx09add*fadd/avefit;
        dx09dad = dx09dad*fdad/avefit;
        dx09ddd = dx09ddd*fddd/avefit;

        paaa = dx09aaa*(1.0-dx09dad-dx09ddd-r2*dx09add-r1*dx09dda+RF1*dx09dad+(1.0-r1)*(1.0-r3)*dx09ddd)+dx09ada*(r2*dx09aad+r1*dx09daa+r1*r3*dx09dad)+dx09aad*((1.0-RF1)*dx09daa+(1.0-r1)*r3*dx09dda)+r1*(1.0-r3)*dx09add*dx09daa;
        pada = dx09ada*(1.0-dx09dad-dx09ddd-r2*dx09aad-r1*dx09daa+(1.0-r1)*(1.0-r3)*dx09dad+RF1*dx09ddd)+ dx09aaa*(r2*dx09add+r1*dx09dda+r1*r3*dx09ddd)+r1*(1.0-r3)*dx09aad*dx09dda+dx09add*((1.0-r1)*r3*dx09daa+(1.0-RF1)*dx09dda);
        paad = dx09aad*(1.0-r1*dx09ddd-dx09daa+RF2*dx09dda+RF1*dx09daa-r2*dx09ada)+dx09aaa*(r2*dx09add+(1.0-RF1)*dx09dad+(1.0-r1)*r3*dx09ddd)+dx09ada*dx09dad*r1*(1.0-r3)+dx09add*(dx09daa*r1*r3+r1*dx09dad);
        padd = dx09add*(1.0-r2*dx09aaa+RF2*dx09daa+(RF1-1.0)*dx09dda-r1*dx09dad)+dx09aaa*dx09ddd*r1*(1.0-r3)+dx09ada*(dx09aad*r2+dx09dad*(1.0-r1)*r3+dx09ddd*(1.0-RF1))+dx09aad*(dx09dda*r1*r3+dx09ddd*r1);
        pdaa = dx09daa*(1.0-r1*dx09ada+RF2*dx09add-r4*dx09ddd+(RF1-1.0)*dx09aad)+dx09aaa*(r1*dx09dda+(1.0-RF1)*dx09dad+r1*(1.0-r3)*dx09ddd) +dx09ada*dx09dad*(1.0-r1)*r3+dx09aad*dx09dda*r1*r3+dx09dda*dx09dad*r4;
        pdda = dx09dda*(1.0-r1*dx09aaa+RF2*dx09aad+(RF1-1.0)*dx09add-r4*dx09dad)+(1.0-r1)*r3*dx09aaa*dx09ddd+dx09ada*(dx09daa*r1+dx09dad*r1*(1.0-r3)+dx09ddd*(1.0-RF1))+dx09daa*(dx09add*r1*r3+dx09ddd*r4);
        pdad = dx09dad*(1.0+(RF1-1.0)*dx09aaa+RF2*dx09ada-r1*dx09add-r4*dx09dda)+dx09aaa*dx09ddd*r1*r3+dx09aad*(dx09daa*(1.0-RF1)+dx09dda*r1*(1.0-r3)+dx09ddd*r1)+dx09daa*(dx09add*(1.0-r1)*r3+dx09ddd*r4);
        dx09aaa = paaa; dx09ada = pada; dx09daa = pdaa; dx09dda = pdda; dx09aad = paad; dx09add = padd; dx09dad = pdad; dx09ddd = 1.00 - dx09aaa - dx09daa - dx09ada - dx09dda - dx09aad - dx09dad - dx09add;

if(m==0){
dx10aaa = eqaaa;
dx10daa = eqdaa;
dx10ada = eqada;
dx10dda = eqdda;
dx10aad = eqaad - epsilon;
dx10dad = eqdad;
dx10add = eqadd;
dx10ddd = 1.00 - dx10aaa - dx10daa - dx10ada - dx10dda - dx10aad - dx10dad - dx10add;
}
        avefit = dx10aaa*faaa + dx10ada*fada + dx10daa *fdaa + dx10dda*fdda + dx10aad*faad   + dx10add*fadd  + dx10dad*fdad  + dx10ddd *fddd;

        dx10aaa = dx10aaa*faaa/avefit;
        dx10ada = dx10ada*fada/avefit;
        dx10daa = dx10daa*fdaa/avefit;
        dx10dda = dx10dda*fdda/avefit;
        dx10aad = dx10aad*faad/avefit;
        dx10add = dx10add*fadd/avefit;
        dx10dad = dx10dad*fdad/avefit;
        dx10ddd = dx10ddd*fddd/avefit;

        paaa = dx10aaa*(1.0-dx10dad-dx10ddd-r2*dx10add-r1*dx10dda+RF1*dx10dad+(1.0-r1)*(1.0-r3)*dx10ddd)+dx10ada*(r2*dx10aad+r1*dx10daa+r1*r3*dx10dad)+dx10aad*((1.0-RF1)*dx10daa+(1.0-r1)*r3*dx10dda)+r1*(1.0-r3)*dx10add*dx10daa;
        pada = dx10ada*(1.0-dx10dad-dx10ddd-r2*dx10aad-r1*dx10daa+(1.0-r1)*(1.0-r3)*dx10dad+RF1*dx10ddd)+ dx10aaa*(r2*dx10add+r1*dx10dda+r1*r3*dx10ddd)+r1*(1.0-r3)*dx10aad*dx10dda+dx10add*((1.0-r1)*r3*dx10daa+(1.0-RF1)*dx10dda);
        paad = dx10aad*(1.0-r1*dx10ddd-dx10daa+RF2*dx10dda+RF1*dx10daa-r2*dx10ada)+dx10aaa*(r2*dx10add+(1.0-RF1)*dx10dad+(1.0-r1)*r3*dx10ddd)+dx10ada*dx10dad*r1*(1.0-r3)+dx10add*(dx10daa*r1*r3+r1*dx10dad);
        padd = dx10add*(1.0-r2*dx10aaa+RF2*dx10daa+(RF1-1.0)*dx10dda-r1*dx10dad)+dx10aaa*dx10ddd*r1*(1.0-r3)+dx10ada*(dx10aad*r2+dx10dad*(1.0-r1)*r3+dx10ddd*(1.0-RF1))+dx10aad*(dx10dda*r1*r3+dx10ddd*r1);
        pdaa = dx10daa*(1.0-r1*dx10ada+RF2*dx10add-r4*dx10ddd+(RF1-1.0)*dx10aad)+dx10aaa*(r1*dx10dda+(1.0-RF1)*dx10dad+r1*(1.0-r3)*dx10ddd) +dx10ada*dx10dad*(1.0-r1)*r3+dx10aad*dx10dda*r1*r3+dx10dda*dx10dad*r4;
        pdda = dx10dda*(1.0-r1*dx10aaa+RF2*dx10aad+(RF1-1.0)*dx10add-r4*dx10dad)+(1.0-r1)*r3*dx10aaa*dx10ddd+dx10ada*(dx10daa*r1+dx10dad*r1*(1.0-r3)+dx10ddd*(1.0-RF1))+dx10daa*(dx10add*r1*r3+dx10ddd*r4);
        pdad = dx10dad*(1.0+(RF1-1.0)*dx10aaa+RF2*dx10ada-r1*dx10add-r4*dx10dda)+dx10aaa*dx10ddd*r1*r3+dx10aad*(dx10daa*(1.0-RF1)+dx10dda*r1*(1.0-r3)+dx10ddd*r1)+dx10daa*(dx10add*(1.0-r1)*r3+dx10ddd*r4);
        dx10aaa = paaa; dx10ada = pada; dx10daa = pdaa; dx10dda = pdda; dx10aad = paad; dx10add = padd; dx10dad = pdad; dx10ddd = 1.00 - dx10aaa - dx10daa - dx10ada - dx10dda - dx10aad - dx10dad - dx10add;

////////////sixth column

if(m==0){
dx11aaa = eqaaa;
dx11daa = eqdaa;
dx11ada = eqada;
dx11dda = eqdda;
dx11aad = eqaad;
dx11dad = eqdad + epsilon;
dx11add = eqadd;
dx11ddd = 1.00 - dx11aaa - dx11daa - dx11ada - dx11dda - dx11aad - dx11dad - dx11add;
}
        avefit = dx11aaa*faaa + dx11ada*fada + dx11daa *fdaa + dx11dda*fdda + dx11aad*faad   + dx11add*fadd  + dx11dad*fdad  + dx11ddd *fddd;

        dx11aaa = dx11aaa*faaa/avefit;
        dx11ada = dx11ada*fada/avefit;
        dx11daa = dx11daa*fdaa/avefit;
        dx11dda = dx11dda*fdda/avefit;
        dx11aad = dx11aad*faad/avefit;
        dx11add = dx11add*fadd/avefit;
        dx11dad = dx11dad*fdad/avefit;
        dx11ddd = dx11ddd*fddd/avefit;

        paaa = dx11aaa*(1.0-dx11dad-dx11ddd-r2*dx11add-r1*dx11dda+RF1*dx11dad+(1.0-r1)*(1.0-r3)*dx11ddd)+dx11ada*(r2*dx11aad+r1*dx11daa+r1*r3*dx11dad)+dx11aad*((1.0-RF1)*dx11daa+(1.0-r1)*r3*dx11dda)+r1*(1.0-r3)*dx11add*dx11daa;
        pada = dx11ada*(1.0-dx11dad-dx11ddd-r2*dx11aad-r1*dx11daa+(1.0-r1)*(1.0-r3)*dx11dad+RF1*dx11ddd)+ dx11aaa*(r2*dx11add+r1*dx11dda+r1*r3*dx11ddd)+r1*(1.0-r3)*dx11aad*dx11dda+dx11add*((1.0-r1)*r3*dx11daa+(1.0-RF1)*dx11dda);
        paad = dx11aad*(1.0-r1*dx11ddd-dx11daa+RF2*dx11dda+RF1*dx11daa-r2*dx11ada)+dx11aaa*(r2*dx11add+(1.0-RF1)*dx11dad+(1.0-r1)*r3*dx11ddd)+dx11ada*dx11dad*r1*(1.0-r3)+dx11add*(dx11daa*r1*r3+r1*dx11dad);
        padd = dx11add*(1.0-r2*dx11aaa+RF2*dx11daa+(RF1-1.0)*dx11dda-r1*dx11dad)+dx11aaa*dx11ddd*r1*(1.0-r3)+dx11ada*(dx11aad*r2+dx11dad*(1.0-r1)*r3+dx11ddd*(1.0-RF1))+dx11aad*(dx11dda*r1*r3+dx11ddd*r1);
        pdaa = dx11daa*(1.0-r1*dx11ada+RF2*dx11add-r4*dx11ddd+(RF1-1.0)*dx11aad)+dx11aaa*(r1*dx11dda+(1.0-RF1)*dx11dad+r1*(1.0-r3)*dx11ddd) +dx11ada*dx11dad*(1.0-r1)*r3+dx11aad*dx11dda*r1*r3+dx11dda*dx11dad*r4;
        pdda = dx11dda*(1.0-r1*dx11aaa+RF2*dx11aad+(RF1-1.0)*dx11add-r4*dx11dad)+(1.0-r1)*r3*dx11aaa*dx11ddd+dx11ada*(dx11daa*r1+dx11dad*r1*(1.0-r3)+dx11ddd*(1.0-RF1))+dx11daa*(dx11add*r1*r3+dx11ddd*r4);
        pdad = dx11dad*(1.0+(RF1-1.0)*dx11aaa+RF2*dx11ada-r1*dx11add-r4*dx11dda)+dx11aaa*dx11ddd*r1*r3+dx11aad*(dx11daa*(1.0-RF1)+dx11dda*r1*(1.0-r3)+dx11ddd*r1)+dx11daa*(dx11add*(1.0-r1)*r3+dx11ddd*r4);
        dx11aaa = paaa; dx11ada = pada; dx11daa = pdaa; dx11dda = pdda; dx11aad = paad; dx11add = padd; dx11dad = pdad; dx11ddd = 1.00 - dx11aaa - dx11daa - dx11ada - dx11dda - dx11aad - dx11dad - dx11add;

if(m==0){
dx12aaa = eqaaa;
dx12daa = eqdaa;
dx12ada = eqada;
dx12dda = eqdda;
dx12aad = eqaad;
dx12dad = eqdad - epsilon;
dx12add = eqadd;
dx12ddd = 1.00 - dx12aaa - dx12daa - dx12ada - dx12dda - dx12aad - dx12dad - dx12add;
}
        avefit = dx12aaa*faaa + dx12ada*fada + dx12daa *fdaa + dx12dda*fdda + dx12aad*faad   + dx12add*fadd  + dx12dad*fdad  + dx12ddd *fddd;

        dx12aaa = dx12aaa*faaa/avefit;
        dx12ada = dx12ada*fada/avefit;
        dx12daa = dx12daa*fdaa/avefit;
        dx12dda = dx12dda*fdda/avefit;
        dx12aad = dx12aad*faad/avefit;
        dx12add = dx12add*fadd/avefit;
        dx12dad = dx12dad*fdad/avefit;
        dx12ddd = dx12ddd*fddd/avefit;

        paaa = dx12aaa*(1.0-dx12dad-dx12ddd-r2*dx12add-r1*dx12dda+RF1*dx12dad+(1.0-r1)*(1.0-r3)*dx12ddd)+dx12ada*(r2*dx12aad+r1*dx12daa+r1*r3*dx12dad)+dx12aad*((1.0-RF1)*dx12daa+(1.0-r1)*r3*dx12dda)+r1*(1.0-r3)*dx12add*dx12daa;
        pada = dx12ada*(1.0-dx12dad-dx12ddd-r2*dx12aad-r1*dx12daa+(1.0-r1)*(1.0-r3)*dx12dad+RF1*dx12ddd)+ dx12aaa*(r2*dx12add+r1*dx12dda+r1*r3*dx12ddd)+r1*(1.0-r3)*dx12aad*dx12dda+dx12add*((1.0-r1)*r3*dx12daa+(1.0-RF1)*dx12dda);
        paad = dx12aad*(1.0-r1*dx12ddd-dx12daa+RF2*dx12dda+RF1*dx12daa-r2*dx12ada)+dx12aaa*(r2*dx12add+(1.0-RF1)*dx12dad+(1.0-r1)*r3*dx12ddd)+dx12ada*dx12dad*r1*(1.0-r3)+dx12add*(dx12daa*r1*r3+r1*dx12dad);
        padd = dx12add*(1.0-r2*dx12aaa+RF2*dx12daa+(RF1-1.0)*dx12dda-r1*dx12dad)+dx12aaa*dx12ddd*r1*(1.0-r3)+dx12ada*(dx12aad*r2+dx12dad*(1.0-r1)*r3+dx12ddd*(1.0-RF1))+dx12aad*(dx12dda*r1*r3+dx12ddd*r1);
        pdaa = dx12daa*(1.0-r1*dx12ada+RF2*dx12add-r4*dx12ddd+(RF1-1.0)*dx12aad)+dx12aaa*(r1*dx12dda+(1.0-RF1)*dx12dad+r1*(1.0-r3)*dx12ddd) +dx12ada*dx12dad*(1.0-r1)*r3+dx12aad*dx12dda*r1*r3+dx12dda*dx12dad*r4;
        pdda = dx12dda*(1.0-r1*dx12aaa+RF2*dx12aad+(RF1-1.0)*dx12add-r4*dx12dad)+(1.0-r1)*r3*dx12aaa*dx12ddd+dx12ada*(dx12daa*r1+dx12dad*r1*(1.0-r3)+dx12ddd*(1.0-RF1))+dx12daa*(dx12add*r1*r3+dx12ddd*r4);
        pdad = dx12dad*(1.0+(RF1-1.0)*dx12aaa+RF2*dx12ada-r1*dx12add-r4*dx12dda)+dx12aaa*dx12ddd*r1*r3+dx12aad*(dx12daa*(1.0-RF1)+dx12dda*r1*(1.0-r3)+dx12ddd*r1)+dx12daa*(dx12add*(1.0-r1)*r3+dx12ddd*r4);
        dx12aaa = paaa; dx12ada = pada; dx12daa = pdaa; dx12dda = pdda; dx12aad = paad; dx12add = padd; dx12dad = pdad; dx12ddd = 1.00 - dx12aaa - dx12daa - dx12ada - dx12dda - dx12aad - dx12dad - dx12add;

////////////seventh column

if(m==0){
dx13aaa = eqaaa;
dx13daa = eqdaa;
dx13ada = eqada;
dx13dda = eqdda;
dx13aad = eqaad;
dx13dad = eqdad;
dx13add = eqadd + epsilon;
dx13ddd = 1.00 - dx13aaa - dx13daa - dx13ada - dx13dda - dx13aad - dx13dad - dx13add;
}

        avefit = dx13aaa*faaa + dx13ada*fada + dx13daa *fdaa + dx13dda*fdda + dx13aad*faad   + dx13add*fadd  + dx13dad*fdad  + dx13ddd *fddd;

        dx13aaa = dx13aaa*faaa/avefit;
        dx13ada = dx13ada*fada/avefit;
        dx13daa = dx13daa*fdaa/avefit;
        dx13dda = dx13dda*fdda/avefit;
        dx13aad = dx13aad*faad/avefit;
        dx13add = dx13add*fadd/avefit;
        dx13dad = dx13dad*fdad/avefit;
        dx13ddd = dx13ddd*fddd/avefit;

        paaa = dx13aaa*(1.0-dx13dad-dx13ddd-r2*dx13add-r1*dx13dda+RF1*dx13dad+(1.0-r1)*(1.0-r3)*dx13ddd)+dx13ada*(r2*dx13aad+r1*dx13daa+r1*r3*dx13dad)+dx13aad*((1.0-RF1)*dx13daa+(1.0-r1)*r3*dx13dda)+r1*(1.0-r3)*dx13add*dx13daa;
        pada = dx13ada*(1.0-dx13dad-dx13ddd-r2*dx13aad-r1*dx13daa+(1.0-r1)*(1.0-r3)*dx13dad+RF1*dx13ddd)+ dx13aaa*(r2*dx13add+r1*dx13dda+r1*r3*dx13ddd)+r1*(1.0-r3)*dx13aad*dx13dda+dx13add*((1.0-r1)*r3*dx13daa+(1.0-RF1)*dx13dda);
        paad = dx13aad*(1.0-r1*dx13ddd-dx13daa+RF2*dx13dda+RF1*dx13daa-r2*dx13ada)+dx13aaa*(r2*dx13add+(1.0-RF1)*dx13dad+(1.0-r1)*r3*dx13ddd)+dx13ada*dx13dad*r1*(1.0-r3)+dx13add*(dx13daa*r1*r3+r1*dx13dad);
        padd = dx13add*(1.0-r2*dx13aaa+RF2*dx13daa+(RF1-1.0)*dx13dda-r1*dx13dad)+dx13aaa*dx13ddd*r1*(1.0-r3)+dx13ada*(dx13aad*r2+dx13dad*(1.0-r1)*r3+dx13ddd*(1.0-RF1))+dx13aad*(dx13dda*r1*r3+dx13ddd*r1);
        pdaa = dx13daa*(1.0-r1*dx13ada+RF2*dx13add-r4*dx13ddd+(RF1-1.0)*dx13aad)+dx13aaa*(r1*dx13dda+(1.0-RF1)*dx13dad+r1*(1.0-r3)*dx13ddd) +dx13ada*dx13dad*(1.0-r1)*r3+dx13aad*dx13dda*r1*r3+dx13dda*dx13dad*r4;
        pdda = dx13dda*(1.0-r1*dx13aaa+RF2*dx13aad+(RF1-1.0)*dx13add-r4*dx13dad)+(1.0-r1)*r3*dx13aaa*dx13ddd+dx13ada*(dx13daa*r1+dx13dad*r1*(1.0-r3)+dx13ddd*(1.0-RF1))+dx13daa*(dx13add*r1*r3+dx13ddd*r4);
        pdad = dx13dad*(1.0+(RF1-1.0)*dx13aaa+RF2*dx13ada-r1*dx13add-r4*dx13dda)+dx13aaa*dx13ddd*r1*r3+dx13aad*(dx13daa*(1.0-RF1)+dx13dda*r1*(1.0-r3)+dx13ddd*r1)+dx13daa*(dx13add*(1.0-r1)*r3+dx13ddd*r4);
        dx13aaa = paaa; dx13ada = pada; dx13daa = pdaa; dx13dda = pdda; dx13aad = paad; dx13add = padd; dx13dad = pdad; dx13ddd = 1.00 - dx13aaa - dx13daa - dx13ada - dx13dda - dx13aad - dx13dad - dx13add;

if(m==0){
dx14aaa = eqaaa;
dx14daa = eqdaa;
dx14ada = eqada;
dx14dda = eqdda;
dx14aad = eqaad;
dx14dad = eqdad;
dx14add = eqadd - epsilon;
dx14ddd = 1.00 - dx14aaa - dx14daa - dx14ada - dx14dda - dx14aad - dx14dad - dx14add;
}
        avefit = dx14aaa*faaa + dx14ada*fada + dx14daa *fdaa + dx14dda*fdda + dx14aad*faad   + dx14add*fadd  + dx14dad*fdad  + dx14ddd *fddd;

        dx14aaa = dx14aaa*faaa/avefit;
        dx14ada = dx14ada*fada/avefit;
        dx14daa = dx14daa*fdaa/avefit;
        dx14dda = dx14dda*fdda/avefit;
        dx14aad = dx14aad*faad/avefit;
        dx14add = dx14add*fadd/avefit;
        dx14dad = dx14dad*fdad/avefit;
        dx14ddd = dx14ddd*fddd/avefit;

        paaa = dx14aaa*(1.0-dx14dad-dx14ddd-r2*dx14add-r1*dx14dda+RF1*dx14dad+(1.0-r1)*(1.0-r3)*dx14ddd)+dx14ada*(r2*dx14aad+r1*dx14daa+r1*r3*dx14dad)+dx14aad*((1.0-RF1)*dx14daa+(1.0-r1)*r3*dx14dda)+r1*(1.0-r3)*dx14add*dx14daa;
        pada = dx14ada*(1.0-dx14dad-dx14ddd-r2*dx14aad-r1*dx14daa+(1.0-r1)*(1.0-r3)*dx14dad+RF1*dx14ddd)+ dx14aaa*(r2*dx14add+r1*dx14dda+r1*r3*dx14ddd)+r1*(1.0-r3)*dx14aad*dx14dda+dx14add*((1.0-r1)*r3*dx14daa+(1.0-RF1)*dx14dda);
        paad = dx14aad*(1.0-r1*dx14ddd-dx14daa+RF2*dx14dda+RF1*dx14daa-r2*dx14ada)+dx14aaa*(r2*dx14add+(1.0-RF1)*dx14dad+(1.0-r1)*r3*dx14ddd)+dx14ada*dx14dad*r1*(1.0-r3)+dx14add*(dx14daa*r1*r3+r1*dx14dad);
        padd = dx14add*(1.0-r2*dx14aaa+RF2*dx14daa+(RF1-1.0)*dx14dda-r1*dx14dad)+dx14aaa*dx14ddd*r1*(1.0-r3)+dx14ada*(dx14aad*r2+dx14dad*(1.0-r1)*r3+dx14ddd*(1.0-RF1))+dx14aad*(dx14dda*r1*r3+dx14ddd*r1);
        pdaa = dx14daa*(1.0-r1*dx14ada+RF2*dx14add-r4*dx14ddd+(RF1-1.0)*dx14aad)+dx14aaa*(r1*dx14dda+(1.0-RF1)*dx14dad+r1*(1.0-r3)*dx14ddd) +dx14ada*dx14dad*(1.0-r1)*r3+dx14aad*dx14dda*r1*r3+dx14dda*dx14dad*r4;
        pdda = dx14dda*(1.0-r1*dx14aaa+RF2*dx14aad+(RF1-1.0)*dx14add-r4*dx14dad)+(1.0-r1)*r3*dx14aaa*dx14ddd+dx14ada*(dx14daa*r1+dx14dad*r1*(1.0-r3)+dx14ddd*(1.0-RF1))+dx14daa*(dx14add*r1*r3+dx14ddd*r4);
        pdad = dx14dad*(1.0+(RF1-1.0)*dx14aaa+RF2*dx14ada-r1*dx14add-r4*dx14dda)+dx14aaa*dx14ddd*r1*r3+dx14aad*(dx14daa*(1.0-RF1)+dx14dda*r1*(1.0-r3)+dx14ddd*r1)+dx14daa*(dx14add*(1.0-r1)*r3+dx14ddd*r4);
        dx14aaa = paaa; dx14ada = pada; dx14daa = pdaa; dx14dda = pdda; dx14aad = paad; dx14add = padd; dx14dad = pdad; dx14ddd = 1.00 - dx14aaa - dx14daa - dx14ada - dx14dda - dx14aad - dx14dad - dx14add;

/////////// now print the Jacobian

if(m==(PHASE - 1)){
e11 = (dx01aaa - dx02aaa)/(2.0* epsilon);
e21 = (dx01daa - dx02daa)/(2.0* epsilon);
e31 = (dx01ada - dx02ada)/(2.0* epsilon);
e41 = (dx01dda - dx02dda)/(2.0* epsilon);
e51 = (dx01aad - dx02aad)/(2.0* epsilon);
e61 = (dx01dad - dx02dad)/(2.0* epsilon);
e71 = (dx01add - dx02add)/(2.0* epsilon);

e12 = (dx03aaa - dx04aaa)/(2.0* epsilon);
e22 = (dx03daa - dx04daa)/(2.0* epsilon);
e32 = (dx03ada - dx04ada)/(2.0* epsilon);
e42 = (dx03dda - dx04dda)/(2.0* epsilon);
e52 = (dx03aad - dx04aad)/(2.0* epsilon);
e62 = (dx03dad - dx04dad)/(2.0* epsilon);
e72 = (dx03add - dx04add)/(2.0* epsilon);

e13 = (dx05aaa - dx06aaa)/(2.0* epsilon);
e23 = (dx05daa - dx06daa)/(2.0* epsilon);
e33 = (dx05ada - dx06ada)/(2.0* epsilon);
e43 = (dx05dda - dx06dda)/(2.0* epsilon);
e53 = (dx05aad - dx06aad)/(2.0* epsilon);
e63 = (dx05dad - dx06dad)/(2.0* epsilon);
e73 = (dx05add - dx06add)/(2.0* epsilon);

e14 = (dx07aaa - dx08aaa)/(2.0* epsilon);
e24 = (dx07daa - dx08daa)/(2.0* epsilon);
e34 = (dx07ada - dx08ada)/(2.0* epsilon);
e44 = (dx07dda - dx08dda)/(2.0* epsilon);
e54 = (dx07aad - dx08aad)/(2.0* epsilon);
e64 = (dx07dad - dx08dad)/(2.0* epsilon);
e74 = (dx07add - dx08add)/(2.0* epsilon);

e15 = (dx09aaa - dx10aaa)/(2.0* epsilon);
e25 = (dx09daa - dx10daa)/(2.0* epsilon);
e35 = (dx09ada - dx10ada)/(2.0* epsilon);
e45 = (dx09dda - dx10dda)/(2.0* epsilon);
e55 = (dx09aad - dx10aad)/(2.0* epsilon);
e65 = (dx09dad - dx10dad)/(2.0* epsilon);
e75 = (dx09add - dx10add)/(2.0* epsilon);

e16 = (dx11aaa - dx12aaa)/(2.0* epsilon);
e26 = (dx11daa - dx12daa)/(2.0* epsilon);
e36 = (dx11ada - dx12ada)/(2.0* epsilon);
e46 = (dx11dda - dx12dda)/(2.0* epsilon);
e56 = (dx11aad - dx12aad)/(2.0* epsilon);
e66 = (dx11dad - dx12dad)/(2.0* epsilon);
e76 = (dx11add - dx12add)/(2.0* epsilon);

e17 = (dx13aaa - dx14aaa)/(2.0* epsilon);
e27 = (dx13daa - dx14daa)/(2.0* epsilon);
e37 = (dx13ada - dx14ada)/(2.0* epsilon);
e47 = (dx13dda - dx14dda)/(2.0* epsilon);
e57 = (dx13aad - dx14aad)/(2.0* epsilon);
e67 = (dx13dad - dx14dad)/(2.0* epsilon);
e77 = (dx13add - dx14add)/(2.0* epsilon);

/*            // this output can be directly copied into Mathematica
fprintf(out2,"Jacob = {{%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf},{%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf},{%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf},{%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf},{%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf},{%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf},{%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf}}; \n", e11, e12, e13, e14, e15, e16, e17, e21, e22, e23, e24, e25, e26, e27, e31, e32, e33, e34, e35, e36, e37, e41, e42, e43, e44, e45, e46, e47, e51, e52, e53, e54, e55, e56, e57, e61, e62, e63, e64, e65, e66, e67, e71, e72, e73, e74, e75, e76, e77);

fprintf(out2, "Eigenvalues[Jacob]\nAbs[Eigenvalues[Jacob]]  \n ");
*/
             // this part of the output can be directly copied into R to obtain eigenvalues
fprintf(out2, "Jacob <- c(%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf,%.16Lf)\n", e11, e12, e13, e14, e15, e16, e17, e21, e22, e23, e24, e25, e26, e27, e31, e32, e33, e34, e35, e36, e37, e41, e42, e43, e44, e45, e46, e47, e51, e52, e53, e54, e55, e56, e57, e61, e62, e63, e64, e65, e66, e67, e71, e72, e73, e74, e75, e76, e77);
fprintf(out2, "mat<-matrix(Jacob, nrow=7, byrow=TRUE)\nprint(eigen(mat)$values[1:3])\n");

}// end print jacob

}}
jacobdone = 1;
}
if(jacobdone == 1) break;
}  // end time

if(jacobdone == 0)fprintf(out2, "\n print('P %d S %Lf      r-old %Lf r-invader %Lf      NOT DONe NOT DOne NOT DONe NOT DOne NOT DONe NOT DONe') \n", PHASE, swc, rec2, rec2new);

}}  // closing allele freq combo loops
WCONT:;
}}  // closing rec combos loops

fclose(inp);
fclose(out2);
fclose(par);
return 0;
}  // end main
	
