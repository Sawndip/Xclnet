// Taken in part from Numerical Recipes in C
#include "NumericalTools.h"

#ifndef PI
#define PI 3.141592654
#define CHANGED_PI
#endif /*PI*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define NITER 4

#define IM1 2147483563
#define IM2 2147483399
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)


#define NUMBER_INTERVALS (1000)


// Define shape of potential well for double-well invivo simulations
double U(double x){
	double Ux;
	double alphaD = 0.0001171457;
	//double alphaP = 0.0000137076;
	double Gp = 0.00993916;
	double Gd = 0.03888171;
	
	//double Gx = (0.5) + Gp + Gd;
	
	Ux = ((1 - alphaD) * ((0.25 * x * x) - (0.5 * pow(x,3)) + (0.25 * pow(x,4))) - Gp * (x - 0.5 * x * x) + Gd * (0.5 * x * x));
	
	return Ux;
}

// Return variable with double-well 1Hz, invivo distribution for synapses
double invivo_double_well_distribution(long *uni_seed){
	static int first_run_flag = 0;
	
	double interval_width = 1.0/(double)NUMBER_INTERVALS;
	
	static double x[NUMBER_INTERVALS+1]; //x-values
	static double P[NUMBER_INTERVALS+1]; //PDF
	static double C[NUMBER_INTERVALS+1]; //CDF
	
	double uni_rnd = ran2(uni_seed);
	
	if(first_run_flag == 0){ // generate pdf and cdf
		double alphaD = 0.0001171457;
		double alphaP = 0.0000137076;
		double sigma = 3.35;
		
		double sigma_rho_sq = (sigma*sigma * (alphaD + alphaP));
		
		double sum_P = 0; // for normalisation
		
		for(int i = 0; i < NUMBER_INTERVALS+1; i++){
			x[i] = i * interval_width;
			P[i] = exp( (-2.0 * U(x[i])) / sigma_rho_sq);
			
			sum_P += P[i];
			
			//printf("%f %f %g %g\n", x[i], U(x[i]), P[i], sum_P);
		}
		
		double normaliser = sum_P * interval_width;
		for(int i = 0; i < NUMBER_INTERVALS+1; i++){
			P[i] = P[i]/normaliser;
			
			if (i > 0){
				C[i] = C[i-1] + P[i]*interval_width;
			}
			else{
				C[i] = P[i]*interval_width;
			}
			
			//printf("%f %g %g\n", x[i], P[i], C[i]);
		}
		
		first_run_flag = 1;
	} //end first run generation of pdf and cdf
	
	// now use uni_rnd to lookup x value from C
	double new_rand;
	double higher_x;
	double higher_y;
	double lower_x;
	double lower_y;
	
	//printf("uni_rnd %f\n", uni_rnd);
	for (int i = 0; i < NUMBER_INTERVALS+1; i++){
		//find sign change in difference
		double diff = (uni_rnd - C[i]);
		if (diff < 0){
			higher_x = x[i];
			lower_x = x[i-1];
			higher_y = C[i];
			lower_y = C[i-1];
			break;
		}
	}
	//printf("%g %g %g %g\n", higher_x, lower_x, higher_y, lower_y);
	// Don't forget swapping x's and y's here (hence the reversed formulae)
	double slope = (higher_x - lower_x) / (higher_y - lower_y);
	new_rand = slope * (uni_rnd - lower_y) + lower_x;
	
	return new_rand;
}


// ran2_resettable returns a Uniform(0,1) distributed value
// L'Ecuyer combination method, with particularly long period
// Combination: two linear congruential generators and Bays-Durham shuffle
// Usage: set (reset != 0) in order to reset internal variables of generator
//      this will make generator behave as if this is the first call
//  set original_seed_value with the originally defined seed value (not location)
float ran2_resettable(long *idum, int reset, long original_seed_value)
{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    if(reset){
        idum2 = 123456789;
        iy = 0;
        *idum = original_seed_value;
        //printf("Resetting ran2_resettable()\n");
    }
    float temp;
    
    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) { // prepare shuffle table
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1; // compute idum=(IA1*idum)%IM1 without overflow
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2; // compute idum2=(IA2*idum2)%IM2 without overflow
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV; // in range 0..NTAB-1
    iy=iv[j]-idum2; // shuffle idum and combine with idum2 for output
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX; // do not return end-point value
    else return temp;
}


// Modified version for expdev() which is coupled to ran2_resettable()
float expdev_resettable(long *idum, int reset, long original_seed_value)
{
    float ran2_resettable(long *idum, int reset, long original_seed_value);
    float dum;
    
    do{
        dum=ran2_resettable(idum, reset, original_seed_value);
        reset = 0;
    }while (dum == 0.0);
    return -log(dum);
}



//// Below here adapted from Numerical Recipes in C

// ran0 returns a Uniform(0,1) value
// Minimal Park and Miller method: multiplicative congruential algorithm (I_{j+1} = a I_j \mod m )
float ran0(long *idum)
{
	long k;
	float ans;

	*idum ^= MASK; // Masking XOR allows use of 0 in seed
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k; // compute idum=(IA*idum)%IM without overflow
	if (*idum < 0) *idum += IM;
	ans=AM*(*idum); // Convert to float
	*idum ^= MASK; // Unmasking XOR
	return ans;
}


// ran1 returns a Uniform(0,1) value
// Minimal Park, Miller alg with Bays-Durham shuffle and extras
float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) { // prepare shuffle table
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k; // compute idum=(IA*idum)%IM without overflow
	if (*idum < 0) *idum += IM;
	j=iy/NDIV; // in range 0..NTAB-1
	iy=iv[j]; // output previously stored value and refill shuffle table
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX; // do not return end-point value
	else return temp;
}


// ran2 returns a Uniform(0,1) distributed value
// L'Ecuyer combination method, with particularly long period
// Combination: two linear congruential generators and Bays-Durham shuffle
float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) { // prepare shuffle table
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1; // compute idum=(IA1*idum)%IM1 without overflow
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2; // compute idum2=(IA2*idum2)%IM2 without overflow
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV; // in range 0..NTAB-1
	iy=iv[j]-idum2; // shuffle idum and combine with idum2 for output
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX; // do not return end-point value
	else return temp;
}


// ran3 returns Uniform(0,1) value
// Knuth suggested subtractive method
// (what is: ? 55 seed addition Lagged Fibonacci generator ?)
float ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) { // initialise
		iff=1;
		mj=labs(MSEED-labs(*idum));
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) { // initialise table with not so random numbers
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++) // randomise the table elements a little more
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;  // indices
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp]; // subtractively generate new random number
	if (mj < MZ) mj += MBIG; // make sure it's in range
	ma[inext]=mj; // store it
	return mj*FAC;
}


// ran4 returns Uniform(0,1) value
// random deviates from DES-like hashing
float ran4(long *idum)
{
	void psdes(unsigned long *lword, unsigned long *irword);
	unsigned long irword,itemp,lword;
	static long idums = 0;
#if defined(vax) || defined(_vax_) || defined(__vax__) || defined(VAX)
	static unsigned long jflone = 0x00004080;
	static unsigned long jflmsk = 0xffff007f;
#else
	static unsigned long jflone = 0x3f800000;
	static unsigned long jflmsk = 0x007fffff;
#endif

	if (*idum < 0) {
		idums = -(*idum);
		*idum=1;
	}
	irword=(*idum);
	lword=idums;
	psdes(&lword,&irword); // pseudo-DES encode words
	itemp=jflone | (jflmsk & irword); // Mask to a floating number between 1 and 2
	++(*idum);
	return (*(float *)&itemp)-1.0; // Move range to (0,1)
}


// pseudo-DES hash function
void psdes(unsigned long *lword, unsigned long *irword)
{
	unsigned long i,ia,ib,iswap,itmph=0,itmpl=0;
	static unsigned long c1[NITER]={
		0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
	static unsigned long c2[NITER]={
		0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};

	for (i=0;i<NITER;i++) {
		ia=(iswap=(*irword)) ^ c1[i];
		itmpl = ia & 0xffff;
		itmph = ia >> 16;
		ib=itmpl*itmpl+ ~(itmph*itmph);
		*irword=(*lword) ^ (((ia = (ib >> 16) |
			((ib & 0xffff) << 16)) ^ c2[i])+itmpl*itmph);
		*lword=iswap;
	}
}


// gasdev returns a Normal(0,1) distributed value
float gasdev(long *idum)
{
	float ran2(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if (*idum < 0) iset=0;
	if (iset == 0){ // we don't have a value saved from last time
		do{
			v1=2.0*ran2(idum)-1.0; // Pick 2 Uni(-1,1) values
			v2=2.0*ran2(idum)-1.0;
			rsq=v1*v1+v2*v2; // See if they're inside unit circle
		}while (rsq >= 1.0 || rsq == 0.0); // if not, do again
		//CONSIDER: is it faster to repeat this while loop, or to use polar coordinates to guarantee membership of unit circle?
		
		fac=sqrt(-2.0*log(rsq)/rsq); // Box-Muller transform
		gset=v1*fac; // save second value for next time
		
		iset=1;
		return v2*fac;
	}
	else{ // return the value we saved from last time
		iset=0; 
		return gset; 
	}
}


float expdev(long *idum)
{
	float ran2(long *idum);
	float dum;

	do
		dum=ran2(idum);
	while (dum == 0.0);
	return -log(dum);
}


float gamdev(int ia, long *idum)
{
	float ran2(long *idum);
	//void nrerror(char error_text[]);
	int j;
	float am,e,s,v1,v2,x,y;

	if (ia < 1) nrerror("Error in routine gamdev");
	if (ia < 6) {
		x=1.0;
		for (j=1;j<=ia;j++) x *= ran2(idum);
		x = -log(x);
	} else {
		do {
			do { 
				do { 
					v1=ran2(idum);
					v2=2.0*ran2(idum)-1.0;
				} while (v1*v1+v2*v2 > 1.0);
				y=v2/v1;
				am=ia-1;
				s=sqrt(2.0*am+1.0);
				x=s*y+am;
			} while (x <= 0.0);
			e=(1.0+y*y)*exp(am*log(x/am)-s*y);
		} while (ran2(idum) > e);
	}
	return x;
}


float poidev(float xm, long *idum)
{
	float gammln(float xx);
	float ran2(long *idum);
	static float sq,alxm,g,oldm=(-1.0);
	float em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran2(idum);
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
				y=tan(PI*ran2(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran2(idum) > t);
	}
	return em;
}


float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
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


float bnldev(float pp, int n, long *idum)
{
	float gammln(float xx);
	float ran2(long *idum);
	int j;
	static int nold=(-1);
	float am,em,g,angle,p,bnl,sq,t,y;
	static float pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (ran2(idum) < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran2(idum);
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran2(idum);
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
				-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran2(idum) > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	printf("big fat error in numerical code");
	exit(1);
}


#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

#undef IM1
#undef IM2
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

#undef NITER

#ifdef CHANGED_PI
#undef PI
#undef CHANGED_PI
#endif /*CHANGED_PI*/
