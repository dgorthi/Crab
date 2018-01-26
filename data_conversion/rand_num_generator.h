#include <math.h>

#define IA 16807
#define IM 2147483647 
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB) 
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float gasdev(long *idum);
float ran1(long *idum);

float ran1(long *idum){
/*“Minimal” random number generator of Park and Miller with Bays-Durham shuffle and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter idum between successive deviates in a sequence. RNMX should approximate the largest floating value that is less than 1.*/
	int j;
	long k;
	static long iy=0; 
	static long iv[NTAB]; 
	float temp;
	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1; 
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
		//Be sure to prevent idum = 0.
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

float gasdev(long *idum){
/*Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum) as the source of uniform deviates.*/
	float ran1(long *idum); 
	static int iset=0; 
	static float gset; 
	float fac,rsq,v1,v2;
	if (*idum < 0) iset=0; 
	if (iset == 0) {
		do { 
			v1=2.0*ran1(idum)-1.0; 
			v2=2.0*ran1(idum)-1.0; 
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
