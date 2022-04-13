#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC




double gasdev(long *idum)
{
	static int iset=0;
	static float gset;
	float fac, r, v1, v2;
	//     	double ran3();

	if(iset == 0) {
		do {
			v1 = 2.0*ran3(idum) - 1.0;
			v2 = 2.0*ran3(idum) - 1.0;
			r = v1*v1+v2*v2;
		} while(r >= 1.0);

		fac = sqrt(-2.0*log(r)/r);
		gset = v1*fac;
		iset = 1;
		return v2*fac;

	}
	else {
		iset = 0;
		return gset;
	}
}



float expdev(long *idum)
/* Returns an exponentially distributed, positive, random deviate of unit mean,
using ran3(idum) as the source of uniform deviates. */

{
  //    float ran3(long *idum);
    float dum;
    
    do
        dum = ran3(idum);
    while(dum == 0.0);
    return -log(dum);
}


double gammln(double xx)
{
	double x, tmp, ser;
	static double cof[6] = {76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003e-2, -0.536382e-5};
	int j;
	
	x = xx-1.0;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.0;
	for(j=0; j<=5; j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}





#define PI 3.141592654

double poidev(float xm,long *idum)
{
	static double sq, alxm, g, oldm=(-1.0);
	double em, t, y;

	if(xm < 12.0) {
		if(xm != oldm) {
			oldm = xm;
			g = exp(-xm);
		}
		em = -1;
		t = 1.0;
		do {
			em += 1.0;
			t *= ran3(idum);
		} while(t>g);
	}
	else
	{
		if(xm != oldm) {
			oldm = xm;
			sq = pow(2.0*xm, 0.5);
			alxm = log(xm);
			g = xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y = tan(PI*ran3(idum));
				em = sq*y+xm;
			} while(em < 0.0);
		em = floor(em);
		t = 0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0) - g);
		} while(ran3(idum) > t);
	}
	return em;
}

double bnldev(float pp, int n, long *idum)
{
        int j;
        static int nold=(-1);
        float am,em,g,angle,p,bnl,sq,t,y;
        static float pold=(-1.0),pc,plog,pclog,en,oldg;

        p=(pp <= 0.5 ? pp : 1.0-pp);
        am=n*p;
        if (n < 25) {
                bnl=0.0;
                for (j=1;j<=n;j++)
                        if (ran3(idum) < p) ++bnl;
        } else if (am < 1.0) {
                g=exp(-am);
                t=1.0;
                for (j=0;j<=n;j++) {
                        t *= ran3(idum);
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
                                angle=PI*ran3(idum);
                                y=tan(angle);
                                em=sq*y+am;
                        } while (em < 0.0 || em >= (en+1.0));
                        em=floor(em);
                        t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
                                -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
                } while (ran3(idum) > t);
                bnl=em;
        }
        if (p != pp) bnl=n-bnl;
        return bnl;
}
#undef PI

