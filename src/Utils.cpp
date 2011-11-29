// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#########################################################################################
#include "Utils.h"

double Utils::dMod(double *x, double *y)
{
/*	
	double xa, ya, z;
	xa = *x;
	ya = *y;
	if ((ya = *y) < 0.)
		ya = -ya;
	z = drem(xa = *x, ya);
	if (xa > 0)
	{
		if (z < 0)
			z += ya;
	}
	else if (z > 0)
		z -= ya;
	return z;
*/
	double quotient;
	if( (quotient = *x / *y) >= 0)
		quotient = floor(quotient);
	else
		quotient = -floor(-quotient);
	return(*x - (*y) * quotient );
}

double Utils::dSign(double *a, double *b)
{
	double x;
	x = (*a >= 0 ? *a : -*a);
	return (*b >= 0 ? x : -x);
}

double Utils::ran2(long int *idum)
{
	static long int idum2 = 123456789;
	static long int iv[32] =
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0 };
	static long int iy = 0;

	long int i__1;
	double ret_val, d__1;
	static long int j, k;

	// This is a random number generator. See the article on
	// p522 of COMPUTERS IN PHYSICS, SEP/OCT 1992, in which the
	// authors claim that ran2 is a 'perfect' random generator.
	// Long period(>2x10^18) random generator of L'Ecuyer with
	// Bays-Durhamc shuffle and added safeguards. Returns a uniform
	// random deviate between 0.0 and 1.0 (exclusive of the endpoint */
	// values). Call with idum a negative long int to initialize; */
	// therefore, do not alter idum between successive deviates in a */
	// sequence. RNMX should approximate the largest floating value */
	// that is less than 1.0

	if (*idum <= 0)
	{
		// Computing MAX
		i__1 = -(*idum);
		*idum = MAX(i__1, 1);
		// be sure to prevent idum=0
		idum2 = *idum;
		for (j = 40; j >= 1; --j)
		{
			// load the shuffle table (after 8 warm-up
			k = *idum / 53668;
			*idum = (*idum - k * 53668) * 40014 - k * 12211;
			if (*idum < 0)
			{
				*idum += 2147483563;
			}
			if (j <= 32)
			{
				iv[j - 1] = *idum;
			}
			/* L11: */
		}
		iy = iv[0];
	}

	k = *idum / 53668;
	// start here when no initializing
	*idum = (*idum - k * 53668) * 40014 - k * 12211;

	// compute idum=mod(IA1*idum,IM1) without overflows by Schrage's
	// method

	if (*idum < 0)
	{
		*idum += 2147483563;
	}
	k = idum2 / 52774;
	idum2 = (idum2 - k * 52774) * 40692 - k * 3791;

	// compute idum2=mod(IA2*idum2,IM2) likewise

	if (idum2 < 0)
	{
		idum2 += 2147483399;
	}
	j = iy / 67108862 + 1;
	// will be in the range 1:NTAB
	iy = iv[j - 1] - idum2;

	// here idum is shuffled, idum and idum2 are combined to grnerate
	/*        output */

	iv[j - 1] = *idum;
	if (iy < 1)
	{
		iy += 2147483562;
	}
	/* Computing MIN */
	d__1 = iy * 4.6566130573917691e-10;
	ret_val = MIN(d__1, .99999987999999995);

	/* because users don't expect endpoints */
	return ret_val;
}

double Utils::yyGas(long int *idum)
{
	static long int usey = 1;
	double ret_val;

	static double s, x, y;
	static double root, ranbyx, ranbyy;

	// This function YYGas use the polar method from George Marsaglia to
	// generate a pair of independent standard normal random variables
	// by choosing random points (x, y) in the square -1 < x < 1,
	// -1 < y < 1. Ref: George Marsaglia, Normal (Gaussian) random
	// variables for supercomputers, The Journal of Supercomputing,
	// Volume 5, Number 1, 49 55, DOI: 10.1007/BF00155857.
	// The below implemention is performed by Yang according to
	// http://en.wikipedia.org/wiki/Marsaglia_polar_method
	if (usey)
	{
		L200: x = ran2(idum) * 2. - 1.;
		y = ran2(idum) * 2. - 1.;
		s = x * x + y * y;
		if (s <= 0. || s >= 1.)
		{
			goto L200;
		}
		root = sqrt(-2. / s * log(s));
		ranbyx = x * root;
		ranbyy = y * root;
		ret_val = ranbyy;
		usey = 0;
		// next time use x
	}
	else
	{
		ret_val = ranbyx;
		usey = 1;
		// next time recalculate ranByX and ranByY
	}
	return ret_val;
}

double Utils::yyGaml(double *xx)
{
	return lgammafn(*xx);
}


