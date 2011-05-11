// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#########################################################################################

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>


//#define GENSA_DBG
#define TRACEMAT_NBDATA		6
#define NBOUND  300000
#define IBOUND 10000
#define NMAX 10000
#define NSTPMX 300000
#define nmaxPara_x 100000
#define TRACEMATSIZE 6000000

#ifndef TRUE_
#define TRUE_ 1
#endif

#ifndef FALSE_
#define FALSE_ 0
#endif

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

double d_mod(double *x, double *y)
{
	double q;
	if( (q = *x / *y) >= 0)
		q = floor(q);
	else
		q = -floor(-q);
	return(*x - (*y) * q );
}




typedef struct opt_struct
{
	SEXP R_fn; /* function */
	SEXP R_jc; /* judge constraint function */
	SEXP R_env; /* where to evaluate the calls */
	SEXP xNames; /* names for x */
	SEXP paramNames; /* names for param */
	SEXP params;
	int verbose;
} opt_struct, *OptStruct;

int isVerbose(void* ex)
{
	OptStruct OS = (OptStruct) ex;
	return OS->verbose;
}

