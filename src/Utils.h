// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#########################################################################################

#ifndef UTILS_H_
#define UTILS_H_

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <map>
#include <string>

//#define GENSA_DBG

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))

typedef std::vector<double> dVec;
typedef std::vector<double>::iterator dVecIt;
typedef std::vector<std::string> strVec;
typedef std::vector<std::string>::const_iterator strVecIt;
typedef std::map<std::string,dVec> Map;
typedef std::map<std::string,dVec>::iterator MapIt;

typedef struct opt_struct
{
	SEXP R_fn; /* function */
	SEXP R_jc; /* judge constraint function */
	SEXP R_env; /* where to evaluate the calls */
	SEXP xNames; /* names for x */
	int verbose;
} opt_struct, *OptStruct;


class Utils
{
public:
	static double dSign(double *a, double *b) ;
	static double dMod(double *x, double *y) ;
	static double ran2(long int *idum) ;
	static double yyGas(long int *idum) ;
	static double yyGaml(double *xx) ;
};


#endif // UTILS_H_
