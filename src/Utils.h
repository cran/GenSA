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
#include <R_ext/Linpack.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <map>
#include <string>

//#define GENSA_DBG

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define FALSE_ 0
#define TRUE_ 1
static int c__1 = 1;
static int c__11 = 11;
static void timer(double * ttime)
{
    *ttime = 0.0;
}


typedef std::vector<double> dVec;
typedef std::vector<double>::iterator dVecIt;
typedef std::vector<std::string> strVec;
typedef std::vector<std::string>::const_iterator strVecIt;
typedef std::map<std::string, dVec> Map;
typedef std::map<std::string, dVec>::iterator MapIt;

typedef struct opt_struct
{
	SEXP R_fn; /* function */
	SEXP R_jc; /* judge constraint function */
	SEXP R_env; /* where to evaluate the calls */
	SEXP xNames; /* names for x */
	int verbose;
} opt_struct, *OptStruct;

static void mainlb(int n, int m, double *x,
		   double *l, double *u, int *nbd, double *f, double *g,
		   double factr, double *pgtol, double *ws, double * wy,
		   double *sy, double *ss, double *wt, double *wn,
		   double *snd, double *z, double *r, double *d,
		   double *t, double *wa, int *indx, int *iwhere,
		   int *indx2, char *task, int iprint,
		   char *csave, int *lsave, int *isave, double *dsave) ;

static void active(int n, double *l, double *u,
		   int *nbd, double *x, int *iwhere, int iprint,
		   int *prjctd, int *cnstnd, int *boxed) ;

static void bmv(int m, double *sy, double *wt,
		int *col, double *v, double *p, int *info);

static void cauchy(int n, double *x, double *l, double *u, int *nbd,
		   double *g, int *iorder, int * iwhere, double *t,
		   double *d, double *xcp, int m,
		   double *wy, double *ws, double *sy, double *wt,
		   double *theta, int *col, int *head, double *p,
		   double *c, double *wbp, double *v, int *nint,
		   int iprint, double *sbgnrm, int *info, double * epsmch);

static void cmprlb(int n, int m, double *x,
		   double *g, double *ws, double *wy, double *sy,
		   double *wt, double *z, double *r, double *wa,
		   int *indx, double *theta, int *col, int *head,
		   int *nfree, int *cnstnd, int *info);

static void errclb(int n, int m, double factr, double *l, double *u,
		   int *nbd, char *task, int *info, int *k);

static void formk(int n, int *nsub, int *ind, int * nenter, int *ileave,
		  int *indx2, int *iupdat, int * updatd, double *wn,
		  double *wn1, int m, double *ws, double *wy, double *sy,
		  double *theta, int *col, int *head, int *info);

static void formt(int m, double *wt, double *sy, double *ss,
		  int *col, double *theta, int *info);

static void freev(int n, int *nfree, int *indx,
		  int *nenter, int *ileave, int *indx2, int *iwhere,
		  int *wrk, int *updatd, int *cnstnd, int iprint,
		  int *iter);

static void hpsolb(int n, double *t, int *iorder, int iheap);

static void lnsrlb(int n, double *l, double *u,
		   int *nbd, double *x, double *f, double *fold,
		   double *gd, double *gdold, double *g, double *d,
		   double *r, double *t, double *z, double *stp,
		   double *dnorm, double *dtd, double *xstep,
		   double *stpmx, int *iter, int *ifun, int *iback, int *nfgv,
		   int *info, char *task, int *boxed, int *cnstnd,
		   char *csave, int *isave, double *dsave);

static void matupd(int n, int m, double *ws,
		   double *wy, double *sy, double *ss, double *d,
		   double *r, int *itail, int *iupdat, int *col,
		   int *head, double *theta, double *rr, double *dr,
		   double *stp, double *dtd);

static void projgr(int n, double *l, double *u,
		   int *nbd, double *x, double *g, double *sbgnrm);

static void subsm(int n, int m, int *nsub, int *ind,
		  double *l, double *u, int *nbd, double *x,
		  double *d, double *ws, double *wy, double *theta,
		  int *col, int *head, int *iword, double *wv,
		  double *wn, int iprint, int *info);

static void dcsrch(double *f, double *g, double *stp,
		   /*Chgd: the next five are no longer pointers:*/
		   double ftol, double gtol, double xtol,
		   double stpmin, double stpmax,
		   char *task, int *isave, double *dsave);

static void dcstep(double *stx, double *fx, double *dx,
		   double *sty, double *fy, double *dy, double *stp,
		   double *fp, double *dp, int *brackt, double *stpmin,
		   double *stpmax);

static void pvector(char *title, double *x, int n);

static void prn1lb(int n, int m, double *l, double *u, double *x,
		   int iprint, double epsmch);

static void prn2lb(int n, double *x, double *f, double *g, int iprint,
		   int iter, int nfgv, int nact, double sbgnrm,
		   int nint, char *word, int iword, int iback,
		   double stp, double xstep);

static void prn3lb(int n, double *x, double *f, char *task, int iprint,
		   int info, int iter, int nfgv, int nintol, int nskip,
		   int nact, double sbgnrm, int nint,
		   char *word, int iback, double stp, double xstep,
		   int k);

class Utils
{
public:
	static double dSign(double *a, double *b);
	static double dMod(double *x, double *y);
	static double ran2(long int *idum);
	static double yyGas(long int *idum);
	static double yyGaml(double *xx);
	static void setulb(int n, int m, double *x, double *l, double *u, int *nbd,
			double *f, double *g, double factr, double *pgtol, double *wa,
			int * iwa, char *task, int iprint, int *lsave, int *isave,
			double *dsave);
};

#endif // UTILS_H_
