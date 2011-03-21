// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#########################################################################################
#include "GenSAFunction.h"

//================================================================================================
/* mainGSAFun_.f -- translated by f2c (version 20090411).
 You must link the resulting object file with libf2c:
 on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm
 or, if you install libf2c.a in a standard place, with -lf2c -lm
 -- in that order, at the end of the command line, as in
 cc *.o -lf2c -lm
 Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

 http://www.netlib.org/f2c/libf2c.zip
 */
// #include "f2c.h"

/* Common Block Declarations */

struct
{
	double xlow[IBOUND], xhigh[IBOUND], putin1[NBOUND], putin2[NBOUND],
			xrange[IBOUND];
} xy_;

#define xy_1 xy_

struct myconstants_1_
{
	double e_up_tot__, s_up_tot__, k2, kd, i_up_tot__;
	long int nbig;
};

#define myconstants_1 (*(struct myconstants_1_ *) &myconstants_)

struct
{
	long int num_x__, nputin;
} use_func__;

#define use_func__1 use_func__

struct
{
	long int num_function_call__;
} function_call__;

#define function_call__1 function_call__

struct
{
	long int have_constraint__;
} constraint_;

#define constraint_1 constraint_

struct
{
	double para_x__[nmaxPara_x];
	long int npara_x__;
} com_para_x__;

#define com_para_x__1 com_para_x__

struct
{
	long int ls;
} com_ls__;

#define com_ls__1 com_ls__

struct
{
	double x_real_energy__, pseudo_real_energy__,
			x_stop_errorpercent_accordingtoenergy__;
	long int x_know_real_energy__;
} yypw_energy__;

#define yypw_energy__1 yypw_energy__

struct
{
	// was [300000][3000000000]
	double xx[NMAX], y[NSTPMX] /* was [10000][300000] */;
} path_;

#define path_1 path_

struct
{
	double pcom[IBOUND], xicom[IBOUND];
	long int ncom;
} f1com_;

#define f1com_1 f1com_

/* Initialized data */

struct
{
	double e_1[5];
	long int e_2;
} myconstants_ = { 20., 1e6, 354.78, 3.3e5, 800., 2 };



/* This procedure is programed by Yang Xiang for fitting */

/* Subroutine */int maingsafun_(long int *idum1, long int *idum, long int *n,
		double *xmini, double *emini, double *tracemat,
		long int * funin_num_function_call__, long int *istep, long int *interval,
		long int *know_real_energy__,
		double * stop_errorpercent_accordingtoenergy__,
		double *real_energy__, long int *funin_have_constraint__,
		double *lb, double *ub, long int *funin_npara_x__,
		double *funin_para_x__, double * x_initial__,
		long int *funin_ls__, double *temsta, double *qv,
		double *qa, long int *lmarkov, void* ex)
{
	/* Initialized data */

	static double real_x__[IBOUND] =
	{ 577.339299, .00832700907, .109027225 };
	static double xstop_error__ = 1e-6;

	/* System generated locals */
	long int i__1, i__2, i__3;
	double d__1;



	/* Local variables */
	static long int i_putin1__[NBOUND];
	static double a, b;
	extern double func_ave_error_x__(long int *, double *, double *);
	static long int i__, j;
	extern long int judge_constraint__(long int *, double *, void* ex);
	static double r__, s, x[IBOUND];
	static long int ind_trace__;
	extern /* Subroutine */int fsavemini_(long int *, double *,
			double *, double *, double *);
	static double s1, t1;
	extern /* Subroutine */int ls_energy__(long int *, long int *, double *,
			double *, void*);
	static double t2;
	static long int ip, it;
	static double fx[IBOUND];
	static long int useuserini;
	static double pqa;
	static long int i_read_real__;
	static double tem, tmp;
	extern double ran2_(long int *);
	static double pqa1;
	static long int know_real_x__;
	static double ave_error_x__, etot;
	static long int just;
	extern /* Subroutine */int save1_(long int *, double *, double *),
			save2_(long int *, double *, double *);
	static double etot0, temqa;
	static long int itdev, itnew;
	static double visit;
	static long int in_constraint__;
	static long int iorder;
	extern /* Subroutine */int energy_(long int *, long int *, double *,
			double *, void* ex);
	extern double visita_(double *, double *, long int *);
	extern /* Subroutine */int coordin_(long int *, long int *, double *);
	static long int canstop;

	
	/* =-1333333  ! idum1 is used for generating different initi */
	/* =-100377  ! parameter for function input */
	/* =2, ! ! parameter for function input */

	/* funIn_num_function_call is an long int */
	/* =5000, ! ! parameter for function input */
	/* =1, ! ! parameter for function input */
	/* = .false. ! parameter for function input */
	/* =0.01d0, ! ! parameter for */
	/* =0.d0 ! parameter for function input */
	/* funIn_have_constraint is a long int variab */
	/* ! (1)=0.d0 ! parameter for function input */
	/* ! (2)=10.d0 ! parameter for function input */
	/* =1 ! funIn_nPara_x is an long int */
	/* funIn_para_x(nmaxPara_x) is a real*8 matrix */

	/* funIn_LS is an long int */
	/* =3000.d0 ! parameter for function input */
	/* =2.62d0 ! parameter for function input */
	/* =-5.d0 ! parameter for function input */
	/* Parameter adjustments */
	--x_initial__;
	--funin_para_x__;
	--ub;
	--lb;
	//tracemat -= 100001;
	--xmini;

	/* Function Body */
	/* data lb,ub/0.d0,0.d0,20.d0,20.d0/  ! parameter for function input */
	/* open(11,file='fitz.in',status='old',form='formatted') */
	/* open(16,file='fitz.in6',status='old',form='formatted') */
	/* open(15,file='fitz.in5',status='old',form='formatted') */
	/* open(12,file='fitz.ou2',status='unknown',form='formatted') */
	/* open(13,file='fitz.ou3',status='unknown',form='formatted') */
	/* open(14,file='fitz.ou4',status='unknown',form='formatted') */
	/* open(17,file='fitz.ou7',status='unknown',form='formatted') */
	/* open(18,file='fitz.ou8',status='unknown',form='formatted') */
	/* open(19,file='fitz.ou9',status='unknown',form='formatted') */
	/* open(20,file='fitz.o20',status='unknown',form='formatted') */
	/* open(21,file='fitz.o21',status='unknown',form='formatted') */
	/* open(22,file='fitz.o22',status='unknown',form='formatted') */
	/*      call dostim(ih,im,is,il) */
	/* assign the variable in function to variable in common block. */

	// Do some initialization
	itnew = 0;

	function_call__1.num_function_call__ = *funin_num_function_call__;

	constraint_1.have_constraint__ = *funin_have_constraint__;
	com_para_x__1.npara_x__ = *funin_npara_x__;
	i__1 = com_para_x__1.npara_x__;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		com_para_x__1.para_x__[i__ - 1] = funin_para_x__[i__];
	}
	com_ls__1.ls = *funin_ls__;
	/* below is for early stop in YYPW. */
	yypw_energy__1.x_know_real_energy__ = *know_real_energy__;
	yypw_energy__1.x_real_energy__ = *real_energy__;
	yypw_energy__1.pseudo_real_energy__ = abs(*real_energy__);
	/* later will still process */
	yypw_energy__1.x_stop_errorpercent_accordingtoenergy__
			= *stop_errorpercent_accordingtoenergy__;
	if (*funin_ls__ == 1 && *lmarkov % *n != 0)
	{
		Rprintf("LMarkov should be size of 'x' (recommended) or 2*n or 3*n ... since component.change is 1\n") ;
		return -1;
	}
	if (isVerbose(ex))
	{
		Rprintf("LMarkov= %f\n", *lmarkov) ;
	}
	ind_trace__ = -1;
	/* idum1=-1333333  ! idum1 is used for generating different initial */
	/*     read(11,*)just */
	/* ind_trace is for traceMat */
	just = 1;
	/*     read(16,*)n,nputin,istep,interval,know_real_x,have_constraint */
	/* n=2 ! parameter for function input */
	/* just could be parameter for function input */
	use_func__1.nputin = 2;
	/* istep=5000 ! parameter for function input */
	/* interval=1 ! parameter for function input */
	know_real_x__ = FALSE_;
	/* know_real_energy = .false. ! parameter for function input */
	/* stop_errorPercent_accordingToEnergy=0.01d0 ! parameter for functi */
	/* real_energy=0.d0 ! parameter for function input */
	/* have_constraint=.false.  ! parameter for function input */
	/* lb(1)=0.d0 ! parameter for function input */
	/* lb(2)=0.d0 ! parameter for function input */
	/* ub(1)=10.d0 ! parameter for function input */
	/* ub(2)=10.d0 ! parameter for function input */
	/* nPara_x=1 ! parameter for function input */
	/* para_x(1)=2.d0 ! parameter for function input */
	/* x_initial(1)=2.33d0 */
	/* x_initial(2)=2.33d0 */
	/* LS=2 */
	/* never use real_x to judge stop since there m */
	if (*know_real_energy__)
	{
		/* get pseudo_real_energy which is positi */
		yypw_energy__1.pseudo_real_energy__ = abs(*real_energy__);
		if (abs(yypw_energy__1.pseudo_real_energy__) <= 1e-10)
		{
			yypw_energy__1.pseudo_real_energy__ = 1e-5;
			/* if contain 0, substitute 0 wi */
		}
	}
	use_func__1.num_x__ = *n;
	if (*n > IBOUND)
	{
		if (isVerbose(ex))
		{
			Rprintf("number of parameters is over, please increase ibound in the program\n") ;
		}
		return -1;
	}
	if (use_func__1.nputin > NBOUND)
	{
		if (isVerbose(ex))
		{
			Rprintf("number of fitting data is over, please increase nbound in the program\n");
		}
		return -1;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		/* read(16,*)xlow(i),xhigh(i) */
		xy_1.xlow[i__ - 1] = lb[i__];
		xy_1.xhigh[i__ - 1] = ub[i__];
		/* L300: */
	}
	for (i__ = 1; i__ <= IBOUND; ++i__)
	{
		xy_1.xrange[i__ - 1] = xy_1.xhigh[i__ - 1] - xy_1.xlow[i__ - 1];
		/* L302: */
	}
	i_read_real__ = 2;
	/* means don't add normal noise. */
	if (isVerbose(ex))
	{
		Rprintf("xrange: ") ;
		for (ip = 1; ip <= 2; ++ip)
		{
			Rprintf("%f ", xy_1.xrange[ip - 1]) ;
		}
	}
	i__1 = just;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = use_func__1.nputin;
		for (j = 1; j <= i__2; ++j)
		{
			/* waste=YYGas(idum1) ! to run YYGas just*nputin times */
		}
	}
	if (i_read_real__ == 1)
	{
		i__1 = use_func__1.nputin;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			/* read(15,*)putin1(i),putin2(i) */
			/* putin2(i)=putin2(i)+YYGas(idum1)*2.d0*0.01d0*40000d0 ! add n */
			/* write(18,1900)putin1(i),putin2(i) */
		}
	}
	else
	{
		i__1 = use_func__1.nputin;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			/* read(15,*)I_putin1(i),putin2(i) */
			i_putin1__[i__ - 1] = 1;
			xy_1.putin2[i__ - 1] = 2.;
			xy_1.putin1[i__ - 1] = (double) i_putin1__[i__ - 1];
			/* write(18,1800)I_putin1(i),putin2(i) */
		}
	}
	for (i__ = use_func__1.nputin + 1; i__ <= NBOUND; ++i__)
	{
		xy_1.putin1[i__ - 1] = 0.;
		xy_1.putin2[i__ - 1] = 0.;
		/* L201: */
	}
	/* idum=-100377  ! parameter for function input */
	/* idum1=-1333333 ! reset idum1 to make fitz5.f is the same as fitz4 */
	/* temsta=3000.d0 ! parameter for function input */
	/* qv=2.62d0 ! parameter for function input */
	/* qa=-5.d0 ! parameter for function input */
	/* use the initial values provided by the user */
	useuserini = TRUE_;
	if (useuserini)
	{
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			x[i__ - 1] = x_initial__[i__];
			/* parameter for function input */
		}
	}
	in_constraint__ = TRUE_;
	/* continue here, 18:34, 2010, 10, 21 */
	/* set in_constraint=.true. to make us not to */
	L400: if (!in_constraint__)
	{
		i__1 = just;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			coordin_(n, idum1, x);
			/* L10: */
		}
	}
	if (constraint_1.have_constraint__)
	{
		/* if have constraint */
		in_constraint__ = judge_constraint__(n, x, ex);
		/* judge_constraint(n,x) a gi */
		if (!in_constraint__)
		{
			goto L400;
		}
		/* if x are not in constraint */
	}
	if (isVerbose(ex))
	{
		Rprintf("Now this program start: just=%d\n", just) ;
		Rprintf("The random intial x coordinates:\n") ;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			Rprintf("%f ", x[j - 1]) ;
		}
	}
	energy_(n, &use_func__1.nputin, &etot, x, ex);
	if (isVerbose(ex))
	{
		Rprintf("The energy of initial x = %f\n", etot) ;
	}
	if (etot != etot)
	{
		/* judge if etot is NAN */
		if (isVerbose(ex))
		{
			i__1 = *n;
			for (j = 1; j <= i__1; ++j)
			{
				Rprintf("%f ", x[j - 1]) ;
			}
			Rprintf(" give NaN, generate new starting point\n") ;
		}
		goto L400;
	}
	fsavemini_(n, &etot, emini, x, &xmini[1]);
	/* record below variables. */
	++ind_trace__;
	//tracemat[ind_trace__ + 100000] = (double) ind_trace__;
	//tracemat[ind_trace__ + 200000] = (double) itnew;
	//tracemat[ind_trace__ + 300000] = tem;
	//tracemat[ind_trace__ + 400000]
	//		= (double) function_call__1.num_function_call__;
	//tracemat[ind_trace__ + 500000] = etot0;
	//tracemat[ind_trace__ + 600000] = *emini;

	// Initialize etot0 and temp
	etot0 = etot;
	tem = *temsta;

	tracemat[TRACEMAT_NBDATA * ind_trace__] = (double) ind_trace__;
	tracemat[TRACEMAT_NBDATA * ind_trace__ + 1] = (double) itnew;
	tracemat[TRACEMAT_NBDATA * ind_trace__ + 2] = tem;
	tracemat[TRACEMAT_NBDATA * ind_trace__ + 3]
			= (double) function_call__1.num_function_call__;
	tracemat[TRACEMAT_NBDATA * ind_trace__ + 4] = etot0;
	tracemat[TRACEMAT_NBDATA * ind_trace__ + 5] = *emini;

	if (com_ls__1.ls == 20020121)
	{
		if (isVerbose(ex))
		{
			Rprintf("LS is 20020121\n") ;
		}
		ls_energy__(n, &use_func__1.nputin, &etot, x, ex);
	}
	else
	{
		/* call energy(n,nputin,etot,x) ! don't want to repeat energy cal. */
		if (isVerbose(ex))
		{
			Rprintf("do not want to repeat energy call before start\n") ;
		}
	}
	if (isVerbose(ex))
	{
		Rprintf("The transformed xinitial x: \n") ;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			Rprintf("%f ", x[j - 1]) ;
		}
		Rprintf("\n") ;
		Rprintf("The energy of  transformed initial x = %f\n", etot) ;
	}
	etot0 = etot;
	if (etot < *emini)
	{
		fsavemini_(n, &etot, emini, x, &xmini[1]);
	}
	/* write(14,*)itnew,' ',tem,' ',etot0,' ',Emini */
	if (isVerbose(ex))
	{
		Rprintf("At the beginning, etot0= %f\n", etot0) ;
		Rprintf("Emini= %f\n", emini) ;
		Rprintf("Current x:\n") ;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			Rprintf("%f ", x[j - 1]) ;
		}
		Rprintf("\nCurrent xmini:\n") ;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			Rprintf("%f ", xmini[j]) ;
		}
		Rprintf("\n") ;
	}
	if (etot < *emini)
	{
		fsavemini_(n, &etot, emini, x, &xmini[1]);
	}
	/* ===> If you know the real x, below give a stop criterion according to the average */
	/* ===> error between current x and real x. */
	if (know_real_x__)
	{
		ave_error_x__ = func_ave_error_x__(n, x, real_x__);
		if (ave_error_x__ < xstop_error__)
		{
			if (isVerbose(ex))
			{
				Rprintf("Get accurate x with precison < %f\n", xstop_error__) ;
			}
			goto L2000;
		}
	}
	/* If you know the real_energy, below give a stop criterion acc */
	/* error between current energy and Emini. */
	if (*know_real_energy__)
	{
		canstop = (d__1 = *emini - *real_energy__, abs(d__1))
				<= yypw_energy__1.pseudo_real_energy__
						* *stop_errorpercent_accordingtoenergy__ * .01;
		if (canstop)
		{
			if (isVerbose(ex))
			{
				Rprintf("Have got accurate energy with relative precison < %f\n", *stop_errorpercent_accordingtoenergy__);
			}
			goto L2000;
		}
	}
	i__1 = *istep;
	for (it = 1; it <= i__1; ++it)
	{
		itnew = it;
		s1 = (double) itnew;
		s = s1 + 1.;
		t1 = exp((*qv - 1.) * log(2.)) - 1.;
		t2 = exp((*qv - 1.) * log(s)) - 1.;
		tem = *temsta * t1 / t2;
		temqa = tem / (double) itnew;
		/* 111     continue */
		/* do 60 j=1,n ! here j used for markov chain length, and order of */
		i__2 = *lmarkov;
		for (j = 1; j <= i__2; ++j)
		{
			/* default of LMarkov is n. */
			save1_(n, x, fx);
			/* save x to fx every time when try to chan */
			L111:
			/* LS: 1. one component change at one time; 2. all components c */
			/* 20020121. all components change at one time and YYPW method */
			if (com_ls__1.ls == 2 || com_ls__1.ls == 20020121)
			{
				/* change all the components of x at same time. */
				i__3 = *n;
				for (iorder = 1; iorder <= i__3; ++iorder)
				{
					/* iorder means the order of component in x */
					visit = visita_(qv, &tem, idum);
					x[iorder - 1] = visit + fx[iorder - 1];
					a = x[iorder - 1] - xy_1.xlow[iorder - 1];
					b = d_mod(&a, &xy_1.xrange[iorder - 1])
							+ xy_1.xrange[iorder - 1];
					x[iorder - 1] = d_mod(&b, &xy_1.xrange[iorder - 1])
							+ xy_1.xlow[iorder - 1];
					tmp = d_mod(&b, &xy_1.xrange[iorder - 1]);
					/*              write(*,*)'after, iorder:',iorder,', x are',(x(ip),ip=1,n)! remove */
					/*     $         ,', visit: ',visit,', a=',a,', b:',b,' tmp=',tmp! remove */
				}
			}
			else
			{
				visit = visita_(qv, &tem, idum);
				x[j - 1] = visit + fx[j - 1];
				a = x[j - 1] - xy_1.xlow[j - 1];
				b = d_mod(&a, &xy_1.xrange[j - 1]) + xy_1.xrange[j - 1];
				x[j - 1] = d_mod(&b, &xy_1.xrange[j - 1]) + xy_1.xlow[j - 1];
				tmp = d_mod(&b, &xy_1.xrange[j - 1]);
				/*              write(*,*)'after, j:',j,', x are',(x(ip),ip=1,n)! remove */
				/*     $         ,', visit: ',visit,', a=',a,', b:',b,' tmp=',tmp! remove */
			}
			/* 60       continue */
			/*          write(*,*)'\n' */
			/*          write(*,*)'tem: ',tem,', temqa: ',temqa,', itnew: ',itnew, */
			/*     $     ', markov j=',j */
			/*          write(*,*)'fx are',(fx(ip),ip=1,2) */
			/*          write(*,*)'x are',(x(ip),ip=1,2) */
			/*          !read(*,*)iiii */
			if (constraint_1.have_constraint__)
			{
				/* if have constraint */
				in_constraint__ = judge_constraint__(n, x, ex);
				/* judge_constraint(n,x) */
				if (!in_constraint__)
				{
					goto L111;
				}
			}
			if (com_ls__1.ls == 20020121)
			{
				if (isVerbose(ex))
				{
					Rprintf("LS is 20020121\n");
				}
				ls_energy__(n, &use_func__1.nputin, &etot, x, ex);
			}
			else
			{
				energy_(n, &use_func__1.nputin, &etot, x, ex);
			}
			/* call LS_energy(n,nputin,etot,x) */
			if (isVerbose(ex))
			{
				Rprintf("Before judge, etot0= %f etot= %f\n",etot0, etot);
			}
			if (etot < etot0)
			{
				etot0 = etot;
				if (isVerbose(ex))
				{
					Rprintf("etot.lt.etot0\n");
				}
				/* remove */
				if (etot < *emini)
				{
					fsavemini_(n, &etot, emini, x, &xmini[1]);
				}
			}
			else
			{
				r__ = ran2_(idum);
				pqa1 = (*qa - 1.) * (etot - etot0) / temqa + 1.;
				/* write(*,*)' etot0=',etot0,', etot=',etot,', pqa1=',pqa1 ! r */
				if (pqa1 < 0.)
				{
					pqa = 0.;
				}
				else
				{
					pqa = exp(log(pqa1) / (1. - *qa));
				}
				if (isVerbose(ex))
				{
					Rprintf("pqa= %f r= %f \n", pqa, r__) ;
				}
				/* remove */
				if (r__ > pqa)
				{
					save2_(n, x, fx);
				}
				else
				{
					etot0 = etot;
				}
			}
			/* record below variables. */
			++ind_trace__;
			//tracemat[ind_trace__ + 100000] = (double) ind_trace__;
			//tracemat[ind_trace__ + 200000] = (double) itnew;
			//tracemat[ind_trace__ + 300000] = tem;
			//tracemat[ind_trace__ + 400000]
			//		= (double) function_call__1.num_function_call__;
			//tracemat[ind_trace__ + 500000] = etot0;
			//tracemat[ind_trace__ + 600000] = *emini;
			tracemat[TRACEMAT_NBDATA * ind_trace__] = (double) ind_trace__;
			tracemat[TRACEMAT_NBDATA * ind_trace__ + 1] = (double) itnew;
			tracemat[TRACEMAT_NBDATA * ind_trace__ + 2] = tem;
			tracemat[TRACEMAT_NBDATA * ind_trace__ + 3]
					= (double) function_call__1.num_function_call__;
			tracemat[TRACEMAT_NBDATA * ind_trace__ + 4] = etot0;
			tracemat[TRACEMAT_NBDATA * ind_trace__ + 5] = *emini;

			/* ===>     If you know the real x, below give a stop criterion according to the average */
			/* ===>     error between current x and real x. */
			if (know_real_x__)
			{
				ave_error_x__ = func_ave_error_x__(n, x, real_x__);
				if (ave_error_x__ < xstop_error__)
				{
					if (isVerbose(ex))
					{
						Rprintf("get accurate x with precison < %f\n", xstop_error__) ;
					}
					goto L2000;
				}
			}
			/* If you know the real_energy, below give a stop criterion acc */
			/* error between current energy and Emini. */
			if (*know_real_energy__)
			{
				canstop = (d__1 = *emini - *real_energy__, abs(d__1))
						<= yypw_energy__1.pseudo_real_energy__
								* *stop_errorpercent_accordingtoenergy__ * .01;
				if (canstop)
				{
					if (isVerbose(ex))
					{
						Rprintf("Have got accurate energy with relative precison < %f\n", *stop_errorpercent_accordingtoenergy__);
					}
					goto L2000;
				}
			}
			/* L60: */
		}
		itdev = itnew % *interval;
		if (itdev == 0)
		{
			/*          write(14,*)itnew,' ',tem,' ',etot0,' ',Emini */
			if (isVerbose(ex))
			{
				Rprintf("> %d %f %f %f <-- Emini\n", itnew, tem, etot0, *emini );
				Rprintf("Current x:\n");
				i__2 = *n;
				for (j = 1; j <= i__2; ++j)
				{
					Rprintf("%f ", x[j - 1]) ;
				}
				Rprintf("\nCurrent xmini:\n") ;
				i__2 = *n;
				for (j = 1; j <= i__2; ++j)
				{
					Rprintf("%f ", xmini[j]) ;
				}
				Rprintf("\n") ;
			}
		}
		/* L40: */
	}

	/* 2000  write(*,*)'Sum of squred deviation is:',Emini */
	L2000: if (isVerbose(ex))
	{
		Rprintf("Emini is: %f\n",*emini);
		Rprintf("Just is: %d\n", just);
		Rprintf("xmini are:\n");
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			Rprintf("%f ", xmini[j]) ;
		}
		Rprintf("\n") ;
	}
	if (com_ls__1.ls == 20020121)
	{
		/* because LS_energy call function one more time at begining. */
		--function_call__1.num_function_call__;
	}
	if (isVerbose(ex))
	{
		Rprintf("No. of function call is: %d\n", function_call__1.num_function_call__);
	}
	*funin_num_function_call__ = function_call__1.num_function_call__;
	/*      write(17,*)num_function_call */
	/* write(13,*)'just is:',just,'     the fitting parameters are:' */
	/*      write(13,*)(xmini(j),j=1,n),Emini */
	/* write(13,*) */
	/* L100: */
	/* L110: */
	/* L1800: */
	/* L1900: */
	/*      call dostim(ih1,im1,is1,il1) */
	/*       write(12,*)just,etot0 */
	/* just=just+1 */
	/* rewind 11 */
	/*      write(11,*)just */
	/*     write(12,*)ih,im,is,il */
	/*     write(12,*)ih1,im1,is1,il1 */
	return 0;
} /* maingsafun_ */

/* Subroutine */int fsavemini_(long int *n, double *etot,
		double * emini, double *x, double *xmini)
{
	/* System generated locals */
	long int i__1;

	/* Local variables */
	static long int i__;

	/* Parameter adjustments */
	--xmini;
	--x;

	/* Function Body */
	*emini = *etot;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		xmini[i__] = x[i__];
		/* L10: */
	}
	return 0;
} /* fsavemini_ */

double func_ave_error_x__(long int *n, double *x, double *real_x__)
{
	/* System generated locals */
	long int i__1;
	double ret_val, d__1;

	/* Local variables */
	static long int i__;
	static double ave_error_x__, error;

	/* Parameter adjustments */
	--real_x__;
	--x;

	/* Function Body */
	ave_error_x__ = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		error = (d__1 = x[i__] - real_x__[i__], abs(d__1));
		ave_error_x__ += error;
	}
	ret_val = ave_error_x__ / (double)(*n);
	return ret_val;
} /* func_ave_error_x__ */

/* Subroutine */int coordin_(long int *n, long int *idum1, double *x)
{
	/* System generated locals */
	long int i__1;

	/* Local variables */
	static long int j;
	extern double ran2_(long int *);

	/* Parameter adjustments */
	--x;

	/* Function Body */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j)
	{
		x[j] = ran2_(idum1) * xy_1.xrange[j - 1] + xy_1.xlow[j - 1];
		/* L100: */
	}
	return 0;
} /* coordin_ */

long int ori_judge_constraint__(long int *n, double *x)
{
	/* System generated locals */
	long int ret_val;

	/* Local variables */
	static double k3, k_3__, factor1;

	/* judge_constraint(x) a */
	/*      long int judge_constraint */
	/* =====> below is the constrain */
	/* Parameter adjustments */
	--x;

	/* Function Body */
	k3 = x[2];
	k_3__ = x[3];
	factor1 = (k_3__ + k3 * (double) myconstants_1.nbig
			* (myconstants_1.i_up_tot__ + myconstants_1.e_up_tot__)) * (k_3__
			+ k3 * (double) myconstants_1.nbig * (myconstants_1.i_up_tot__
					+ myconstants_1.e_up_tot__)) - k3 * 4. * k3 * (double)
	myconstants_1.nbig * (double) myconstants_1.nbig
			* myconstants_1.i_up_tot__ * myconstants_1.e_up_tot__;
	if (factor1 < 0.)
	{
		ret_val = FALSE_;
	}
	if (factor1 >= 0.)
	{
		ret_val = TRUE_;
	}
	return ret_val;
} /* ori_judge_constraint__ */

/* Subroutine */int energy_(long int *num_x__, long int *nputin,
		double * etot, double *x, void* ex)
{
	/* Initialized data */

	/* System generated locals */
	long int i__1;

	/* Local variables */
	extern double yangfunc_(long int *, long int *, double *,
			double *, void* ex);
	static long int i__;
	extern long int judge_constraint__(long int *, double *, void* ex);
	static double delta_energy__;
	static long int in_constraint__;
	static double penalty;

	/* may be deleted for different energy */
	/* may be deleted for different en */
	/* Parameter adjustments */
	--x;

	/* Function Body */
	if (constraint_1.have_constraint__)
	{
		/* if have constraint */
		in_constraint__ = judge_constraint__(num_x__, &x[1], ex);
		/* judge_constraint(n,x) */
		if (!in_constraint__)
		{
			*etot = 1e10;
			/*          write(*,*)'in ',number_function_call+1,' times calling energy' */
			/*     $ ,'subtoutine, x not in constraint, return very big etot = ',etot */
			return 0;
		}
	}
	/* ######> below could be changed according to different fitting function. */
	/* below is the function provided by user. */
	/* etot=userFunc(x,para_x) */
	*etot = yangfunc_(num_x__, nputin, etot, &x[1], ex);
	/* #####> the above can be changed according to different fitting function. */
	if (*etot != *etot)
	{
		/* judge if etot is NaN */
		if (function_call__1.num_function_call__ == 1)
		{
			/* if initial NaN, don't change i */
		}
		else
		{
			*etot = function_call__1.num_function_call__ * 1e5 + 1e10;
			/* write(*,*)'num_function_call= ',num_function_call+1 */
			/*          write(*,*)num_function_call+1,(x(j),j=1,num_x), */
			/*     $     '->NAN, etot= ',etot */
			/* add num_function_call to m */
		}
	}
	/* =====> below give penalty */
	penalty = 0.;
	i__1 = *num_x__;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		if (x[i__] >= xy_1.xlow[i__ - 1] && x[i__] <= xy_1.xhigh[i__ - 1])
		{
			delta_energy__ = 0.;
		}
		else
		{
			if (x[i__] < xy_1.xlow[i__ - 1])
			{
				delta_energy__ = (x[i__] - xy_1.xlow[i__ - 1]) * 1e11 * (x[i__]
						- xy_1.xlow[i__ - 1]);
			}
			if (x[i__] > xy_1.xhigh[i__ - 1])
			{
				delta_energy__ = (x[i__] - xy_1.xhigh[i__ - 1]) * 1e11
						* (x[i__] - xy_1.xhigh[i__ - 1]);
			}
			/* write(*,*)'out of boundary x: ',(x(j),j=1,num_x) */
			/* write(*,*)'delta_energy= ',delta_energy */
		}
		penalty += delta_energy__;
	}
	*etot += penalty;
	function_call__1.num_function_call__ += 1;
	/* L200: */
	/* L300: */
	return 0;
} /* energy_ */

/* Subroutine */int ls_energy__(long int *num_x__, long int *nputin,
		double *etot, double *x, void* ex)
{
	/* System generated locals */
	long int i__1, i__2;

	/* Local variables */
	static long int i__, j, n;
	static double p[10000];
	static long int np;
	static double xi[100000000] /* was [10000][10000] */, fret, ftol;
	//static long int iter;
	extern /* Subroutine */int yypw_(double *, double *, long int *,
			long int *, double *, long int *, double *, void*);

	/*      REAL*8 my_parameter(NMAX),vstart(NMAX),xx(NSTPMX),y(NMAX,NSTPMX), */
	/*     $ production(NSTPMX),curve(NSTPMX) */
	/*      real*8 k3,k_3,k2,Kd,I_up_tot */
	/*      EXTERNAL derivs */
	/*      COMMON /path/ xx,y */
	/*      common/myconstants/E_up_tot,S_up_tot,k2,Kd,I_up_tot,Nbig */
	/* U    USES rk4 */
	/*      REAL*8 h,x,dv(NMAX),v(NMAX) */
	/* below is the initial setting of YYPW */
	/* Parameter adjustments */
	--x;

	/* Function Body */
	n = *num_x__;
	np = 10000;
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		p[i__ - 1] = x[i__];
	}
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = n;
		for (j = 1; j <= i__2; ++j)
		{
			xi[i__ + j * 10000 - 10001] = 0.;
		}
	}
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		xi[i__ + i__ * 10000 - 10001] = 1.;
	}
	ftol = 1e-10;
	//yypw_(p, xi, &n, &np, &ftol, &iter, &fret, ex);
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		x[i__] = p[i__ - 1];
	}
	*etot = fret;
	return 0;
} /* ls_energy__ */

/* Subroutine */int save1_(long int *n, double *x, double *fx)
{
	/* System generated locals */
	long int i__1;

	/* Local variables */
	static long int i__;

	/* Parameter adjustments */
	--fx;
	--x;

	/* Function Body */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		fx[i__] = x[i__];
		/* L10: */
	}
	return 0;
} /* save1_ */

/* Subroutine */int save2_(long int *n, double *x, double *fx)
{
	/* System generated locals */
	long int i__1;

	/* Local variables */
	static long int i__;

	/* Parameter adjustments */
	--fx;
	--x;

	/* Function Body */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		x[i__] = fx[i__];
		/* L10: */
	}
	return 0;
} /* save2_ */

double visita_(double *q, double *temp, long int *idum)
{
	/* System generated locals */
	double ret_val, d__1;

	/* Builtin functions */
	double asin( double), log( double), exp( double), sqrt(
			double), sin( double);

	/* Local variables */
	static double x, y, pi, den;
	extern double yygas_(long int *);
	static double fator1, fator2, fator3, fator4, fator5, fator6, sigmax;
	extern double yygaml_(double *);

	pi = asin(1.) * 2.;
	fator1 = exp(log(*temp) / (*q - 1.));
	fator2 = exp((4. - *q) * log(*q - 1.));
	fator3 = exp((2. - *q) * log(2.) / (*q - 1.));
	fator4 = sqrt(pi) * fator1 * fator2 / (fator3 * (3. - *q));
	fator5 = 1. / (*q - 1.) - .5;
	/*      calculates the gamma function using the reflection formula for */
	/*      0<arg<1 */
	d__1 = 2. - fator5;
	fator6 = pi * (1. - fator5) / sin(pi * (1. - fator5)) / exp(yygaml_(&d__1));
	sigmax = exp(-(*q - 1.) * log(fator6 / fator4) / (3. - *q));
	x = sigmax * yygas_(idum);
	y = yygas_(idum);
	den = exp((*q - 1.) * log((abs(y))) / (3. - *q));
	ret_val = x / den;
	return ret_val;
} /* visita_ */

double yygas_(long int *idum)
{
	//    /* Initialized data */
	//
	//    static long int usey = TRUE_;
	//
	//    /* System generated locals */
	//    double ret_val;
	//
	//    /* Builtin functions */
	//    double log(double), sqrt(double);
	//
	//    /* Local variables */
	//    static double s, x, y;
	//    extern double ran2_(long int *);
	//    static double ranbyx, ranbyy;
	//
	///* This function YYGas use the polar method from George Marsaglia t */
	///* generate a pair of independent standard normal random variables */
	///* by choosing random points (x, y) in the square -1 < x < 1, */
	///* -1 < y < 1. Ref: George Marsaglia, Normal (Gaussian) random */
	///* variables for supercomputers, The Journal of Supercomputing, */
	///* Volume 5, Number 1, 49�55, DOI: 10.1007/BF00155857. */
	///* The below implemention is performed by Yang according to */
	///* http://en.wikipedia.org/wiki/Marsaglia_polar_method */
	//    if (usey) {
	//L100:
	//	x = ran2_(idum) * 2. - 1.;
	//	y = ran2_(idum) * 2. - 1.;
	//	s = x * x + y * y;
	//	if (s >= 1. || s <= 0.) {
	//	    goto L100;
	//	}
	//	ranbyx = x * sqrt(log(s) * -2. / s);
	//	ranbyy = y * sqrt(log(s) * -2. / s);
	//	ret_val = ranbyy;
	//	usey = FALSE_;
	///* next time use x */
	//    } else {
	//	ret_val = ranbyx;
	//	usey = TRUE_;
	///* next time recalculate ranByX and ranByY */
	//    }
	//    return ret_val;

	/* Initialized data */

	static long int usey = TRUE_;

	/* System generated locals */
	double ret_val;

	/* Builtin functions */
	double log( double), sqrt( double);

	/* Local variables */
	static double s, x, y;
	extern double ran2_(long int *);
	static double root, ranbyx, ranbyy;

	/* This function YYGas use the polar method from George Marsaglia t */
	/* generate a pair of independent standard normal random variables */
	/* by choosing random points (x, y) in the square -1 < x < 1, */
	/* -1 < y < 1. Ref: George Marsaglia, Normal (Gaussian) random */
	/* variables for supercomputers, The Journal of Supercomputing, */
	/* Volume 5, Number 1, 49�55, DOI: 10.1007/BF00155857. */
	/* The below implemention is performed by Yang according to */
	/* http://en.wikipedia.org/wiki/Marsaglia_polar_method */
	if (usey)
	{
		L100: x = ran2_(idum) * 2. - 1.;
		y = ran2_(idum) * 2. - 1.;
		s = x * x + y * y;
		if (s >= 1. || s <= 0.)
		{
			goto L100;
		}
		root = sqrt(log(s) * -2. / s);
		ranbyx = x * root;
		ranbyy = y * root;
		ret_val = ranbyy;
		usey = FALSE_;
		/* next time use x */
	}
	else
	{
		ret_val = ranbyx;
		usey = TRUE_;
		/* next time recalculate ranByX and ranByY */
	}
	return ret_val;

} /* yygas_ */

double yygaml_(double *xx)
{
	return lgammafn(*xx);
} /* yygaml_ */

double ran2_(long int *idum)
{
	/* Initialized data */

	static long int idum2 = 123456789;
	static long int iv[32] =
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0 };
	static long int iy = 0;

	/* System generated locals */
	long int i__1;
	double ret_val, d__1;

	/* Local variables */
	static long int j, k;

	/*        this is a random number generator. See the article on */
	/*        p522 of COMPUTERS IN PHYSICS, SEP/OCT 1992, in which the */
	/*        authors claim that ran2 is a 'perfect' random generator. */

	/*        Long period(>2x10^18) random generator of L'Ecuyer with */
	/*        Bays-Durhamc shuffle and added safeguards. Returns a uniform */
	/*        random deviate between 0.0 and 1.0 (exclusive of the endpoint */
	/*        values). Call with idum a negative long int to initialize; */
	/*        therefore, do not alter idum between successive deviates in a */
	/*        sequence. RNMX should approximate the largest floating value */
	/*        that is less than 1.0 */

	if (*idum <= 0)
	{
		/* Computing MAX */
		i__1 = -(*idum);
		*idum = max(i__1, 1);
		/* be sure to prevent idum=0 */
		idum2 = *idum;
		for (j = 40; j >= 1; --j)
		{
			/* load the shuffle table (after 8 warm-up */
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
	/* start here when no initializing */
	*idum = (*idum - k * 53668) * 40014 - k * 12211;

	/*        compute idum=mod(IA1*idum,IM1) without overflows by Schrage's */
	/*        method */

	if (*idum < 0)
	{
		*idum += 2147483563;
	}
	k = idum2 / 52774;
	idum2 = (idum2 - k * 52774) * 40692 - k * 3791;

	/*        compute idum2=mod(IA2*idum2,IM2) likewise */

	if (idum2 < 0)
	{
		idum2 += 2147483399;
	}
	j = iy / 67108862 + 1;
	/* will be in the range 1:NTAB */
	iy = iv[j - 1] - idum2;

	/*        here idum is shuffled, idum and idum2 are combined to grnerate */
	/*        output */

	iv[j - 1] = *idum;
	if (iy < 1)
	{
		iy += 2147483562;
	}
	/* Computing MIN */
	d__1 = iy * 4.6566130573917691e-10;
	ret_val = min(d__1, .99999987999999995);

	/* because users don't expect endpoints */
	return ret_val;
} /* ran2_ */

double func_(double *p, void* ex)
{
	/* System generated locals */
	long int i__1;
	double ret_val;

	/* Local variables */
	static long int i__;
	static double x[IBOUND], etot;
	extern /* Subroutine */int energy_(long int *, long int *, double *,
			double *, void* ex);

	/* Parameter adjustments */
	--p;

	/* Function Body */
	i__1 = use_func__1.num_x__;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		x[i__ - 1] = p[i__];
	}
	energy_(&use_func__1.num_x__, &use_func__1.nputin, &etot, x, ex);
	ret_val = etot;
	return ret_val;
} /* func_ */

double yangfunc_(long int * xSize, long int * notused1,
		double * notused2, double * xData, void* ex)
{
	SEXP x4R, val;
	int i;
	double res = 0;
	OptStruct OS = (OptStruct) ex;
	//PROTECT_INDEX ipx;

	// Allocate vector for R which is size of the vector in the R context.
	PROTECT(x4R = allocVector(REALSXP, *xSize));
	if (!isNull(OS->xNames))
		setAttrib(x4R, R_NamesSymbol, OS->xNames);

	for (i = 0; i < *xSize; i++)
	{
		if (!R_FINITE(xData[i]))
		{
			Rprintf("non-finite value supplied by GenSAFunction");
			xData[i] = 0;
		}
		REAL(x4R)[i] = xData[i];
	}

	SETCADR(OS->R_fn, x4R);
	val = eval(OS->R_fn, OS->R_env);
	res = REAL(val)[0];
	UNPROTECT(1);

	return res;
}

long int judge_constraint__(long int * xSize, double * xData, void* ex)
{
	SEXP x4R, val;
	int i;
	int res;
	OptStruct OS = (OptStruct) ex;
	//PROTECT_INDEX ipx;

	// Allocate vector for R which is size of the vector in the R context.
	PROTECT(x4R = allocVector(REALSXP, *xSize));
	if (!isNull(OS->xNames))
		setAttrib(x4R, R_NamesSymbol, OS->xNames);

	for (i = 0; i < *xSize; i++)
	{
		if (!R_FINITE(xData[i]))
		{
			Rprintf("non-finite value supplied by GenSAFunction");
			xData[i] = 0;
		}

		REAL(x4R)[i] = xData[i];
	}

	SETCADR(OS->R_jc, x4R);
	val = eval(OS->R_jc, OS->R_env);
	res = LOGICAL(val)[0];
	UNPROTECT(1);

	return res;
}

