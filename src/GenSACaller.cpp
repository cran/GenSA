// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#########################################################################################
#include "GenSACaller.h"

// Mixed C++/C code for GenSAFunction call (converted from fortran to C)

GenSACaller::GenSACaller()
{
	// Allocate memory for xMini
	// xmini is the vector which contains the results of minimization
	xMini_ = 0;
	traceMat_ = 0;

	// funin_num_function_call__ value
	nbFnCalls_ = 0;
	n_ = 0;
	iSteps_ = 0;
	// traceMat for results
	traceMat_ = (double*) malloc(TRACEMATSIZE * sizeof(double));
}

GenSACaller::~GenSACaller()
{
}

void GenSACaller::release()
{
	if (xMini_)
	{
		free( xMini_);
		xMini_ = 0;
	}
	if (traceMat_)
	{
		free( traceMat_);
		traceMat_ = 0;
	}
}

SEXP GenSACaller::getListElement(SEXP list, char* elementName)
{
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;

	for (i = 0; i < length(list); i++)
		if (strcmp(CHAR(STRING_ELT(names, i)), elementName) == 0)
		{
			elmt = VECTOR_ELT(list, i);
			break;
		}
	return elmt;
}

void GenSACaller::execute(SEXP x_R, SEXP par_R, SEXP lb_R, SEXP ub_R,
		SEXP fn_R, SEXP jc_R, SEXP controls_R, SEXP genSAEnvironment)
{
	// Allocate members which will be used in C section, changed an returned to R context
	// **********************************************************************************
	xMini_ = (double*) malloc(IBOUND * sizeof(double));
	for (int i = 0; i < IBOUND; ++i)
		xMini_[i] = 0;

	// lmarkov
	markovLength_ = (long int) (asInteger(getListElement(controls_R,
			(char*) "markov.length")));

	// Number of maximum step to iterate
	iSteps_ = (long int) (asInteger(getListElement(controls_R,
			(char*) "max.step")));

	for (int i = 0; i < TRACEMATSIZE; ++i)
		traceMat_[i] = 0;
	// **********************************************************************************

	nbFnCalls_ = 0;

	// set the size of the x vector (n)
	n_ = (long int) LENGTH(x_R);

	// emini (final result)
	eMini_ = -DBL_MAX;


	// Get the idum1 parameter (seed.init)
	long int seedInit = (long int) (asInteger(getListElement(controls_R,
			(char*) "seed.init")));

	// Get the idum parameter (seed.random)
	long int seedRandom = (long int) (asInteger(getListElement(controls_R,
			(char*) "seed.random")));

	// interval
	long int interval = (long int) (asInteger(getListElement(controls_R,
			(char*) "interval")));

	// know_real_energy (logical is an int for f2c converter)
	long int knowRealEnergy = (long int) (asInteger(getListElement(controls_R,
			(char*) "know.real.energy")));

	// stop_errorPercent_accordingToEnergy
	double stopErrorPercent = (double) (asReal(getListElement(controls_R,
			(char*) "error.real.energy")));

	// real_energy
	double realEnergy = (double) (asReal(getListElement(controls_R,
			(char*) "real.energy")));

	// funIn_have_constraint (logical is an int for f2c converter)
	long int hasJudgeFunction = (long int) (asInteger(getListElement(controls_R,
			(char*) "has.judge.function")));

	// lb: allocate memory and copy the R given values to the Fortran buffer
	double* lb = (double*) malloc(IBOUND * sizeof(double));
	for (int i = 0; i < IBOUND; ++i)
		lb[i] = 0;
	for (int i = 0; i < n_; ++i)
		lb[i] = REAL(lb_R)[i];

	// ub: allocate memory and copy the R given values to the Fortran buffer
	double* ub = (double*) malloc(IBOUND * sizeof(double));
	for (int i = 0; i < IBOUND; ++i)
		ub[i] = 0;
	for (int i = 0; i < n_; ++i)
		ub[i] = REAL(ub_R)[i];

	// funIn_nPara_x
	long int nPara = (long int) LENGTH(par_R);

	// funIn_para_x
	double* para = (double*) malloc(nmaxPara_x * sizeof(double));
	for (int i = 0; i < nmaxPara_x; ++i)
		para[i] = 0;
	for (int i = 0; i < nPara; i++)
		para[i] = REAL(par_R)[i];

	// x_initial: allocate memory and copy the R given values to the Fortran buffer
	double* x = (double*) malloc(IBOUND * sizeof(double));
	for (int i = 0; i < IBOUND; ++i)
		x[i] = 0;
	for (int i = 0; i < n_; ++i)
		x[i] = REAL(x_R)[i];

	// temsta
	double temperature = (double) (asReal(getListElement(controls_R,
			(char*) "temp.init")));

	// qv
	double visitingParam = (double) (asReal(getListElement(controls_R,
			(char*) "visiting.param")));

	// qa
	double acceptanceParam = (double) (asReal(getListElement(controls_R,
			(char*) "acceptance.param")));

	// funin_ls__
	long int componentChange = (long int) (asInteger(getListElement(controls_R,
			(char*) "component.change")));

#ifdef GENSA_DBG
	Rprintf("C++->seedInit:%i\n", seedInit);
	Rprintf("C++->seedRandom:%i\n", seedRandom);
	Rprintf("C++->Length of x:%i\n", n_);
	Rprintf("C++->maxStep:%i\n", iSteps_);
	Rprintf("C++->interval:%i\n", interval);
	Rprintf("C++->knowRealEnergy:%i\n", knowRealEnergy);
	Rprintf("C++->stopErrorPercent:%f\n", stopErrorPercent);
	Rprintf("C++->realEnergy:%f\n", realEnergy);
	Rprintf("C++->hasJudgeFunction:%i\n", hasJudgeFunction);
	Rprintf("C++->temperature:%f\n", temperature);
	Rprintf("C++->visitingParam:%f\n", visitingParam);
	Rprintf("C++->acceptanceParam:%f\n", acceptanceParam);
	Rprintf("C++->nb param:%i\n", nPara);
	Rprintf("C++->nb markov length:%i\n", markovLength_);

#endif

	// OptStruct will be the structure which contains the external R functions handlers, env where to call them and R names for vectors
	OptStruct OS;
	OS = (OptStruct) R_alloc(1, sizeof(opt_struct));
	OS->R_env = genSAEnvironment;
	OS->xNames = getAttrib(x_R, R_NamesSymbol);
	OS->paramNames = getAttrib(par_R, R_NamesSymbol);
	OS->params = par_R;
	OS->verbose = (int) (asInteger(
			getListElement(controls_R, (char*) "verbose")));

	PROTECT(OS->R_fn = lang2(fn_R, R_NilValue));

	if (!isNull(jc_R))
	{
		PROTECT(OS->R_jc = lang2(jc_R, R_NilValue));
	}
	else
	{
		PROTECT(OS->R_jc = R_NilValue);
	}
	int res = 0;

	res = maingsafun_(&seedInit, &seedRandom,
			&n_, xMini_, &eMini_, traceMat_,
			&nbFnCalls_, &iSteps_, &interval,
			&knowRealEnergy, &stopErrorPercent, &realEnergy, &hasJudgeFunction,
			lb, ub, &nPara, para, x, &componentChange,
			&temperature, &visitingParam, &acceptanceParam,
			&markovLength_, (void*) OS);

	UNPROTECT(2);
#ifdef GENSA_DBG
	Rprintf("C++->****** Call of workhorse done.********\n");
	Rprintf("C++->Freeing memory...\n");
	Rprintf("temperature: %f\n", temperature);
	Rprintf("nbFnCalls: %i\n", nbFnCalls_);
	Rprintf("Markov length: %i\n",markovLength_);
#endif

	// Free memory
	free(lb);
	free(ub);
	free(para);
	free(x);

#ifdef GENSA_DBG
	Rprintf("C++->It is Done.\n");
#endif
}

SEXP GenSACaller::getEnergy()
{
	SEXP returnValue = R_NilValue;
	double* doubleValuePtr = 0;
	Rf_protect(returnValue = allocVector(REALSXP, 1));
	doubleValuePtr = NUMERIC_POINTER(returnValue);
	doubleValuePtr[0] = eMini_;
	Rf_unprotect(1);
	return returnValue;
}

SEXP GenSACaller::getXMiniVector()
{
	SEXP returnValue = R_NilValue;
	double* doubleValuePtr = 0;
	Rf_protect(returnValue = allocVector(REALSXP, n_));
	doubleValuePtr = NUMERIC_POINTER(returnValue);
	memcpy(doubleValuePtr, xMini_, n_ * sizeof(double));
	Rf_unprotect(1);
	// Free xMini_ here
	free( xMini_);
	xMini_ = 0;
	return returnValue;
}

SEXP GenSACaller::getTraceMat()
{
	SEXP returnValue = R_NilValue;
	double* doubleValuePtr = 0;

	//	Rf_protect(returnValue = allocMatrix(REALSXP, TRACEMAT_NBDATA,
	//			markovLength_ * iSteps_));
	Rf_protect(returnValue = allocMatrix(REALSXP, TRACEMAT_NBDATA, TRACEMATSIZE
			/ TRACEMAT_NBDATA));
	doubleValuePtr = NUMERIC_POINTER(returnValue);
	//memcpy(doubleValuePtr, traceMat_, TRACEMAT_NBDATA * iSteps_ * markovLength_
	//		* sizeof(double));
	memcpy(doubleValuePtr, traceMat_, TRACEMATSIZE * sizeof(double));
	Rf_unprotect(1);
	// Free traceMat_ here
	//free( traceMat_);
	//traceMat_ = 0;
	return returnValue;
}

SEXP GenSACaller::getNbFuncCall()
{
	SEXP returnValue = R_NilValue;
	int* intValuePtr = 0;
	Rf_protect(returnValue = NEW_INTEGER(1));
	intValuePtr = INTEGER_POINTER(returnValue);
	intValuePtr[0] = nbFnCalls_;
	Rf_unprotect(1);
	return returnValue;
}

// Functions for RWrapping of C++ code
SEXP createInstance()
{
	SEXP Rptr = R_NilValue;
	GenSACaller* instancePtr = new GenSACaller();
	Rf_protect(Rptr = R_MakeExternalPtr((void*) instancePtr, R_NilValue,
			R_NilValue));
	R_RegisterCFinalizer(Rptr, (R_CFinalizer_t) releaseInstance);
	Rf_unprotect(1);
	return Rptr;
}

SEXP releaseInstance(SEXP R_instancePtr)
{
	GenSACaller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return errR(-1);
	instancePtr = (GenSACaller*) R_ExternalPtrAddr(R_instancePtr);
	if (instancePtr)
	{
		instancePtr->release();
		delete instancePtr;
		R_ClearExternalPtr(R_instancePtr);
		return errR(0);
	}
	else
	{
		return errR(-1);
	}
}

SEXP execute(SEXP x_R, SEXP par_R, SEXP lb_R, SEXP ub_R, SEXP fn_R, SEXP jc_R,
		SEXP controls_R, SEXP genSAEnvironment, SEXP R_instancePtr)
{
	GenSACaller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return R_NilValue;
	instancePtr = (GenSACaller*) R_ExternalPtrAddr(R_instancePtr);

	if (!instancePtr)
		return R_NilValue;

	instancePtr->execute(x_R, par_R, lb_R, ub_R, fn_R, jc_R, controls_R,
			genSAEnvironment);

	return errR(0);
}

SEXP getREnergy(SEXP R_instancePtr)
{
	GenSACaller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return R_NilValue;
	instancePtr = (GenSACaller*) R_ExternalPtrAddr(R_instancePtr);

	if (!instancePtr)
		return R_NilValue;
	return instancePtr->getEnergy();
}

SEXP getRXMiniVector(SEXP R_instancePtr)
{
	GenSACaller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return R_NilValue;
	instancePtr = (GenSACaller*) R_ExternalPtrAddr(R_instancePtr);

	if (!instancePtr)
		return R_NilValue;
	return instancePtr->getXMiniVector();
}

SEXP getRTraceMat(SEXP R_instancePtr)
{
	GenSACaller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return R_NilValue;
	instancePtr = (GenSACaller*) R_ExternalPtrAddr(R_instancePtr);

	if (!instancePtr)
		return R_NilValue;
	return instancePtr->getTraceMat();
}

SEXP getRNbFuncCall(SEXP R_instancePtr)
{
	GenSACaller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return R_NilValue;
	instancePtr = (GenSACaller*) R_ExternalPtrAddr(R_instancePtr);

	if (!instancePtr)
		return R_NilValue;
	return instancePtr->getNbFuncCall();
}

