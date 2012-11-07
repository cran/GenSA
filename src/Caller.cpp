// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#########################################################################################
#include "Caller.h"

SEXP Caller::getEnergy()
{
	SEXP returnValue = R_NilValue;
	double* doubleValuePtr = 0;
	Rf_protect(returnValue = allocVector(REALSXP, 1));
	doubleValuePtr = NUMERIC_POINTER(returnValue);
	doubleValuePtr[0] = engine_.getEmini();
	Rf_unprotect(1);
	return returnValue;
}

SEXP Caller::getXMiniVector()
{
	SEXP returnValue = R_NilValue;
	double* doubleValuePtr = 0;
	Rf_protect(returnValue = allocVector(REALSXP, engine_.getX().size()));
	doubleValuePtr = NUMERIC_POINTER(returnValue);
	memcpy(doubleValuePtr, &(engine_.getXMini())[0],
			engine_.getXMini().size() * sizeof(double));
	Rf_unprotect(1);
	return returnValue;
}

SEXP Caller::getTraceMatSize()
{
	SEXP returnValue = R_NilValue;
	int* intPtr = 0;
	int res = 0;
	Tracer t = engine_.getTracer();
	res = t.getTracerLength();
	PROTECT(returnValue = NEW_INTEGER(1));
	intPtr = INTEGER_POINTER(returnValue);
	intPtr[0] = res;
	UNPROTECT(1);
	return returnValue;
}

SEXP Caller::getTraceMat(const char* key)
{
	SEXP returnValue = R_NilValue;
	std::string k = key;
	double* doubleValuePtr = 0;
	Tracer tracer = engine_.getTracer();
	unsigned int size = tracer.getTracerLength();
	if (0 == size)
	{
		return returnValue;
	}
	const double* vPtr = tracer.getVectorPtr(k);
	if (!vPtr)
	{
		return returnValue;
	}

	Rf_protect(returnValue = allocVector(REALSXP, size));
	doubleValuePtr = NUMERIC_POINTER(returnValue);
	memcpy(doubleValuePtr, vPtr, size * sizeof(double));
	Rf_unprotect(1);
	return returnValue;
}

SEXP Caller::getNbFuncCall()
{
	SEXP returnValue = R_NilValue;
	int* intValuePtr = 0;
	Rf_protect(returnValue = NEW_INTEGER(1));
	intValuePtr = INTEGER_POINTER(returnValue);
	intValuePtr[0] = engine_.getNbFctCall();
	Rf_unprotect(1);
	return returnValue;
}

SEXP Caller::getListElement(SEXP list, char* elementName)
{
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;

	for (i = 0; i < length(list); ++i)
		if (strcmp(CHAR(STRING_ELT(names, i)), elementName) == 0)
		{
			elmt = VECTOR_ELT(list, i);
			break;
		}
	return elmt;
}

void Caller::execute(SEXP x_R, SEXP lb_R, SEXP ub_R, SEXP fn_R, SEXP jc_R,
		SEXP controls_R, SEXP genSAEnvironment)
{
	int xSize = LENGTH(x_R);

	// Markov chain length
	engine_.setMarkovLength(
			asInteger(getListElement(controls_R, (char*) "markov.length")));

	// Number of maximum step to iterate
	engine_.setMaxStep(asInteger(getListElement(controls_R, (char*) "maxit")));

	// Get the idum parameter (seed.random)
	long int seedRandom = -100377;
	if (!isNull(getListElement(controls_R, (char*) "internal.random")))
	{
		seedRandom = (long int) (asInteger(
				getListElement(controls_R, (char*) "internal.random")));
	}
	else
		seedRandom = -100377;

	engine_.setInterval(
			asInteger(getListElement(controls_R, (char*) "REPORT")));

	// Real energy threshold
	if (!isNull(getListElement(controls_R, (char*) "threshold.stop")))
	{

		engine_.setKnowRealEnergy(true);
		engine_.setRealEnergyThreshold(
				asReal(getListElement(controls_R, (char*) "threshold.stop")));
	}
	else
	{
		engine_.setKnowRealEnergy(false);
	}

	// Tem restart for re-annealing
	if (!isNull(getListElement(controls_R, (char*) "tem.restart")))
	{

		engine_.setTemRestart(
				asReal(getListElement(controls_R, (char*) "tem.restart")));
	}
	else
	{
		engine_.setTemRestart(1.);
	}

	// Maximum time allowed for calculation
	if (!isNull(getListElement(controls_R, (char*) "max.time")))
	{
		engine_.setMaxTime(
				asReal(getListElement(controls_R, (char*) "max.time")));
	}
	else
	{
		engine_.setMaxTime(DBL_MAX);
	}

	// Number of improvement for stopping
	engine_.setNoImprovementStop(
			asInteger(
					getListElement(controls_R, (char*) "nb.stop.improvement")));

	// Method choice

	if (asLogical(getListElement(controls_R, (char*) "smooth")))
	{
		engine_.setMethod(Engine::SMOOTH);
	}
	else
	{
		engine_.setMethod(Engine::HARD);
	}

	// Maximum function calls
	engine_.setMaxFctCall(
			asInteger(getListElement(controls_R, (char*) "max.call")));

	// Is there any constraint function defined
	if (!isNull(jc_R))
	{
		engine_.setHasConstraint(true);
	}
	else
	{
		engine_.setHasConstraint(false);
	}

	// Is it a simple function or not
	if (!isNull(getListElement(controls_R, (char*) "simple.function")))
	{
		engine_.setIsSimpleFunction(
				(bool)asInteger(getListElement(controls_R, (char*) "simple.function")));
	}
	else
	{
		engine_.setMaxTime(DBL_MAX);
	}


	// Lower bounds
	engine_.setLower(xSize, REAL(lb_R));

	// Upper bounds
	engine_.setUpper(xSize, REAL(ub_R));

	// Initial x vector
	engine_.setX(xSize, REAL(x_R));

	// Initial temperature
	engine_.setInitialTemp(asReal(getListElement(controls_R, (char*) "temp")));

	// Visiting parameter
	engine_.setVisitingParam(
			asReal(getListElement(controls_R, (char*) "visiting.param")));

	// Acceptance parameter
	engine_.setAcceptanceParam(
			asReal(getListElement(controls_R, (char*) "acceptance.param")));

	// LSEnd param
	engine_.setLSEnd(
			(bool) (asInteger(getListElement(controls_R, (char*) "high.dim"))));

#ifdef GENSA_DBG
	Rprintf("C++->Length of x:%i\n", xSize);
	Rprintf("C++->maxStep:%i\n", engine_.getMaxStep());
	Rprintf("C++->interval:%i\n", engine_.getInterval());
	Rprintf("C++->knowRealEnergy:%i\n", engine_.knowRealEnergy());
	Rprintf("C++->thresholdStop:%f\n", engine_.getRealEnergyThreshold());
	Rprintf("C++->hasJudgeFunction:%i\n", engine_.hasConstraint());
	Rprintf("C++->temperature:%f\n", engine_.getInitialTemp());
	Rprintf("C++->visitingParam:%f\n", engine_.getVisitingParam());
	Rprintf("C++->acceptanceParam:%f\n", engine_.getAcceptanceParam());
	Rprintf("C++->markov length:%i\n", engine_.getMarkovLenght());
#endif

	// OptStruct will be the structure which contains the external R functions handlers, env where to call them and R names for vectors
	OptStruct OS;
	OS = (OptStruct) R_alloc(1, sizeof(opt_struct));
	OS->R_env = genSAEnvironment;
	OS->xNames = getAttrib(x_R, R_NamesSymbol);
	OS->verbose = (int) (asInteger(
			getListElement(controls_R, (char*) "verbose")));
	engine_.setREnv(OS);

	PROTECT(OS->R_fn = lang2(fn_R, R_NilValue));

	if (!isNull(jc_R))
	{
		PROTECT(OS->R_jc = lang2(jc_R, R_NilValue));
	}
	else
	{
		PROTECT(OS->R_jc = R_NilValue);
	}
	int err = 0;

	err = engine_.initialize();
	if (err)
	{
		UNPROTECT(2);
		return;
	}
	else
	{
		err = engine_.startSearch();
	}
	UNPROTECT(2);

#ifdef GENSA_DBG
	Rprintf("C++->It is Done.\n");
#endif
}

// Functions for RWrapping of C++ code
SEXP createInstance()
{
	SEXP Rptr = R_NilValue;
	Caller* instancePtr = new Caller();
	Rf_protect(
			Rptr = R_MakeExternalPtr((void*) instancePtr, R_NilValue,
					R_NilValue));
	R_RegisterCFinalizer(Rptr, (R_CFinalizer_t) releaseInstance);
	Rf_unprotect(1);
	return Rptr;
}

SEXP releaseInstance(SEXP R_instancePtr)
{
	//Rprintf("Releasing GenSA instance.\n") ;
	Caller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return errR(-1);
	instancePtr = (Caller*) R_ExternalPtrAddr(R_instancePtr);
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

SEXP execute(SEXP x_R, SEXP lb_R, SEXP ub_R, SEXP fn_R, SEXP jc_R,
		SEXP controls_R, SEXP genSAEnvironment, SEXP R_instancePtr)
{
	Caller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return R_NilValue;
	instancePtr = (Caller*) R_ExternalPtrAddr(R_instancePtr);

	if (!instancePtr)
		return R_NilValue;

	instancePtr->execute(x_R, lb_R, ub_R, fn_R, jc_R, controls_R,
			genSAEnvironment);

	return errR(0);
}

SEXP getREnergy(SEXP R_instancePtr)
{
	Caller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return R_NilValue;
	instancePtr = (Caller*) R_ExternalPtrAddr(R_instancePtr);

	if (!instancePtr)
		return R_NilValue;
	return instancePtr->getEnergy();
}

SEXP getRXMiniVector(SEXP R_instancePtr)
{
	Caller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return R_NilValue;
	instancePtr = (Caller*) R_ExternalPtrAddr(R_instancePtr);

	if (!instancePtr)
		return R_NilValue;
	return instancePtr->getXMiniVector();
}

SEXP getRTraceMatSize(SEXP R_instancePtr)
{
	Caller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return R_NilValue;
	instancePtr = (Caller*) R_ExternalPtrAddr(R_instancePtr);

	if (!instancePtr)
		return R_NilValue;
	return instancePtr->getTraceMatSize();
}

SEXP getRTraceMat(SEXP R_instancePtr, SEXP R_str)
{

	char* strPtr[1];
	Rf_protect(R_str = AS_CHARACTER(R_str));
	strPtr[0] = R_alloc(strlen(CHAR(STRING_ELT(R_str, 0))), sizeof(char));
	strcpy(strPtr[0], CHAR(STRING_ELT(R_str, 0)));
	Rf_unprotect(1);
	Caller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return R_NilValue;
	instancePtr = (Caller*) R_ExternalPtrAddr(R_instancePtr);

	if (!instancePtr)
		return R_NilValue;
	return instancePtr->getTraceMat(strPtr[0]);
}

SEXP getRNbFuncCall(SEXP R_instancePtr)
{
	Caller* instancePtr = 0;
	if (R_instancePtr == R_NilValue)
		return R_NilValue;
	instancePtr = (Caller*) R_ExternalPtrAddr(R_instancePtr);

	if (!instancePtr)
		return R_NilValue;
	return instancePtr->getNbFuncCall();
}

