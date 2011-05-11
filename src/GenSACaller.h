// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#########################################################################################

#ifndef GENSACALLER_H
#define GENSACALLER_H

#include "GenSAFunction.h"

class GenSACaller
{
public:
	GenSACaller();
	~GenSACaller();
	void execute(SEXP x_R, SEXP par_R, SEXP lb_R, SEXP ub_R, SEXP fn_R,
			SEXP jc_R, SEXP controls_R, SEXP genSAEnvironment);
	SEXP getEnergy();
	SEXP getXMiniVector();
	SEXP getTraceMat();
	SEXP getNbFuncCall() ;
	void release();

protected:
	static SEXP getListElement(SEXP list, char* elementName);

private:
	double* xMini_;
	double* traceMat_;
	long int nbFnCalls_;
	double eMini_;
	long int n_ ; //real R size of the x vector
	long int iSteps_ ;
	long int markovLength_ ;
};

// C functions for wrapping C++ code for R

extern "C"
{
int maingsafun_(long int *idum1, long int *idum, long int *n,
		double *xmini, double *emini, double *tracemat,
		long int * funin_num_function_call__, long int *istep, long int *interval,
		long int *know_real_energy__,
		double * stop_errorpercent_accordingtoenergy__,
		double *real_energy__, long int *funin_have_constraint__,
		double *lb, double *ub, long int *funin_npara_x__,
		double *funin_para_x__, double * x_initial__,
		long int *funin_ls__, double *temsta, double *qv,
		double *qa, long int *lmarkov, void* ex) ;


SEXP createInstance();
SEXP releaseInstance(SEXP R_instancePtr);
SEXP execute(SEXP x_R, SEXP par_R, SEXP lb_R, SEXP ub_R, SEXP fn_R, SEXP jc_R,
		SEXP controls_R, SEXP genSAEnvironment, SEXP R_instancePtr);
SEXP getREnergy(SEXP R_instancePtr);
SEXP getRXMiniVector(SEXP R_instancePtr);
SEXP getRTraceMat(SEXP R_instancePtr);
SEXP getRNbFuncCall(SEXP R_instancePtr) ;

SEXP errR(int err)
{
	SEXP errorCode; // Error code
	int* errorCodePtr = 0;
	PROTECT(errorCode = NEW_INTEGER(1));
	errorCodePtr = INTEGER_POINTER(errorCode);
	errorCodePtr[0] = err;
	UNPROTECT(1);
	return (errorCode);
}

}

#endif //GENSACALLER_H
