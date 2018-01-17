// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// Author: Sylvain Gubian, PMP SA
//
//#############################################################################
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP createInstance();
extern SEXP execute(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getREnergy(SEXP);
extern SEXP getRNbFuncCall(SEXP);
extern SEXP getRTraceMat(SEXP, SEXP);
extern SEXP getRTraceMatSize(SEXP);
extern SEXP getRXMiniVector(SEXP);
extern SEXP releaseInstance(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"createInstance",   (DL_FUNC) &createInstance,   0},
    {"execute",          (DL_FUNC) &execute,          8},
    {"getREnergy",       (DL_FUNC) &getREnergy,       1},
    {"getRNbFuncCall",   (DL_FUNC) &getRNbFuncCall,   1},
    {"getRTraceMat",     (DL_FUNC) &getRTraceMat,     2},
    {"getRTraceMatSize", (DL_FUNC) &getRTraceMatSize, 1},
    {"getRXMiniVector",  (DL_FUNC) &getRXMiniVector,  1},
    {"releaseInstance",  (DL_FUNC) &releaseInstance,  1},
    {NULL, NULL, 0}
};

void R_init_GenSA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

