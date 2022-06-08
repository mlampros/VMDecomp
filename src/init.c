#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _VMDecomp_vmd_1d(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _VMDecomp_vmd_2d(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_VMDecomp_vmd_1d", (DL_FUNC) &_VMDecomp_vmd_1d, 8},
    {"_VMDecomp_vmd_2d", (DL_FUNC) &_VMDecomp_vmd_2d, 8},
    {NULL, NULL, 0}
};

void R_init_VMDecomp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
