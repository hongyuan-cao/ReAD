#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ReAD_em_hmm(void *, void *, void *, void *);

// static const R_CallMethodDef CallEntries[] = {
//     {"_ReAD_em_hmm", (DL_FUNC) &_ReAD_em_hmm, 4},
//     {NULL, NULL, 0}
// };

// void R_init_ReAD(DllInfo *dll)
// {
//     R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
//     R_useDynamicSymbols(dll, FALSE);
// }
