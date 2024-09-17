#include "DiscreteDLM.h"

// C function declarations
extern int sampleC(int *x, double *p, int len_p);
extern void PSI_FUN(double *xi, int *y, int *N, int *y_seq, int *ymax, double R[][*ymax], int *psi);

// Define the CallEntries array to map R functions to C functions
static const R_CMethodDef CEntries[] = {
  {"PSI_FUN", (DL_FUNC) &PSI_FUN, 7},
  {NULL, NULL, 0}
};

// Initialize the package with registered routines
void R_init_myPackage(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
