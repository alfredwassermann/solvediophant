#ifndef _LLLCONST_H
#define _LLLCONST_H

#define TRUE 1
#define FALSE 0

#define DEEPINSERT 1
#define DEEPINSERT_CONST 30

#define ITERATE 0
#define PROGBKZ 1
#define BKZ     2

#define CLASSIC_LLL -1
#define POT_LLL -2

#define VERBOSE 1

#define GIVENS 1
#define LASTLINESFACTOR "1000000" /* "100000000" */
#define EPSILON 0.00001      /* 0.0001  */
#define LLLCONST_LOW  0.60 /* 0.75*/
#define LLLCONST_MED 0.82
#define LLLCONST_HIGH 0.99    /* 0.99 */
#define LLLCONST_HIGHER 0.999
#define ETACONST 0.505

#define TWOTAUHALF 67108864.0 /* $2^{\tau/2}$*/

#define PRUNE_NO 0
#define PRUNE_HOERNER 1
#define PRUNE_BKZ 2
#define PRUNE_BKZ2 3

/* Global variable used in stop_program */
extern int SILENT;
extern int PRINT_REQUIRED;
extern int DUMP_REQUIRED;

#define ROUND(r) ceil(r-0.5)
#define SQRT sqrt
#define DOUBLE double

#endif
