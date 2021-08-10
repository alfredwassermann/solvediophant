#ifndef _DIOPHANT_H
#define _DIOPHANT_H
#include <gmp.h>
#include "lgs.h"
#include "const.h"
#include "datastruct.h"
#include "lattice.h"

/* -------------------------------------------------------------------- */
extern long diophant(lgs_t *LGS, lattice_t *lattice, FILE* solfile, int restart, char *restart_filename);

/* Basic subroutines */
extern int cutlattice(lattice_t *lattice);
extern int solutiontest(lattice_t *lattice, int position);
extern int solutiontest_long(lattice_t *lattice, int position);

extern void lll(lattice_t *lattice, int s, int z, DOUBLE quality, int reduction_type);
extern DOUBLE iteratedlll(lattice_t *lattice, int s, int z, int no_iterates, DOUBLE quality, int reduction_type);
extern DOUBLE block_reduce(lattice_t *lattice, int s, int z, int block_size, DOUBLE quality, int reduction_type);

extern void print_NTL_lattice(lattice_t *lattice);

extern void print_lattice_sig(int sig);
extern void dump_lattice_sig(int sig);
#endif
