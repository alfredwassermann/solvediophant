#ifndef _GLS_H
#define _GLS_H
#include "dio2.h"

#define ZLENGTH 16000

extern void gls_allocate_mem(gls_t *GLS);
extern void gls_free_mem(gls_t *GLS);
extern void read_upper_bounds(char *file_name, gls_t *GLS);
extern void read_selected_cols(char *file_name, gls_t *GLS);
extern void read_linear_system(FILE *txt, gls_t *GLS);
extern void incorrect_input_file();

#endif
