#include <stdlib.h>
#include <string.h>
#include "lgs.h"

/**
 * Allocate memnory for input matrix
 */
void lgs_allocate_mem(lgs_t *LGS) {
    int i, j;

    LGS->matrix = (mpz_t**)calloc(LGS->num_rows, sizeof(mpz_t*));

    for (j = 0; j < LGS->num_rows; j++) {
        LGS->matrix[j] = (mpz_t*)calloc(LGS->num_cols, sizeof(mpz_t));
        for (i = 0; i < LGS->num_cols; i++)
           mpz_init(LGS->matrix[j][i]);
    }

    LGS->rhs = (mpz_t*)calloc(LGS->num_rows, sizeof(mpz_t));
    for (i = 0; i < LGS->num_rows; i++)
       mpz_init(LGS->rhs[i]);
}

void lgs_free_mem(lgs_t *LGS) {
    int i, j;

    for (j = 0; j < LGS->num_rows; j++) {
        for (i = 0; i < LGS->num_cols; i++)
			mpz_clear(LGS->matrix[j][i]);
        free(LGS->matrix[j]);
    }
    free(LGS->matrix);

    for (i = 0; i < LGS->num_rows; i++)
		mpz_clear(LGS->rhs[i]);
    free(LGS->rhs);

    if (LGS->upperbounds != NULL) {
        for (i = 0; i < LGS->num_cols; i++) {
            mpz_clear(LGS->upperbounds[i]);
        }
        free(LGS->upperbounds);
    }
}

/**
 * Read linear system from input file.
 * File must be open at this point.
 * @param txt [description]
 * @param LGS [description]
 */
void read_linear_system(FILE *txt, lgs_t *LGS) {
    int i, j, res;
    for (j = 0; j < LGS->num_rows; j++) {
        for (i = 0; i < LGS->num_cols; i++) {
            res = mpz_inp_str(LGS->matrix[j][i], txt, 10);
            if (res == 0) {
                incorrect_input_file();
            }
        }
        res = mpz_inp_str(LGS->rhs[j], txt, 10);
        if (res == 0) {
            incorrect_input_file();
        }
    }
}

/**
 * Read upper bounds and allocate memory for LGS->upperbounds
 * @param file_name [description]
 * @param LGS       [description]
 */
void read_upper_bounds(char *file_name, lgs_t *LGS) {
    int i;
    FILE *txt;
    char zeile[ZLENGTH];
    char *rowp;
    char detectstring[1024];

    LGS->upperbounds = NULL;
    txt = fopen(file_name, "r");
    if (txt == NULL) {
        printf("Could not open file %s !\n", file_name);
        fflush(stdout);
        exit(1);
    }

    zeile[0] = '\0';
    sprintf(detectstring, "BOUNDS");
    do {
        rowp = fgets(zeile, ZLENGTH, txt);
    } while ((rowp != NULL) && (strstr(zeile,detectstring) == NULL));

    LGS->num_boundedvars = LGS->num_cols;
    if (rowp == NULL) {
        LGS->upperbounds = NULL;
        printf("No %s \n",detectstring);
        fflush(stdout);
    } else {
        sscanf(zeile,"BOUNDS %d", &(LGS->num_boundedvars));
        if (LGS->num_boundedvars > 0) {
            fprintf(stderr, "Nr. bounded variables=%d\n", LGS->num_boundedvars);
        } else {
            LGS->num_boundedvars = 0;
        }

        LGS->upperbounds = (mpz_t*)calloc(LGS->num_cols, sizeof(mpz_t));
        for (i = 0; i < LGS->num_boundedvars; i++) {
            mpz_init(LGS->upperbounds[i]);
            mpz_inp_str(LGS->upperbounds[i], txt, 10);
        }
    }
    fclose(txt);
}

/**
 * Search for pre-selected variables
 * @param file_name [description]
 * @param LGS       [description]
 */
void read_selected_cols(char *file_name, lgs_t *LGS) {
    FILE *txt;
    char zeile[ZLENGTH];
    char *rowp;
    char detectstring[1024];
    int i, res;

    txt = fopen(file_name, "r");
    if (txt == NULL) {
        printf("Could not open file %s !\n", file_name);
        fflush(stdout);
        exit(1);
    }

    sprintf(detectstring, "SELECTEDCOLUMNS");
    do {
        rowp = fgets(zeile, ZLENGTH, txt);
    } while ((rowp != NULL) && (strstr(zeile, detectstring) == NULL));

     if (rowp != NULL) {
         fprintf(stderr, "SELECTEDCOLUMNS detected\n");
         fflush(stderr);
         res = fscanf(txt, "%d" , &(LGS->num_original_cols));
         if (res == (long)NULL || res == (long)EOF) {
             incorrect_input_file();
         }
     } else {
         LGS->num_original_cols = LGS->num_cols;
     }

     LGS->original_cols = (int*)calloc(LGS->num_original_cols, sizeof(int));

     if (rowp != NULL) {
         for (i = 0; i < LGS->num_original_cols; i++) {
             res = fscanf(txt, "%d", &(LGS->original_cols[i]));
             if (res == (long)NULL || res == (long)EOF) {
                 incorrect_input_file();
             }
         }
     } else {
         for (i = 0; i < LGS->num_original_cols; i++)
            LGS->original_cols[i] = 1;
     }
     fclose(txt);
}

/**
 * Content of input file can not be read.
 */
void incorrect_input_file() {
    fprintf(stderr,"Incomplete input file -> exit\n");
    fflush(stderr);
    exit(1);
}
