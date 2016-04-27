#include <stdlib.h>
#include <string.h>
#include "gls.h"

/**
 * Allocate memnory for input matrix
 */
void gls_allocate_mem(gls_t *GLS) {
    int i, j;

    GLS->matrix = (mpz_t**)calloc(GLS->num_rows, sizeof(mpz_t*));

    for (j = 0; j < GLS->num_rows; j++) {
        GLS->matrix[j] = (mpz_t*)calloc(GLS->num_cols, sizeof(mpz_t));
        for (i = 0; i < GLS->num_cols; i++)
           mpz_init(GLS->matrix[j][i]);
    }

    GLS->rhs = (mpz_t*)calloc(GLS->num_rows, sizeof(mpz_t));
    for (i = 0; i < GLS->num_rows; i++)
       mpz_init(GLS->rhs[i]);
}

void gls_free_mem(gls_t *GLS) {
    int i, j;

    for (j = 0; j < GLS->num_rows; j++) {
        for (i = 0; i < GLS->num_cols; i++)
			mpz_clear(GLS->matrix[j][i]);
        free(GLS->matrix[j]);
    }
    free(GLS->matrix);

    for (i = 0; i < GLS->num_rows; i++)
		mpz_clear(GLS->rhs[i]);
    free(GLS->rhs);

    if (GLS->upperbounds != NULL) {
        for (i = 0; i < GLS->num_cols; i++) {
            mpz_clear(GLS->upperbounds[i]);
        }
        free(GLS->upperbounds);
    }
}

/**
 * Read linear system from input file.
 * File must be open at this point.
 * @param txt [description]
 * @param GLS [description]
 */
void read_linear_system(FILE *txt, gls_t *GLS) {
    int i, j, res;
    for (j = 0; j < GLS->num_rows; j++) {
        for (i = 0; i < GLS->num_cols; i++) {
            res = mpz_inp_str(GLS->matrix[j][i], txt, 10);
            if (res == 0) {
                incorrect_input_file();
            }
        }
        res = mpz_inp_str(GLS->rhs[j], txt, 10);
        if (res == 0) {
            incorrect_input_file();
        }
    }
}

/**
 * Read upper bounds and allocate memory for GLS->upperbounds
 * @param file_name [description]
 * @param GLS       [description]
 */
void read_upper_bounds(char *file_name, gls_t *GLS) {
    int i;
    FILE *txt;
    char zeile[ZLENGTH];
    char *rowp;
    char detectstring[1024];

    GLS->upperbounds = NULL;
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

    GLS->num_boundedvars = GLS->num_cols;
    if (rowp == NULL) {
        GLS->upperbounds = NULL;
        printf("No %s \n",detectstring);
        fflush(stdout);
    } else {
        sscanf(zeile,"BOUNDS %d", &(GLS->num_boundedvars));
        if (GLS->num_boundedvars > 0) {
            fprintf(stderr, "Nr. bounded variables=%d\n", GLS->num_boundedvars);
        } else {
            GLS->num_boundedvars = 0;
        }

        GLS->upperbounds = (mpz_t*)calloc(GLS->num_cols, sizeof(mpz_t));
        for (i = 0; i < GLS->num_boundedvars; i++) {
            mpz_init(GLS->upperbounds[i]);
            mpz_inp_str(GLS->upperbounds[i], txt, 10);
        }
    }
    fclose(txt);
}

/**
 * Search for pre-selected variables
 * @param file_name [description]
 * @param GLS       [description]
 */
void read_selected_cols(char *file_name, gls_t *GLS) {
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
         res = fscanf(txt, "%d" , &(GLS->num_original_cols));
         if (res == (long)NULL || res == (long)EOF) {
             incorrect_input_file();
         }
     } else {
         GLS->num_original_cols = GLS->num_cols;
     }

     GLS->original_cols = (int*)calloc(GLS->num_original_cols, sizeof(int));

     if (rowp != NULL) {
         for (i = 0; i < GLS->num_original_cols; i++) {
             res = fscanf(txt, "%d", &(GLS->original_cols[i]));
             if (res == (long)NULL || res == (long)EOF) {
                 incorrect_input_file();
             }
         }
     } else {
         for (i = 0; i < GLS->num_original_cols; i++)
            GLS->original_cols[i] = 1;
     }
     fclose(txt);
}

/**
 * Content of inpout file can not be read.
 */
void incorrect_input_file() {
    fprintf(stderr,"Incomplete input file -> exit\n");
    fflush(stderr);
    exit(1);
}
