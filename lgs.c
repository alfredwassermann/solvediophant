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
        // LGS->num_boundedvars has to be the number of variables!
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
         for (i = 0; i < LGS->num_original_cols; i++) {
            LGS->original_cols[i] = 1;
         }
     }
     fclose(txt);
}
/**
 * Computes the gcd of two integers
 */
long gcd(long u, long v) {
    if (u < 0) u = -u;
    if (v < 0) v = -v;
    if (v) {
        while ((u %= v) && (v %= u));
    }
    return (u + v);
}

/**
 * Computes the multiplicative inverse of a modulo b.
 */
int mul_inv(long a, long b) {
	long b0 = b, t, q;
    long x0 = 0, x1 = 1;

	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}

/**
 * Returns the rank of the matrix of a linear system
 * (M | rhs) modulo a prime p
 * The return value is the negative of the rank
 * if the rhs column is linearly independent from the other columns.
 * This means, the linear system has no solution in Z_p.
 */
int rank(lgs_t *LGS, long p) {
    int i, j, r;
    int cols = LGS->num_cols;
    int rows = LGS->num_rows;
    long **M;
    int lead;
    long *row_swap;
    long inv, f;
    int rnk;
    int sgn = 1;

    M = (long**)calloc(rows, sizeof(long*));
    for (j = 0; j < rows; j++) {
        M[j] = (long*)calloc(cols + 1, sizeof(long));
        for (i = 0; i < cols; i++) {
            M[j][i] = (long)mpz_tdiv_ui(LGS->matrix[j][i], p);
        }
        M[j][cols] = (long)mpz_tdiv_ui(LGS->rhs[j], p);
    }

    cols++;
    lead = 0;
    rnk = 0;
    for (r = 0; r < rows; r++) {
        if (lead >= cols) {
            goto end_rank;
        }
        i = r;
        while (M[i][lead] == 0) {
            i++;
            if (i == rows) {
                i = r;
                lead++;
                if (lead == cols) {
                    goto end_rank;
                }
            }
        }

        if (lead == cols - 1) {
            // System is not solvable, since
            // the last col, i.e. the RHS is
            // linearly independent from the other columns
            sgn = -1;
        }
        row_swap = M[i];
        M[i] = M[r];
        M[r] = row_swap;

        if (M[r][lead] != 0) {
            inv = mul_inv(M[r][lead], p);
            for (i = lead; i < cols; i++) {
                M[r][i] = (M[r][i] * inv) % p;
            }
        }
        for (i = 0; i < rows; i++) {
            if (i != r) {
                f = M[i][lead];
                for (j = lead; j < cols; j++) {
                    M[i][j] = (M[i][j] - f * M[r][j]) % p;
                }
            }
        }
        lead++;
    }
end_rank:
    rnk = r;

    return sgn * rnk;
}

void remove_column(lgs_t *LGS, int col_num) {
    int r, s;
    int cols = LGS->num_cols;
    int rows = LGS->num_rows;

    for (r = 0; r < rows; r++) {
        for (s = col_num + 1; s < cols; s++) {
            mpz_set(LGS->matrix[r][s - 1], LGS->matrix[r][s]);
        }
    }
    if (LGS->num_boundedvars > 0) {
        for (s = col_num + 1; s < LGS->num_boundedvars; s++) {
            mpz_set(LGS->upperbounds[s - 1], LGS->upperbounds[s]);
        }
        LGS->num_boundedvars--;
    }

    r = 0;
    for (s = 0; s < LGS->num_original_cols; s++) {
        if (LGS->original_cols[s] != 0) {
            if (r == col_num) {
                LGS->original_cols[s] = 0;
                LGS->num_cols--;
                break;
            }
            r++;
        }
    }

}

int preprocess(lgs_t *LGS) {
    int i, j;
    int cols = LGS->num_cols;
    int rows = LGS->num_rows;
    int rnk1, rnk2;

    // Remove columns whose entries are too large.
    for (i = cols - 1; i >= 0; i--) {
        for (j = 0; j < rows; j++) {
            if (mpz_cmp(LGS->matrix[j][i], LGS->rhs[j]) > 0) {
                // Delete column i
                fprintf(stderr, "Remove column %d\n", i);
                remove_column(LGS, i);
                break;
            }
        }
    }

    // Remove columns whose upper bounds on the variables are zero.
    cols = LGS->num_cols;
    if (LGS->num_boundedvars > 0) {
        for (i = cols - 1; i >= 0; i--) {
            if (mpz_sgn(LGS->upperbounds[i]) == 0) {
                // Delete column i
                fprintf(stderr, "Remove column %d (upper bound = 0)\n", i);
                remove_column(LGS, i);
            }
        }
    }
    printf("SEL ");
    for (i = 0; i < LGS->num_original_cols; i++) {
        printf("%d ", LGS->original_cols[i]);
    }
    printf("\n");

    // Check rank
    rnk1 = rank(LGS, 1073741827);
    if (rnk1 <= 0) {
        fprintf(stderr, "Rank mod p: no solution possible.");
        return 0;
    }
    rnk2 = rank(LGS, 1073741827);
    if (rnk2 <= 0) {
        fprintf(stderr, "Rank mod p: no solution possible.");
        return 0;
    }
    LGS->rank = (rnk1 > rnk2) ? rnk1 : rnk2;
    fprintf(stderr, "Rank=%d\n", LGS->rank);

    return 1;
}

/**
 * Content of input file can not be read.
 */
void incorrect_input_file() {
    fprintf(stderr,"Incomplete input file -> exit\n");
    fflush(stderr);
    exit(1);
}
