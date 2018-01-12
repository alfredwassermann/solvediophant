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

    LGS->rhs = (mpz_t*)calloc((unsigned int)LGS->num_rows, sizeof(mpz_t));
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
    char *rowp;
    char zeile[ZLENGTH];
    int found_bounds = 0;
    int found_selcols = 0;

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

    do {
        rowp = fgets(zeile, ZLENGTH, txt);
        if (strstr(zeile, "BOUNDS") != NULL) {
            read_upper_bounds(txt, zeile, LGS);
            found_bounds = 1;
        } else if (strstr(zeile, "SELECTEDCOLUMNS") != NULL) {
            read_selected_cols(txt, LGS);
            found_selcols = 1;
        }
    } while (rowp != NULL);

    if (!found_bounds) {
        LGS->num_boundedvars = LGS->num_cols;
        LGS->upperbounds = NULL;
        fprintf(stderr, "No BOUNDS found \n");
    }

    if (!found_selcols) {
        fprintf(stderr, "No SELECTCOLUMNS found \n");
        LGS->num_original_cols = LGS->num_cols;
        LGS->original_cols = (int*)calloc(LGS->num_original_cols, sizeof(int));
        for (i = 0; i < LGS->num_original_cols; i++) {
           LGS->original_cols[i] = 1;
        }
    }
}

/**
 * Read upper bounds and allocate memory for LGS->upperbounds
 * @param file_name [description]
 * @param LGS       [description]
 */
void read_upper_bounds(FILE *txt, char *zeile, lgs_t *LGS) {
    int i;

    LGS->upperbounds = NULL;

    // LGS->num_boundedvars has to be the number of variables!
    sscanf(zeile, "BOUNDS %d", &(LGS->num_boundedvars));

    if (LGS->num_boundedvars > 0) {
        fprintf(stderr, "Num. bounded variables = %d\n", LGS->num_boundedvars);
    } else {
        LGS->num_boundedvars = 0;
    }

    LGS->upperbounds = (mpz_t*)calloc(LGS->num_cols, sizeof(mpz_t));
    for (i = 0; i < LGS->num_boundedvars; i++) {
        mpz_init(LGS->upperbounds[i]);
        mpz_inp_str(LGS->upperbounds[i], txt, 10);
    }
}

/**
 * Search for pre-selected variables
 * @param file_name [description]
 * @param LGS       [description]
 */
void read_selected_cols(FILE *txt, lgs_t *LGS) {
    int i, res;

    fprintf(stderr, "SELECTEDCOLUMNS detected\n");
    fflush(stderr);
    res = fscanf(txt, "%d" , &(LGS->num_original_cols));
    if (res == (long)NULL || res == (long)EOF) {
        incorrect_input_file();
    }

    LGS->original_cols = (int*)calloc(LGS->num_original_cols, sizeof(int));

    for (i = 0; i < LGS->num_original_cols; i++) {
        res = fscanf(txt, "%d", &(LGS->original_cols[i]));
        if (res == (long)NULL || res == (long)EOF) {
            incorrect_input_file();
        }
    }
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

    #if 0
        if (lead == cols - 1) {
            // System is not solvable, since
            // the last col, i.e. the RHS is
            // linearly independent from the other columns
            sgn = -1;
        }
    #endif
        row_swap = M[i];
        M[i] = M[r];
        M[r] = row_swap;

        if (M[r][lead] != 0) {
            inv = mul_inv(M[r][lead], p);
            for (i = lead; i < cols; i++) {
                if (M[r][i] != 0 && inv != 1) {
                    M[r][i] = (M[r][i] * inv) % p;
                }
            }
        }
        for (i = 0; i < rows; i++) {
            if (i != r) {
                f = M[i][lead];
                for (j = lead; j < cols; j++) {
                    if (M[r][j] != 0) {
                        M[i][j] = (M[i][j] - f * M[r][j]) % p;
                    }
                }
            }
        }
        lead++;
    }

end_rank:
    rnk = r;

    return sgn * rnk;
}

/**
 * Remove a column from the LGS.
 * Also remove the corresponding entry from upperbounds
 * and decrease num_cols and num_boundedvars
 */
void remove_column(lgs_t *LGS, int col_num) {
    int r, s;
    int cols = LGS->num_cols;
    int rows = LGS->num_rows;

    for (r = 0; r < rows; r++) {
        for (s = col_num + 1; s < cols; s++) {
            mpz_set(LGS->matrix[r][s - 1], LGS->matrix[r][s]);
        }
    }
    if (LGS->upperbounds != NULL) {
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

/**
 * Check if for every row the rhs is an integer multiple
 * of the gcd of the matrix entries.
 */
int check_gcd(lgs_t *LGS) {
    int i, j;
    int cols = LGS->num_cols;
    int rows = LGS->num_rows;
    mpz_t g;
    mpz_init(g);

    if (cols == 0) {
        return 1;
    }
    fprintf(stderr, ">> preprocess: find rows whose rhs can not be reached because of gcd\n");
    for (j = 0; j < rows; j++) {
        mpz_set(g, LGS->matrix[j][0]);
        for (i = 1; i < cols; i++) {
            if (mpz_sgn(LGS->matrix[j][i]) != 0) {
                mpz_gcd(g, g, LGS->matrix[j][i]);
            }
        }
        if (mpz_sgn(g) && !mpz_divisible_p(LGS->rhs[j], g)) {
            fprintf(stderr, "GCD check: contradiction in row %d\n", j);
            return 0;
        }
    }
    return 1;
}

/**
 * Check for every row, if
 *     sum_{i=0}^cols entry_i * upperbound_i >= rhs
 *     where entry_i > 0
 * Otherwise the system has no solution
 */
int check_rows(lgs_t *LGS) {
    int i, j;
    int cols = LGS->num_cols;
    int rows = LGS->num_rows;
    mpz_t g;
    mpz_init(g);

    if (cols == 0) {
        return 1;
    }
    fprintf(stderr, ">> preprocess: find rows whose entries are too small\n");
    for (j = 0; j < rows; j++) {
        mpz_set_ui(g, 0);
        if (LGS->upperbounds != NULL) {
            for (i = 0; i < cols; i++) {
                if (mpz_sgn(LGS->matrix[j][i]) > 0) {
                    mpz_addmul(g,  LGS->matrix[j][i], LGS->upperbounds[i]);
                }
            }
        } else {
            for (i = 0; i < cols; i++) {
                if (mpz_sgn(LGS->matrix[j][i]) > 0) {
                    mpz_addmul_ui(g, LGS->matrix[j][i], 1);
                }
            }
        }
        if (mpz_cmp(LGS->rhs[j], g) > 0) {
            fprintf(stderr, "RHS too large! contradiction in row %d\n", j);
            return 0;
        }
    }
    return 1;
}

int preprocess(lgs_t *LGS) {
    int i, j, k, found_a_column = 0;
    long p;
    int cols = LGS->num_cols;
    int rows = LGS->num_rows;
    int rnk1_a, rnk1_b, rnk2_a, rnk2_b;

    mpz_t sum_neg;

    fprintf(stderr, ">> preprocess: remove zero-forced variables\n");
    // Remove columns whose upper bounds on the variables are zero.
    cols = LGS->num_cols;
    if (LGS->upperbounds != NULL) {
        found_a_column = 0;
        for (i = cols - 1; i >= 0; i--) {
            if (mpz_sgn(LGS->upperbounds[i]) == 0) {
                remove_column(LGS, i);

                found_a_column++;
                if (found_a_column == 1)  {
                    fprintf(stderr, "Remove columns because upper bound = 0:\n");
                }
                fprintf(stderr, "%d ", i);
            }
        }
        if (found_a_column > 0) {
            fprintf(stderr, "\n");
        }
    }
    #if 0
    printf("SEL ");
    for (i = 0; i < LGS->num_original_cols; i++) {
        printf("%d ", LGS->original_cols[i]);
    }
    printf("\n");
    #endif

    // Remove columns whose entries are too large.
    #if 1
    fprintf(stderr, ">> preprocess: remove columns with entries that are too large\n");
    mpz_init(sum_neg);

    for (j = 0; j < rows; j++) {
        cols = LGS->num_cols;

        mpz_set_ui(sum_neg, 0);
        for (k = 0; k < cols; k++) {
            if (mpz_sgn(LGS->matrix[j][k]) < 0) {
                mpz_submul(sum_neg, LGS->matrix[j][k], LGS->upperbounds[k]);
            }
        }
        mpz_add(sum_neg, sum_neg, LGS->rhs[j]);

        for (i = cols - 1; i >= 0; i--) {
            if (mpz_cmp(LGS->matrix[j][i], sum_neg) > 0) {
                remove_column(LGS, i);

                found_a_column++;
                if (found_a_column == 1)  {
                    fprintf(stderr, "Remove columns because an entry is larger then rhs:\n");
                }
                fprintf(stderr, "%d ", i);
            }
        }
    }
    if (found_a_column > 0) {
        fprintf(stderr, "\n");
    }
    #endif

    if (!check_rows(LGS)) {
        return 0;
    }
    if (!check_gcd(LGS)) {
        return 0;
    }

    #if 1
        /*
        fprintf(stderr, ">> preprocess: test solvability over the rationals by rank computation\n");
        rnk1_a = rank(LGS, 7);
        LGS->num_cols--;
        rnk1_b = rank(LGS, 7);
        LGS->num_cols++;
        fprintf(stderr, "Ranks p=7: %d - %d\n", rnk1_b, rnk1_a);

        rnk1_a = rank(LGS, 11);
        LGS->num_cols--;
        rnk1_b = rank(LGS, 11);
        LGS->num_cols++;
        fprintf(stderr, "Ranks p=11: %d - %d\n", rnk1_b, rnk1_a);
        */

        #if 0
        p = 423847;
        rnk1_a = rank(LGS, );
        LGS->num_cols--;
        rnk1_b = rank(LGS, p);
        LGS->num_cols++;
        fprintf(stderr, "Ranks p=%ld: %d - %d\n", p, rnk1_b, rnk1_a);

        if (/*cols > rows &&*/ rnk1_a != rnk1_b) {
            fprintf(stderr, "First rank test mod p failed: no solution possible. Ranks: %d < %d\n", rnk1_b, rnk1_a);
            return 0;
        }
        #endif

        // Check rank
        p = 1073741827;
        rnk1_a = rank(LGS, p);
        LGS->num_cols--;
        rnk1_b = rank(LGS, p);
        LGS->num_cols++;
        fprintf(stderr, "Ranks p=%ld: %d - %d\n", p, rnk1_b, rnk1_a);

        if (/*cols > rows &&*/ rnk1_a != rnk1_b) {
            fprintf(stderr, "First rank test mod p failed: no solution possible. Ranks: %d < %d\n", rnk1_b, rnk1_a);
            return 0;
        }

        p = 2116084177;
        rnk2_a = rank(LGS, p);
        LGS->num_cols--;
        rnk2_b = rank(LGS, p);
        LGS->num_cols++;
        fprintf(stderr, "Ranks p=%ld: %d - %d\n", p, rnk2_b, rnk2_a);

        if (/*cols > rows &&*/ rnk2_a != rnk2_b) {
            fprintf(stderr, "Second rank test mod p failed: no solution possible. Ranks: %d < %d \n", rnk2_b, rnk2_a);
            return 0;
        }
        // LGS->rank = abs((rnk1 > rnk2) ? rnk1 : rnk2);
        LGS->rank = (rnk1_a > rnk2_a) ? rnk1_a : rnk2_a;

    #else
        LGS->rank = rows;
    #endif
    fprintf(stderr, "rank >= %d\n", LGS->rank);

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
