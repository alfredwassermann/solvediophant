/**
Copyright 2024 Alfred Wassermann

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
“AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <signal.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <gmp.h>

#include "const.h"
#include "arith.h"
#include "lgs.h"
#include "datastruct.h"
#include "lattice.h"
#include "enum.h"

#if defined(USE_BLAS)
    #define BLAS 1
    #include <cblas-openblas.h>
#elif defined(USE_BLAS_DEV)
    #define BLAS 1
    #include "common.h"
    #include "cblas.h"
#elif defined(USE_BLAS_OLD)
    #define BLAS 1
    #include <cblas.h>
#else
    #define BLAS 0
#endif


/**
 * Exhaustive enumeration
*/
// Global variable, also used in dio2
long num_solutions;

// Defined in enum.h
solution_t solution;

/**
 * Globals for enumeration
 */
int MAX_DUAL_BOUNDS = 8192; // 8192; //1024;
int level_max;
long loops;

/*|mpz_t *upb,*lowb;|*/
long dual_bound_success;
DOUBLE dum1, dum2;

long only_zeros_no, only_zeros_success, hoelder_no, hoelder_success;
long hoelder2_success;
long cs_success;

typedef struct {
    int num;

    DOUBLE cs;
    DOUBLE us;
    DOUBLE* w;
} enum_node_t;

typedef struct {
    int pos;
    int num;
    int is_leave_count;
    int k;
    enum_node_t* nodes;
} enum_level_t;

void alloc_enum_data(enum_level_t** enum_data, DOUBLE **fipo, int cols, int rows) {
    int i, j;
    long k;
    size_t len;

    (*enum_data) = (enum_level_t*)malloc((cols + 1) * sizeof(enum_level_t));
    for (i = 0; i <= cols; i++) {
        k = 2 * ((long)(fipo[0][i]) + 1);
        k = (k > MAX_DUAL_BOUNDS) ? MAX_DUAL_BOUNDS : k;

        (*enum_data)[i].nodes = (enum_node_t*)malloc(k * sizeof(enum_node_t));

        len = DO_ALIGN(rows * sizeof(DOUBLE));
        for (j = 0; j < k; j++) {
            (*enum_data)[i].nodes[j].w = (DOUBLE*)aligned_alloc(ALIGN_SIZE, len);
            for (int l = 0; l < rows; l++) { (*enum_data)[i].nodes[j].w[l] = 0.0; }
        }
        (*enum_data)[i].k = k;
        (*enum_data)[i].num = 0;
        (*enum_data)[i].pos = 0;
        (*enum_data)[i].is_leave_count = 0;
    }
}

void free_enum_data(enum_level_t *enum_data, int cols) {
    int i, j, k;

    for (i = 0; i <= cols; i++) {
        k = enum_data[i].k;
        for (j = 0; j < k; j++) {
            free(enum_data[i].nodes[j].w);
        }
        free(enum_data[i].nodes);
    }
    free(enum_data);
}

/**
 * @brief Generate all enumeration nodes of a given level `level`.
 * The nodes are stored in `enum_data`.
 *
 * @param enum_data
 * @param lattice
 * @param us
 * @param fipo
 * @param level
 * @param max_steps
 * @return int
 *      0: after all nodes of the level have been generated
 *      -1: if LLL_params.stop_after_loops has been reached.
 */
int enumLevel(enum_level_t* enum_data, lattice_t* lattice,
                DOUBLE *us,
                DOUBLE** fipo,
                int level, int max_steps) {

    int i;
    bool goto_back;
    bool is_good;
    bool is_new_node = true;
    int columns = lattice->num_cols;
    int rows = lattice->num_rows;

    DOUBLE** bd = lattice->decomp.bd;
    DOUBLE* c = lattice->decomp.c;
    DOUBLE Fd = lattice->decomp.Fd;
    DOUBLE Fqeps = lattice->decomp.Fqeps;
    DOUBLE Fq = lattice->decomp.Fq;

    enum_level_t* ed = &(enum_data[level]);
    enum_node_t *node, *parent_node;

    DOUBLE coeff = 0.0;
    DOUBLE y = 0.0;
    DOUBLE u;
    DOUBLE u_previous = 0.0;
    long delta;
    long v;
    int eta;
    int d;
    DOUBLE norm1;

    ed->num = 0;
    node = &(ed->nodes)[ed->num];
    parent_node = &(enum_data[level + 1].nodes[enum_data[level + 1].pos]);

    delta = eta = 0;
    if (level == level_max) {
        u = 1.0;
        v = 1;
        y = 0.0;
        d = 1;
        eta = 1;
    } else {
        y = compute_y(lattice->decomp.mu_trans, us, level, level_max);
        u = ROUND(-y);
        v = (long)u;
        d = (v > -y) ? -1 : 1;
    }

    // {
    //     long lo_bound_loc, up_bound_loc;
    //     lo_bound_loc = (long)ceil(-lattice->decomp.bd_1norm[level] - y);
    //     up_bound_loc =  (long)floor(lattice->decomp.bd_1norm[level] - y);
    //     fprintf(stderr, "%d bounds on u: %ld %ld\n",
    //         level,
    //         lo_bound_loc,
    //         up_bound_loc
    //     );
    // }
    do {
        /* increase loop counter */
        loops++;
        if ((lattice->LLL_params.stop_after_loops > 0) &&
            (lattice->LLL_params.stop_after_loops <= loops)) {
            return -1;
        }

        #if VERBOSE > -1
        if (loops % 100000000 ==0) {
            fprintf(stderr, "%ld loops, solutions: %ld",
                loops, num_solutions);
            fprintf(stderr, ", dual bounds: %ld ", dual_bound_success);
            fprintf(stderr, "\n");
            fflush(stderr);
        }
        #endif

        // handle_signals(lattice, NULL);

        goto_back = false;
        is_good = true;

        /* compute new |cs| */
        coeff = u + y;
        node->cs = parent_node->cs + coeff * coeff * c[level];

        if (node->cs >= Fd)  {
            goto_back = true;
        } else if (fabs(coeff) > lattice->decomp.bd_1norm[level]) {
            /* Use (1, -1, 0, ...) as linear combination in Hoelder pruning */
            goto_back = true;
            ++hoelder2_success;
        } else if (fabs(u) > fipo[0][level] || fabs(u + 1) > fipo[1][level]) {
            dual_bound_success++;
            is_good = false;
        } else {
            if (is_new_node) {
                norm1 = compute_w(node->w, parent_node->w, bd, coeff, level, rows);
                is_new_node = false;
            } else {
                norm1 = compute_w2(node->w, bd, u - u_previous, level, rows);
            }
            u_previous = u;
            if (level > 0) {
                if (node->cs > Fqeps * norm1) {
                    ++hoelder_success;
                    is_good = false;
                    if (eta == 1) {
                        goto_back = true;
                    } else {
                        eta = 1;
                        delta *= -1;
                        if (delta * d >= 0) delta += d;
                        u = v + delta;
                        continue;
                    }
                } else {
                    // fprintf(stderr, "level=%d, cs=%0.5lf, n1=%0.5lf, diff=%0.5lf\n", level, node->cs, norm1, norm1-node->cs);
                    i = prune_only_zeros(lattice, node->w, parent_node->w,
                        level, rows,
                        Fq,
                        bd, y, columns);
                    // if (i < 0) {
                    //     goto_back = true;
                    // } else
                    if (i > 0) {
                        is_good = false;
                    }
                }
            }
        }

        if (goto_back) {
            return 0;
        } else if (is_good) {
            node->us = u;
            ed->num++;
            node = &(ed->nodes)[ed->num];
            is_new_node = true;

            if (max_steps >= 0 && ed->num >= max_steps) {
                return 0;
            }

            if (ed->num >= MAX_DUAL_BOUNDS) {
                fprintf(stderr, "%d vs. %d\n", ed->num, MAX_DUAL_BOUNDS);
                fprintf(stderr, "enum_data too small! Exit\n");
                fflush(stderr);
                exit(EXIT_ERR_INPUT);
            }
        }

        /*
            Side step: the next value in the same level is
            chosen.
        */
        if (eta == 0) {
            delta *= -1;
            if (delta * d >= 0) {
                delta += d;
            }
        } else {
            delta += d * ((delta * d >= 0) ? 1: -1);
        }
        u = v + delta;
    } while (true);

    return 0;
}

int dfs(enum_level_t* enum_data, lattice_t* lattice,
        DOUBLE *us,
        DOUBLE** fipo,
        int level) {

    int pos;
    enum_level_t* ed = &(enum_data[level]);

    if (-1 == enumLevel(enum_data, lattice, us, fipo, level, -1)) {
        return -1;
    }

    for (pos = 0; pos < ed->num; pos++) {
        us[level] = ed->nodes[pos].us;
        ed->pos = pos;

        if (level == 0) {
            // Solution found
            if (final_test(ed->nodes[pos].w, lattice->num_rows, lattice->decomp.Fq, us, lattice) == 1) {
              #if TRUE
                print_solution(lattice, ed->nodes[pos].w, lattice->num_rows, lattice->decomp.Fq, us, lattice->num_cols);
              #else
                int f = 0;
                for (int j = lattice->num_cols - 1 ; j >= 0; j--) {
                    fprintf(stderr, "%d: %d of %d\t%0.0lf\t%d",
                        j,
                        enum_data[j].pos, enum_data[j].num,
                        us[j],
                        enum_data[j].is_leave_count
                    );
                    if (enum_data[j].pos > 1) {
                        // fprintf(stderr, "%d,",enum_data[j].pos);
                        fprintf(stderr, "\t*");
                        f = 1;
                    }
                    fprintf(stderr, "\n");
                    if (j == 0) {
                        fprintf(stderr, "-------------------------------------\n");
                    }
                }
                if (f == 1) fprintf(stderr, " ");
              #endif

                if (lattice->LLL_params.stop_after_solutions > 0 &&
                    lattice->LLL_params.stop_after_solutions <= num_solutions) {
                    return -1;
                }
            }
        } else {
            if (-1 == dfs(enum_data, lattice, us, fipo, level - 1)) {
                return -1;
            }
        }
    }

    level++;
    if (level > level_max) {
        level_max = level;
    }
    return level;
}

/**
 * @brief Recursive procedure for "limited discrepancy search"
 *
 * @param enum_data
 * @param lattice
 * @param us
 * @param fipo
 * @param level
 * @param lds_k
 * @param lds_threshold
 * @return int
 *      -1: if enumLevel() returns -1 (reached max loops.)
 *      -1: if max solutions is reached
 *      -1: if recursive call of lds() return -1
 *      1: solutions are still possible
 *      2: No solution is possible
 */
int lds(enum_level_t* enum_data, lattice_t* lattice,
                DOUBLE *us,
                DOUBLE** fipo,
                int level,
                int lds_k,
                int lds_threshold) {

    int start, end;
    int pos;
    int do_left_branch_last, p;
    int next_lds_k;
    int ret;

    // 1 if we still can reach lds_k = 0, 2 otherwise.
    // In the latter case we can stop enumeration,
    // because there are no enumerations left, but we will never
    // reach lds_k == 0
    int exhausted = 2;

    // Controls the maximum values which are searched for in
    // one level in enumLevel.
    // If lds_k == 0, less interactions in enumLevel are
    // necessary, since only the first node has to be determined.
    int max_steps = -1;

    //int max_height = 0;

    // Short cut for the parameter enum_data
    enum_level_t* ed = &(enum_data[level]);

    if (level >= lds_threshold && lds_k == 0) {
        max_steps = 1;
    }

    // if (level < lds_threshold && lds_k > 0) {
    //     return exhausted;
    // }

    if (-1 == enumLevel(enum_data, lattice,
            us,
            fipo,
            level, max_steps)) {

        // enumLevel reached max. num of solutions
        return -1;
    }

    // ed->num == 0 means that in this level there are no branches to enumerate
    if (ed->num == 0) {
        ed->is_leave_count++;

        // There are still discrepancies available
        if (lds_k > 0) {
            exhausted = 2;
        } else {
            exhausted = 1;
        }
    }
    // We are in ILDS mode

    start = 1;
    do_left_branch_last = 1;
    if (level < lds_threshold) {
        // lds_threshold is typically 0.
        // If it is larger than 0, we do conventional dfs branching below this level.
        //
        start = 0;
        end = ed->num;
        do_left_branch_last = 0;
    } else {
        // lds branching
        //
        // This works for binary search trees only:
        // depth <= k: no more left branches are possible,
        // if all discrepancies should be used.
        // Without this check, search trees for smaller
        // values of lds_k are repeated.
        // We now check this at level 0.
        // if (level - lds_threshold < lds_k) {
        //     //printf("B %d %d, %d\n", level, lds_k, ed->num);
        //     start = 1;
        //     do_left_branch_last = 0;
        // }

        if (lds_k > 0) {
            // Take all nodes per level:
            end = (lds_k < ed->num) ? lds_k + 1 : ed->num;
            // Take only two nodes per level:
            //end = (lds_k < 2) ? lds_k + 1 : 2;
        } else {
            // left-branches only
            end = 1;
        }
    }

    // printf("START lds level=%d end=%d num=%d do_left_branch_last=%d lds_k=%d\n", level, end, ed->num, do_left_branch_last, lds_k);

    for (p = start; p <= ed->num; p++) {
        if (do_left_branch_last) {
            // Right branches first:
            // if (p >= end &&
            //     !(do_left_branch_last && p == ed->num)) {
            //     continue;
            // } // Otherwise we execute the leftmost branch at the end.
            if (end <= p && p < ed->num) {
                continue;
            }
            // At the end, we execute the leftmost branch at the end.
            ed->pos = pos = p % ed->num;
        } else {
            // Not relevant since this happens only if lds_threshold > 0
            pos = p;
            if (p == end) {
                break;
            }
        }

        // printf("pos=%d\n", pos);
        us[level] = ed->nodes[pos].us;
        if (level == 0) {
            // Solution found
            if (lds_k - pos == 0 &&
                final_test(ed->nodes[pos].w, lattice->num_rows, lattice->decomp.Fq, us, lattice) == 1) {
                print_solution(lattice, ed->nodes[pos].w, lattice->num_rows, lattice->decomp.Fq, us,
                    lattice->num_cols);
                #if 0
                int f = 0;
                for (int j = lattice->num_cols - 1 ; j >= 0; j--) {
                    // fprintf(stderr, "%d: %d of %d\t%0.0lf\t%d",
                    //     j,
                    //     enum_data[j].pos, enum_data[j].num,
                    //     us[j],
                    //     enum_data[j].is_leave_count
                    // );
                    if (enum_data[j].pos > 1) {
                        fprintf(stderr, "%d,",enum_data[j].pos);
                        // fprintf(stderr, "\t*");
                        f = 1;
                    }
                    // fprintf(stderr, "\n");
                    // if (j == lds_threshold) {
                    //     fprintf(stderr, "-------------------------------------\n");
                    // }
                }
                if (f == 1) fprintf(stderr, " ");
                #endif
                if (lattice->LLL_params.stop_after_solutions > 0 &&
                    lattice->LLL_params.stop_after_solutions <= num_solutions) {
                    return -1;
                }
            }
        } else {
            next_lds_k = lds_k;
            if (level > lds_threshold) {
                // We are in ILDS mode
                if (pos == 0) {
                    // depth > k, left branch
                    next_lds_k = lds_k;
                } else if (lds_k > 0) {
                    next_lds_k = (lds_k > p) ? lds_k - pos : 0;
                }
            }

            ret = lds(enum_data, lattice,
                    us, fipo, level - 1,
                    next_lds_k, lds_threshold);

            if (ret == -1) {
                return -1;
            }
            if (ret == 1) {
                exhausted = ret;
            }
        }
    }
    if (level == lds_threshold) {
        if (lds_k - ed->num > 0) {
            exhausted = 2;
        } else {
            exhausted = 1;
        }
    }


    level++;
    // if (level >= lattice->num_cols) {
    //     return 0;
    // }

    // level_max is used in enumLevel()
    // for the computation of y.
    if (level_max < level && level < lattice->num_cols) {
        // level_max is global!!!
        level_max = level;
    }
    // if (exhausted == 2) {
    //     printf("EXHAUSTED %d\n", level);
    //     fflush(stdout);

    // }
    return exhausted;
}

void init_dualbounds(lattice_t *lattice, DOUBLE ***fipo) {
    DOUBLE **muinv;
    DOUBLE entry;
    DOUBLE norm_1, norm_2;
    DOUBLE norm_1_1, norm_1_2;

    int i, j, l;
    int cols = lattice->num_cols;
    int rows = lattice->num_rows;

    (*fipo) = (DOUBLE**)calloc(cols + 1, sizeof(DOUBLE*));
    for (i = 0; i <= cols; i++) {
        (*fipo)[i] = (DOUBLE*)calloc(cols + 1, sizeof(DOUBLE));
    }
    muinv = (DOUBLE**)calloc(cols, sizeof(DOUBLE*));
    for(i = 0; i < cols; ++i) {
        muinv[i] = (DOUBLE*)calloc(rows, sizeof(DOUBLE));
    }

    /* determine inverse of mu */
    inverse(lattice->decomp.mu, muinv, cols);

    #if VERBOSE > -1
        fprintf(stderr, "Dual bounds:\n");
        fflush(stderr);
    #endif

    /* Symmetric Fincke-Pohst */
    for (i = 0; i < cols; i++) {
        norm_1 = norm_2 = 0.0;
        norm_1_1 = norm_1_2 = 0.0;
        for (j = 0; j < rows; j++) {
            entry = 0.0;
            for (l = i; l < cols; l++) {
                entry += muinv[i][l] * lattice->decomp.bd[l][j] / lattice->decomp.c[l];
            }
            norm_2 += entry * entry;
            norm_1 += fabs(entry);

            for (l = cols - 1; l < cols; l++) {
                entry += muinv[cols - 1][l] * lattice->decomp.bd[l][j] / lattice->decomp.c[l];
            }
            norm_1_2 += entry * entry;
            norm_1_1 += fabs(entry);

        }
        norm_2 = SQRT(norm_2 * lattice->decomp.Fd);
        norm_1 =  fabs(norm_1 * lattice->decomp.Fq) * (1.0 + EPSILON);
        (*fipo)[0][i] = (norm_1 < norm_2) ? norm_1 : norm_2;
        norm_1_2 = SQRT(norm_1_2 * lattice->decomp.Fd);
        norm_1_1 =  fabs(norm_1_1 * lattice->decomp.Fq) * (1.0 + EPSILON);
        (*fipo)[1][i] = (norm_1_1 < norm_1_2) ? norm_1_1 : norm_1_2;

        #if VERBOSE > -1
            fprintf(stderr, "%0.3lf ", (*fipo)[0][i]);
        #endif
    }

    #if VERBOSE > -1
        fprintf(stderr, "\n\n");
        fflush(stderr);
    #endif

    for (i = 0; i < cols; i++) {
        free(muinv[i]);
    }
    free(muinv);
}

DOUBLE exhaustive_enumeration(lattice_t *lattice) {
    /* local variables for |explicit_enumeration() */
    /*|__attribute((aligned(16)))|*/

    int level;
    int i, j, l, k;
    int result;
    size_t len;

    enum_level_t* enum_data;
    DOUBLE *us;
    mpz_t *swap_vec;

    DOUBLE **fipo;
    /* test the size of the basis */
    fprintf(stderr, "Dimension of solution space (k): %d compared to (columns - rank): %d\n",
                lattice->num_cols, lattice->lgs_cols - lattice->lgs_rank + 1 + lattice->free_RHS);
    fflush(stderr);

    if (lattice->lgs_cols - lattice->lgs_rank + 1 + lattice->free_RHS == 0) {
        fprintf(stderr, "System not solvable.\n");
        return 0;
    }

    if (lattice->num_cols < lattice->lgs_cols - lattice->lgs_rank + 1 + lattice->free_RHS) {
        fprintf(stderr,"LLL didn't succeed in computing a basis of the kernel.\n");
        fprintf(stderr,"Please increase c0 (the first parameter)!\n");
        return 0;
    }

    lattice->decomp.bit_size = get_bit_size(lattice);

    /* Allocate the memory for enumeration */
    // Integer
    lattice->decomp.first_nonzero = (int*)calloc(lattice->num_rows, sizeof(int));
    lattice->decomp.first_nonzero_in_column = (int*)calloc(lattice->num_cols + lattice->num_rows + 1, sizeof(int));
    if (lattice->decomp.first_nonzero_in_column == NULL) {
        return(0);
    }
    lattice->decomp.firstp = (int*)calloc(lattice->num_cols + 1, sizeof(int));

    // Float
    len = DO_ALIGN((lattice->num_cols + 1) * sizeof(DOUBLE));
    lattice->decomp.bd_1norm = (DOUBLE*)aligned_alloc(ALIGN_SIZE, len);
    for (int j = 0; j < lattice->num_cols + 1; j++) { lattice->decomp.bd_1norm[j] = 0.0; }

    len = DO_ALIGN((lattice->num_cols + 1) * sizeof(DOUBLE));
    us = (DOUBLE*)aligned_alloc(ALIGN_SIZE, len);
    for (int j = 0; j < lattice->num_cols + 1; j++) { us[j] = 0.0; }

    len = DO_ALIGN((lattice->num_cols + 1) * sizeof(DOUBLE*));
    lattice->decomp.mu_trans = (DOUBLE**)aligned_alloc(ALIGN_SIZE, len);
    len = DO_ALIGN((lattice->num_cols + 1) * sizeof(DOUBLE));
    for (i = 0; i <= lattice->num_cols; i++) {
        lattice->decomp.mu_trans[i]=(DOUBLE*)aligned_alloc(ALIGN_SIZE, len);
        for (int j = 0; j < lattice->num_cols + 1; j++) { lattice->decomp.mu_trans[i][j] = 0.0; }
    }
    /* End of memory allocation*/


    /* count nonzero entries in the last rows(s) */
    if (lattice->free_RHS) {
        i=0;
        for (j = lattice->num_cols - 1; j >= 0; j--)
            if (mpz_sgn(get_entry(lattice->basis, j, lattice->num_rows - 2)) != 0)
            i++;
        fprintf(stderr, "Number of nonzero entries in the second last row: %d\n", i);
        fflush(stderr);
    }

    i = 0;
    for (j = lattice->num_cols - 1; j >= 0; j--)
        if (mpz_sgn(get_entry(lattice->basis, j, lattice->num_rows - 1)) !=0 )
        i++;
    fprintf(stderr, "Number of nonzero entries in the last row: %d\n", i);
    fprintf(stderr, "Max bit size: %d\n", lattice->decomp.bit_size);
    fflush(stderr);

    // Move basis lattice->num_cols which have a nonzero entry in the last row to the end.
    // This is mandatory for lds!!!
    if (lattice->LLL_params.exhaustive_enum.lds == 1) {
        for (j = lattice->num_cols - 1; j > 0; j--) {
            for (l = j - 1; l >= 0;  l--) {
                if (mpz_cmpabs(get_entry(lattice->basis, l, lattice->num_rows - 1),
                               get_entry(lattice->basis, j, lattice->num_rows - 1)) > 0) {
                    swap_vec = lattice->basis[l];
                    for (i = l + 1; i <= j; i++) lattice->basis[i - 1] = lattice->basis[i];
                    lattice->basis[j] = swap_vec;
                }
            }
        }
        //print_lattice(lattice, stderr);
    }

    /* set the simple pruning bounds */
    lattice->decomp.Fq = (DOUBLE)mpz_get_d(lattice->max_norm);
    lattice->decomp.Fd = (lattice->num_rows * lattice->decomp.Fq * lattice->decomp.Fq) * (1.0 + EPSILON);
    lattice->decomp.Fqeps = (1.0 + EPSILON) * lattice->decomp.Fq;        // Used in prune()
    #if VERBOSE > 0
        fprintf(stderr, "Fq: %f\n", (double)lattice->decomp.Fq);
        fprintf(stderr, "Fd: %f\n", (double)lattice->decomp.Fd);
        fflush(stderr);
    #endif

    /* orthogonalize the basis */
    #if GIVENS
        givens(lattice, lattice->num_cols, lattice->num_rows, lattice->decomp.mu, lattice->decomp.bd, lattice->decomp.c);
    #else
        gramschmidt(lattice, lattice->num_cols, lattice->num_rows, lattice->decomp.mu, lattice->decomp.bd, lattice->decomp.c);
    #endif

    /* compute $mu^\top$, the transpose of $mu$. */
    for (i = 0; i < lattice->num_cols; i++) {
        for (j = 0; j < lattice->num_cols; j++) {
            lattice->decomp.mu_trans[j][i] = lattice->decomp.mu[i][j];
        }
    }

    /* Compute 1-norm of orthogonal basis */
    for (i = 0; i <= lattice->num_cols; ++i) {
        lattice->decomp.bd_1norm[i] = 0.0;
        for (j = 0; j < lattice->num_rows; ++j) {
            lattice->decomp.bd_1norm[i] += fabs(lattice->decomp.bd[i][j]);
        }
        lattice->decomp.bd_1norm[i] *= lattice->decomp.Fqeps / lattice->decomp.c[i];
    }

    dual_bound_success = 0;
    init_dualbounds(lattice, &fipo);

    /**
     * Remove trailing unnecessary lattice->num_cols.
     *
     * Contradiction to sorting columns, see above!
     * That means, columns whose corresponding dual bounds bounds
     * are equal to 0 can be removed.
     * This is important for the Selfdual Bent Functions Problems
     */
    #if 1
    for (i = lattice->num_cols - 1; i >= 0; i--) {
        if (fipo[0][i] < 0.9) {
            fprintf(stderr, "DEL\n");
            lattice->num_cols--;
        } else {
            break;
        }
    }
    #endif

    alloc_enum_data(&enum_data, fipo, lattice->num_cols, lattice->num_rows);

    /* initialize first-nonzero arrays */
    for (l = 0; l < lattice->num_rows; l++) {
        for (i = 0; i < lattice->num_cols; i++) if (mpz_sgn(get_entry(lattice->basis, i, l)) != 0) {
            lattice->decomp.first_nonzero[l] = i;
            break;
        }
    }

    fprintf(stderr, "First non-zero entries:\n");
    j = 0;
    for (l = 0; l < lattice->num_cols; l++) {
        lattice->decomp.firstp[l] = j;
        lattice->decomp.first_nonzero_in_column[j] = 0;
        j++;
        for (i = 0; i < lattice->num_rows; i++) {
            if (lattice->decomp.first_nonzero[i] == l) {
                lattice->decomp.first_nonzero_in_column[j] = i;
                lattice->decomp.first_nonzero_in_column[lattice->decomp.firstp[l]]++;
                j++;
            }
        }
        fprintf(stderr, "%d ", lattice->decomp.first_nonzero_in_column[lattice->decomp.firstp[l]]);
    }
    fprintf(stderr, ": %d\n", lattice->num_rows);
    lattice->decomp.firstp[lattice->num_cols] = j;
    lattice->decomp.first_nonzero_in_column[j] = 0;

    /* more initialization */
    level = lattice->decomp.first_nonzero[lattice->num_rows - 1];
    level = (level >= 0) ? level : 0;

    level = lattice->decomp.first_nonzero[lattice->num_rows - 1];//columns - 1;
    level_max = level;
    // us[level] = 1;
    // v[level] = 1;
    loops = 0;

    only_zeros_no = only_zeros_success = 0;
    hoelder_no = hoelder_success = hoelder2_success = 0;
    cs_success = 0;

    fprintf(stderr, "Start enumeration at level=%d\n", level); fflush(stderr);
    /* the loop of the exhaustive enumeration */
    if (lattice->LLL_params.exhaustive_enum.lds == 1) {
        /* -------- LDS -------- */
        for (k = 0; k < lattice->LLL_params.exhaustive_enum.lds_k_max; k++) {
            level_max = level;
            fprintf(stderr, "lds_k=%d\n", k); fflush(stderr);
            // result = lds(enum_data, lattice, us, fipo, level, k, lattice->num_cols / 6);
            // result = lds(enum_data, lattice, us, fipo, level, k, 3);
            result = lds(enum_data, lattice, us, fipo, level, k, 0);

            if (result == -1) {
                fprintf(stderr, "max_solutions or max_loops reached for lds_k=%d\n\n", k);
                break;
            } else if (result == 2) {
                fprintf(stderr, "No more discrepancies possible than lds_k=%d\n\n", k);
                break;
            }
        }
    } else {
        /* -------- DFS -------- */
        while (0 <= level && level < lattice->num_cols) {
            level = dfs(enum_data, lattice, us, fipo, level);
        }
    }

    /* final output */
    fprintf(stderr, "Prune_cs: %ld\n", cs_success);
    fprintf(stderr, "Prune_only_zeros: %ld of %ld\n", only_zeros_success, only_zeros_no);
    fprintf(stderr, "Prune_hoelder: %ld of %ld\n", hoelder_success, hoelder_no);
    fprintf(stderr, "Prune_hoelder interval: %ld\n", hoelder2_success);
    fprintf(stderr, "Dual bounds: %ld\n", dual_bound_success);
    fprintf(stderr, "Loops: %ld\n", loops);

    if ((lattice->LLL_params.stop_after_solutions <= num_solutions &&
         lattice->LLL_params.stop_after_solutions > 0) ||
        (lattice->LLL_params.stop_after_loops <= loops &&
         lattice->LLL_params.stop_after_loops > 0 )) {
        fprintf(stderr, "Stopped after number of solutions: %ld\n", num_solutions);

        if (lattice->LLL_params.silent)
            print_num_solutions(num_solutions);
        if ((lattice->LLL_params.stop_after_loops <= loops &&
            lattice->LLL_params.stop_after_loops > 0)) {
            exit(EXIT_MAX_LOOPS);
        } else {
            exit(EXIT_MAX_SOLUTION);
        }
    } else {
        fprintf(stderr, "Total number of solutions: %ld\n", num_solutions);
    }
    fprintf(stderr, "\n");
    /* fflush(stdout); */
    fflush(stderr);

    /* Free allocated memory for enumeration */
    // Integer
    free(lattice->decomp.first_nonzero); 
    free(lattice->decomp.first_nonzero_in_column);
    free(lattice->decomp.firstp);

    // Float
    free(lattice->decomp.bd_1norm);
    free(us);
    for (i = 0; i <= lattice->num_cols; i++) { free(lattice->decomp.mu_trans[i]); }
    free(lattice->decomp.mu_trans);

    for (i = 0; i <= lattice->num_cols; i++) {
        free(fipo[i]);
    }
    free(fipo);
    free_enum_data(enum_data, lattice->num_cols);

    // free(cs);
    // free(bd_1norm);
    // free(y);
    // free(delta);
    // free(d);
    //
    // free(eta);
    // free(v);
    // for (l = 0; l <= lattice->num_cols; l++) free(w[l]);
    // free(w);
    // free(original_lattice->num_cols);

    // for (l = 0; l <= lattice->num_cols; l++) free(dual_basis[l]);
    // free(dual_basis);
    // free(dual_bound);
    // for (l = 0; l < lattice->num_cols; l++) free(mu_trans[l]);
    // free(mu_trans);

    return num_solutions;
}

DOUBLE compute_y(DOUBLE **mu_trans, DOUBLE *us, int level, int level_max) {
    #if BLAS
        return cblas_ddot(level_max - level, &(us[level+1]), 1, &(mu_trans[level][level+1]), 1);
    #else
        int i;
        DOUBLE sum;
        i = level_max;
        sum = 0.0;
        while (i >= level + 1) {
            sum += mu_trans[level][i]*us[i];
            i--;
        }
        return sum;
    #endif
}

DOUBLE compute_w2(DOUBLE *w, DOUBLE **bd, DOUBLE alpha, int level, int rows) {
    if (HAS_AVX2) {
        return daxpy_dasum_AVX(alpha, bd[level], w, w, rows);
    } else {
        #if BLAS
            cblas_daxpy(rows, alpha, bd[level], 1, w, 1);
            return cblas_dasum(rows, w, 1);
            // return hiprec_daxpy_dasum_AVX(alpha, bd[level], w, w, rows);
            // return daxpy_dasum_AVX(alpha, bd[level], w, w, rows);
        #else
            int i;
            register DOUBLE norm1 = 0.0;
            DOUBLE *b = &(bd[level][0]);

            for (i = rows - 1; i >= 0; --i) {
                w[i] += alpha * b[i];
                norm1 += fabs(w[i]);
            }
            return norm1;
        #endif
    }
}

DOUBLE compute_w(DOUBLE *w, DOUBLE *w1, DOUBLE **bd, DOUBLE alpha, int level, int rows) {
    if (HAS_AVX2) {
        return daxpy_dasum_AVX(alpha, bd[level], w1, w, rows);
    } else {
        #if BLAS
            cblas_dcopy(rows, w1, 1, w, 1);
            cblas_daxpy(rows, alpha, bd[level], 1, w, 1);
            return cblas_dasum(rows, w, 1);
            // return hiprec_daxpy_dasum_AVX(alpha, bd[level], w1, w, rows); // Seems to be slightly slower
            // return daxpy_dasum_AVX(alpha, bd[level], w1, w, rows);
        #else
            register int i;
            register DOUBLE norm1 = 0.0;
            DOUBLE *b = &(bd[level][0]);

            for (i = 0; i < rows; ++i) {
                w[i] = w1[i] + alpha * b[i];
                norm1 += fabs(w[i]);
            }
            return norm1;
        #endif
    }
}

void gramschmidt(lattice_t *lattice, int columns, int rows, DOUBLE **mu, DOUBLE **bd, DOUBLE *c) {
    int i, l, j;
    DOUBLE sum;

    for (i = 0; i < columns; i++) {
        for (l = 0; l < rows; l++) bd[i][l] = (DOUBLE)mpz_get_d(get_entry(lattice->basis, i, l));
        for (j = 0; j < i; j++) {
            sum = 0.0;
            for (l = 0; l < rows; l++) sum += (DOUBLE)mpz_get_d(get_entry(lattice->basis, i, l)) * bd[j][l];
            mu[i][j] = sum / c[j];
            for (l = 0; l < rows; l++) bd[i][l] -= mu[i][j]*bd[j][l];
        }

        c[i] = dot_double(bd[i], bd[i], rows);
        #if VERBOSE > 0
            fprintf(stderr, "%lf ",(double)c[i]);
        #endif
    }
    #if VERBOSE > 0
        fprintf(stderr, "\n\n");
        fflush(stderr);
    #endif
    return;
}

void givens(lattice_t *lattice, int columns, int rows, DOUBLE **mu,
            DOUBLE **bd, DOUBLE *c) {
    int i,l,j;
    int mm;
    DOUBLE d1,d2;
    DOUBLE gc,gs;
    DOUBLE t;


    /* The matrix |b| is copied to |mu|.
       |bd| is set to a $z\times z$ unity matrix.
    */
    for (i = 0; i < columns; i++) {
        for (l = 0; l < rows; l++) {
            mu[i][l] = (DOUBLE)mpz_get_d(get_entry(lattice->basis, i, l));
        }
    }

    for (i = 0; i < rows; i++) {
        for (l = 0; l < rows; l++) bd[i][l] = 0.0;
        bd[i][i] = 1.0;
    }

    for (j = 1; j < rows; j++) {    /* Givens rotation */
        mm = (j < columns) ? j : columns;
        for (i = 0; i < mm; i++) {
            if (mu[i][j] != 0.0) {
                /* Stable computation of the
                   rotation coefficients.
                */
                if (fabs(mu[i][j]) >= fabs(mu[i][i])) {
                    t = mu[i][i] / mu[i][j];
                    gs = 1.0 / SQRT(1.0 + t*t);
                    gc = gs * t;
                } else {
                    t = mu[i][j] / mu[i][i];
                    gc = 1.0 / SQRT(1.0 + t*t);
                    gs = gc * t;
                }
                /* Rotation of |mu| */
                for (l = i; l < columns; l++) {
                    d1 = mu[l][i];
                    d2 = mu[l][j];
                    mu[l][i] =  gc*d1 + gs*d2;
                    mu[l][j] = -gs*d1 + gc*d2;
                }
                /* Rotation of the matrix $Q^t$ */
                for (l = 0; l < rows; l++) {
                    d1 = bd[i][l];
                    d2 = bd[j][l];
                    bd[i][l] =  gc*d1 + gs*d2;
                    bd[j][l] = -gs*d1 + gc*d2;
                }
            }
        }
    }

    /* Finally some scaling has to be done, since $Q$ is a orthonormal matrix */
    for (i = 0; i < columns; i++) {
        c[i] = mu[i][i] * mu[i][i];
        for (j = 0; j < rows; j++) {
            bd[i][j] *= mu[i][i];
        }
        for (j = columns - 1; j >= i; j--)
            mu[j][i] /= mu[i][i];

        #if VERBOSE > -1
            fprintf(stderr, "%6.3f ",(double)c[i]);
            if (i>0 && i%15==0) fprintf(stderr, "\n");
        #endif
    }
    #if VERBOSE > -1
        fprintf(stderr, "\n\n");
        fflush(stderr);
    #endif

    return;
}

void inverse(DOUBLE **mu, DOUBLE **muinv, int columns) {
    int i, j, k;
    DOUBLE sum;

    for (j = 0; j < columns; j++)
        for (i = j; i >= 0; i--) {
            sum = 0.0;
            for (k = i + 1; k < columns; k++)
                sum += mu[k][i]*muinv[k][j];
            if (i == j)
                muinv[i][j] = 1.0 - sum;
            else
                muinv[i][j] = -sum;
        }
    return;
}

/* There are several pruning methods.*/
int final_test(DOUBLE *v, int rows, DOUBLE Fq, DOUBLE *us, lattice_t *lattice) {
    register int i;
    register int k;

    i = rows - 1;
    do {
        if (fabs(v[i]) > Fq + 0.5 + EPSILON) {
            return 0;
        }
        i--;
    } while (i>=0);

    // If the involved numbers are too big,
    // an exact test is done.
    if (lattice->decomp.bit_size < 27) {
        return 1;
    }

    for (i = 0; i < rows; i++) {
        if (!lattice->is_zero_one) {
            if (mpz_cmp_si(lattice->upperbounds[i], 0) != 0) {
                mpz_divexact(solution.upfac, lattice->upperbounds_max, lattice->upperbounds[i]);
            } else {
                mpz_set(solution.upfac, lattice->upperbounds_max);
            }
        }

        mpz_set_si(solution.u,0);
        for (k = 0; k < lattice->num_cols; k++) {
            if (ROUND(us[k]) > 0) {
                mpz_addmul_ui(solution.u, get_entry(lattice->basis, k, i), ROUND(us[k]));
            } else {
                mpz_submul_ui(solution.u, get_entry(lattice->basis, k, i), -ROUND(us[k]));
            }
        }

        mpz_sub(solution.u, solution.u, solution.s);
        mpz_divexact(solution.u, solution.u, lattice->max_norm_initial);
        mpz_divexact(solution.u, solution.u, solution.upfac);
        mpz_divexact_ui(solution.u, solution.u, lattice->denom);
        mpz_abs(solution.u, solution.u);
        if (!lattice->is_zero_one && (mpz_cmp_si(solution.u, 0) < 0 ||
            mpz_cmp(solution.u, lattice->upperbounds[i]) > 0) ) {
            //fprintf(stderr," rounding error -> this is not a solution!\n");
            return 0;
        }
    }

    return 1;
}

/* Pruning according to H\"olders inequality */
int prune(DOUBLE *w, DOUBLE cs, int rows, DOUBLE Fqeps) {
    #if BLAS
        if (cs < Fqeps * cblas_dasum(rows, w, 1)) {
            return 0;
        }
    #else
        register DOUBLE reseite;
        register int i;

        reseite = 0.0; /*|-cs/Fqeps;|*/ /* | * (1-eps) | */
        i = rows - 1;
        do {
            reseite += fabs(w[i]);
            i--;
        } while (i >= 0);
        if (cs < Fqeps * reseite) return 0;
    #endif

    return 1;
}

int prune_only_zeros(lattice_t *lattice, DOUBLE *w, DOUBLE *w1,
        int level, int rows, DOUBLE Fq,
        //int *first_nonzero_in_column, int *firstp,
        DOUBLE **bd, DOUBLE y, int columns) {

    int i;
    int f;
    // DOUBLE u1, u2;

    only_zeros_no++;
    for (i = 0; i < lattice->decomp.first_nonzero_in_column[lattice->decomp.firstp[level]]; i++) {
        f = lattice->decomp.first_nonzero_in_column[lattice->decomp.firstp[level] + 1 + i];
        // u1 = ( Fq - w1[f]) / bd[level][f] - y;
        // u2 = (-Fq - w1[f]) / bd[level][f] - y;

        if (lattice->is_zero_one) {
            // if (fabs(u1 - round(u1)) > EPSILON && fabs(u2 - round(u2)) > EPSILON) {
            //     only_zeros_success++;
            //     return -1;
            // }

            if ( fabs(fabs(w[f]) - Fq) > EPSILON ) {
                only_zeros_success++;
                return 1;
            }

        } else {  /* Not zero-one */

            /* Here we have to be very conservative */
            // if (u2 - u1 <= 1.0 + EPSILON &&
            //         fabs(w[f]) < UINT32_MAX &&
            //         fabs(w[f] - round(w[f])) > 0.001) {
            //     only_zeros_success++;
            //     return -1;
            // }

            if (fabs(w[f]) > Fq * (1 + EPSILON)) {
                only_zeros_success++;
                return 1;
            }
        }
    }
    return 0;
}

int print_solution(lattice_t *lattice, DOUBLE *w, int rows, DOUBLE Fq, DOUBLE *us, int columns) {
    int i, j, k;
    int upper;
    int end;

    /* Test again, if the vector is really a solution */
    if (fabs(fabs(w[rows-1]) - Fq) > 0.5*Fq*EPSILON)  {
        return 0;
    }
    upper = rows - 1 - lattice->free_RHS;
    if (lattice->free_RHS && fabs(w[upper]) > Fq * (1 + EPSILON)) {
        return 0;
    }

    #if IS_USED
        for (k = 0; k < columns; k++) {
            fprintf(stderr, "%2.0lf ", us[k]);
        }
        fprintf(stderr, "\n");
    #endif

    if (!lattice->LLL_params.silent) {
        mpz_set_si(solution.upfac,1);
        mpz_set_si(solution.s,0);
        for (k = 0; k < columns; k++) {
            if (ROUND(us[k])>0) {
                mpz_addmul_ui(solution.s, get_entry(lattice->basis, k, rows-1), ROUND(us[k]));
            } else {
                mpz_submul_ui(solution.s, get_entry(lattice->basis, k,rows-1), -ROUND(us[k]));
            }
        }

        i = 0;
        end = (lattice->cut_after == -1) ? lattice->num_cols_org : lattice->cut_after;
        for (j = 0; j < end; j++) {
            if (lattice->original_cols[j] == 0) {
                mpz_set_si(solution.u, 0);
            } else {
                if (!lattice->is_zero_one) {
                    if (mpz_cmp_si(lattice->upperbounds[i],0) != 0) {
                        mpz_divexact(solution.upfac, lattice->upperbounds_max, lattice->upperbounds[i]);
                    } else {
                        mpz_set(solution.upfac, lattice->upperbounds_max);
                    }
                }
                mpz_set_si(solution.u, 0);
                for (k = 0; k < columns; k++) {
                    if (ROUND(us[k])>0) {
                        mpz_addmul_ui(solution.u, get_entry(lattice->basis, k, i), ROUND(us[k]));
                    } else {
                        mpz_submul_ui(solution.u, get_entry(lattice->basis, k, i), -ROUND(us[k]));
                    }
                }
                mpz_sub(solution.u, solution.u, solution.s);
                mpz_divexact(solution.u, solution.u, lattice->max_norm_initial);
                mpz_divexact(solution.u, solution.u, solution.upfac);
                mpz_divexact_ui(solution.u, solution.u, lattice->denom);
                mpz_abs(solution.u, solution.u);

                i++;
            }
            mpz_out_str(NULL, 10, solution.u);
            mpz_out_str(solution.fp, 10, solution.u);

            /* Meanwhile, all solution vectors are written with separating blanks. */
            /*|if (!lattice->iszeroone) { }|*/
            printf(" ");
            fprintf(solution.fp, " ");
        }
        if (lattice->free_RHS) {
            mpz_set_d(solution.u, ROUND(w[i]));
            mpz_divexact(solution.u, solution.u, lattice->max_up);
            mpz_abs(solution.u, solution.u);
            printf(" L = ");
            mpz_out_str(NULL, 10, solution.u);
        }
        printf("\n");
        /*fflush(stdout);*/
        fprintf(solution.fp, "\n");
        /*fflush(solution.fp);*/
    }

    num_solutions++;
    if (num_solutions%10000==0) {
        fprintf(stderr, "%ld\n", num_solutions);
        fflush(stderr);
    }

    return 1;
}

void print_num_solutions(long num_solutions) {
    fprintf(solution.fp, "%ld solutions\n", num_solutions);
    fflush(solution.fp);
}

