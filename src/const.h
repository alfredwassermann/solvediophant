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

#ifndef _LLLCONST_H
#define _LLLCONST_H

#include <stdbool.h>

#define TRUE 1
// #define FALSE 0
#define IS_USED 0

#define WORDLEN_LONG 1
#define WORDLEN_MPZ 2

#define DEEPINSERT 1
#define DEEPINSERT_CONST 10

#define ITERATE 0
#define PROGBKZ 1
#define BKZ     2

#define CLASSIC_LLL -1
#define POT_LLL -2
#define DEEP_LLL -3
#define KERNEL_LLL -4
#define SS_LLL -5

#define VERBOSE 1

#define GIVENS 1
#define LASTLINESFACTOR 1024  // =2**10 Before: "1000000" "100000000"
#define EPSILON 0.00001       // 0.0001
#define LLLCONST_LOW 0.55     // 0.50    // 0.75
#define LLLCONST_MED 0.82
#define LLLCONST_HIGH 0.99    // 0.99
#define LLLCONST_HIGHER 0.999
#define ETACONST 0.505

#define TWOTAUHALF 67108864.0 // $2^{\tau/2}$

#define PRUNE_NO 0
#define PRUNE_HOERNER 1
#define PRUNE_BKZ 2
#define PRUNE_BKZ2 3

/* Global variable used in stop_program */
extern int SILENT;
extern int PRINT_REQUIRED;
extern int DUMP_REQUIRED;
extern bool HAS_AVX2;

#define ROUND(r) ceil(r-0.5)
#define SQRT sqrt
#define DOUBLE double

/**
 *Error codes for exit():
 * -- 0: normal program flow (reduction plus exhaustive enumeration)
 * -- 1: Input error or internal error
 * -- 2: Solution not possible, system not solvable over the reals. This may also come from parameter -c being too small
 * -- 3: Program has been called with parameters -? or -h
 * -- 4: Stop because of numerical problems in tricol
 * -- 8: Stopped after finding a random solution in phase one (''\% stopafter: 1'' has been set in the problem file)
 * -- 9: Stopped after the maximum number of solutions (''\% stopafter: n'' has been set in the problem file)
 * -- 10: Stopped after reaching the maximum number of loops (''\% stoploops: n'' has been set in the problem file)
 * -- 11: Stopped after SIGALRM, i.e. max time has been reached
 */
#define EXIT_ERR_INPUT 1
#define EXIT_NOT_SOLVABLE 2
#define EXIT_HELP 3
#define EXIT_ERR_NUMERIC 4
#define EXIT_RANDOM_SOLUTION 8
#define EXIT_MAX_SOLUTION 9
#define EXIT_MAX_LOOPS 10
#define EXIT_MAX_SIGALRM 11


#endif
