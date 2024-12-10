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
#include <stdlib.h>
#include <string.h>
#include <sys/times.h>  /* For run time measurements */
#include <unistd.h>
#include <gmp.h>

#include "lgs.h"
#include "const.h"
#include "datastruct.h"
#include "lattice.h"
#include "dio2.h"

/*
 * Global variables for time measurements
 */
int user_time, time_0, time_1;
char timestring[256];

int os_ticks() {
    struct tms tms_buffer;

    if (-1 == times(&tms_buffer)) {
        return(-1);
    }
    return(tms_buffer.tms_utime);
}

int os_ticks_per_second() {
    int clk_tck = 1;

    clk_tck = sysconf(_SC_CLK_TCK);
    return(clk_tck);
}

int os_ticks_to_hms_tps(int ticks, int tps, int *h, int *m, int *s) {
    int l1;

    l1 = ticks / tps;    /* |l1| is set to overall the number of seconds. */
    *s = l1 % 60;        /* number of seconds */
    l1 -= *s;
    l1 /= 60;
    *m = l1 % 60;        /* number of minutes */
    l1 -= *m;
    l1 /= 60;            /* number of hours */
    *h = l1;
    return(1);
}

int os_ticks_to_hms(int ticks, int *h, int *m, int *s) {
    os_ticks_to_hms_tps(ticks, os_ticks_per_second(), h, m, s);
    return(1);
}

void print_delta_time_tps(int l, int tps, char *str) {
    int h, m, s;

    os_ticks_to_hms_tps(l, tps, &h, &m, &s);
    sprintf(str, "%d:%02d:%02d", h, m, s);
}

void print_delta_time(int l, char *str) {
    print_delta_time_tps(l, os_ticks_per_second(), str);
}

int get_param(int argc, char *argv[], int i, char *name, char *suffix) {
    int len;
    //char suffix[1024];

    if (i >= argc - 1 || strlen(argv[i]) >= 1024) {
        suffix = "";
        return 0;
    }

    len = strlen(name);
    if (strncmp(argv[i], name, len) != 0) {
        return 0;
    }

    strcpy(suffix, argv[i] + len);

    if (strlen(suffix) == 0 && i < argc - 2 && argv[i + 1][0] != '-') {
        // We do not have negative numbers as parameters
        //|| (argv[i + 1][0] == '-' && isdigit(argv[i + 1][1])))) {
        strcpy(suffix, argv[i + 1]);
    }

    return 1;
}

int SILENT;
int PRINT_REQUIRED;
int DUMP_REQUIRED;
bool HAS_AVX2;
bool HAS_AVX512;

/**
 * Main program
 * @param  argc number of command line arguments
 * @param  argv command line arguments
 * @see const.h
 * @return
 * -- 0: normal program flow (reduction plus exhaustive enumeration)
 * -- 1: Input error or internal error
 * -- 2: Solution not possible, system not solvable over the reals. This may also come from parameter -c being too small
 * -- 3: Program has been called with parameters -? or -h
 * -- 4: Stop because of numerical problems in tricol
 * -- 8: Stopped after finding a random solution in phase one (''\% stopafter: 1'' has been set in the problem file)
 * -- 9: Stopped after the maximum number of solutions (''\% stopafter: n'' has been set in the problem file)
 * -- 10: Stopped after reaching the maximum number of loops (''\% stoploops: n'' has been set in the problem file)
 * -- 11: Stopped after SIGALRM, i.e. max time has been reached
 *
 */
int main(int argc, char *argv[]) {
    int i, flag;

    /*
        Variables
     */
	lgs_t LGS;
    lattice_t lattice;

    int maxruntime = 0;

    FILE *txt;
    char *inputfile_name;

    char zeile[ZLENGTH];
    char suffix[1024];
    char *endptr;

    FILE* solfile;
    char sol_filename[1024];
    char restart_filename[1024];
    int restart = 0;
    char *res;

    strcpy(sol_filename, "solutions");

    if (__builtin_cpu_supports("avx512f")) {
        #if USE_AVX
            fprintf(stderr, "CPU supports AVX512\n");
            HAS_AVX512 = true;
        #else
            fprintf(stderr, "AVX512 not enabled during compilation\n");
            HAS_AVX512 = false;
        #endif
    }
    if (__builtin_cpu_supports("avx2")) {
        #if USE_AVX
            fprintf(stderr, "CPU supports AVX2\n");
            HAS_AVX2 = true;
        #else
            fprintf(stderr, "AVX2 not enabled during compilation\n");
            HAS_AVX2 = false;
        #endif
    } else {
        fprintf(stderr, "CPU does not support AVX2\n");
        HAS_AVX2 = false;
    }

    /**
     * Init structs lattice and LLL_params
     */
    lattice.LLL_params.type = lattice.LLL_params.bkz.beta = lattice.LLL_params.bkz.p = -1;
    lattice.LLL_params.exhaustive_enum.lds = -1;
    lattice.LLL_params.exhaustive_enum.lds_k_max = 10;

    lattice.LLL_params.lll.delta_low = LLLCONST_LOW;
    lattice.LLL_params.lll.delta_med = LLLCONST_MED;
    lattice.LLL_params.lll.delta_high = LLLCONST_HIGH;
    lattice.LLL_params.lll.delta_higher = LLLCONST_HIGHER;

    mpz_init(lattice.LLL_params.scalelastlinefactor);
    mpz_init(lattice.matrix_factor);
    mpz_init(lattice.max_norm);

    mpz_set_si(lattice.LLL_params.scalelastlinefactor, -1);
    mpz_set_si(lattice.matrix_factor, -1);
    mpz_set_si(lattice.max_norm, -1);
    lattice.LLL_params.silent = SILENT = 0;
    lattice.LLL_params.print_ntl = 0;
    lattice.LLL_params.kernel = 0;
    PRINT_REQUIRED = 0;
    DUMP_REQUIRED = 0;

    /**
     * Read CLI parameters
     */
    for (i = 1; i < argc; i++) {
        if (get_param(argc, argv, i, "-silent", suffix) != 0) {
            lattice.LLL_params.silent = SILENT = 1;
            fprintf(stderr,"No output of solutions, just counting.\n");

        } else if (get_param(argc, argv, i, "-printntl", suffix) != 0) {
            lattice.LLL_params.print_ntl = 1;

        } else if (get_param(argc, argv, i, "-kernel", suffix) != 0) {
            lattice.LLL_params.kernel = 1;

        } else if (get_param(argc, argv, i, "-iterate", suffix) != 0) {
            lattice.LLL_params.iterate_no  = (int)strtol(suffix, &endptr, 10);
            lattice.LLL_params.type = ITERATE;

        } else if (get_param(argc, argv, i, "-pbkz", suffix) != 0) {
            //fprintf(stderr, "SUFFIX %s\n", suffix);
            lattice.LLL_params.type = PROGBKZ;

        } else if (get_param(argc, argv, i, "-bkz", suffix) != 0) {
            //fprintf(stderr, "SUFFIX %s\n", suffix);
            lattice.LLL_params.type = BKZ;

        } else if (get_param(argc, argv, i, "-beta", suffix) != 0) {
            lattice.LLL_params.bkz.beta = (int)strtol(suffix, &endptr, 10);

        // -p is obsolete
        // } else if (get_param(argc, argv, i, "-p", suffix) != 0) {
        //     strcpy(suffix, argv[i] + 2);
        //     lattice.LLL_params.bkz.p = strtod(suffix, &endptr);

        } else if (get_param(argc, argv, i, "-lds", suffix) != 0) {
            lattice.LLL_params.exhaustive_enum.lds = 1;
            if (strlen(suffix) > 0) {
                lattice.LLL_params.exhaustive_enum.lds_k_max = (int)strtol(suffix, &endptr, 10);
            }

        } else if (get_param(argc, argv, i, "-time", suffix) != 0) {
            maxruntime = (int)strtol(suffix, &endptr, 10);

        } else if (get_param(argc, argv, i, "-c", suffix) != 0) {
            #if 1
                mpz_set_str(lattice.matrix_factor, suffix, 10);  /* Regular version */
            #else
                mpz_ui_pow_ui(lattice.matrix_factor, 10, strtoul(suffix, &endptr, 10)); /* Version for the NTL output */
            #endif

        } else if (get_param(argc, argv, i, "-C", suffix) != 0) {
            mpz_ui_pow_ui(lattice.matrix_factor, 2, strtoul(suffix, &endptr, 10));

        } else if (get_param(argc, argv, i, "-maxnorm", suffix) != 0) {
            mpz_set_str(lattice.max_norm,suffix,10);

        } else if (get_param(argc, argv, i, "-scalelastline", suffix) != 0) {
            mpz_set_str(lattice.LLL_params.scalelastlinefactor, suffix, 10);

        } else if (get_param(argc, argv, i, "-o", suffix) != 0) {
            strcpy(sol_filename, suffix);

        } else if (get_param(argc, argv, i, "-delta_low", suffix) != 0) {
           lattice.LLL_params.lll.delta_low = strtod(suffix, &endptr);

        } else if (get_param(argc, argv, i, "-delta_med", suffix) != 0) {
            lattice.LLL_params.lll.delta_med = strtod(suffix, &endptr);

        } else if (get_param(argc, argv, i, "-delta_higher", suffix) != 0) { // must be before delta_high
            lattice.LLL_params.lll.delta_higher = strtod(suffix, &endptr);

        } else if (get_param(argc, argv, i, "-delta_high", suffix) != 0) {
            lattice.LLL_params.lll.delta_high = strtod(suffix, &endptr);

        } else if (get_param(argc, argv, i, "-restart", suffix) != 0) {
            strcpy(restart_filename, suffix);
            restart = 1;

        } else if (strcmp(argv[i], "-?") == 0 || strcmp(argv[i], "-h") == 0) {
            fprintf(stderr,"sd2 --- multiple precision version --- \n");
            fprintf(stderr,"Usage:\n\tsd2 options inputfile\n");
            fprintf(stderr,"Options:\n");
            fprintf(stderr,"\t inputfile: file name or '-'  for stdin\n");
            fprintf(stderr,"\t-iterate{num} do num LLL calls with delta=delta_high\n");
            fprintf(stderr,"\t-bkz -beta{num} do BKZ with blocksize num\n");
            fprintf(stderr,"\t-pbkz -beta{num} do progressive BKZ with blocksize num\n");
            fprintf(stderr,"\t-c{num} scale equations by num (default=1099511627776=2**40)\n");
            fprintf(stderr,"\t-C{num} scale equations by 2^num (default=1099511627776=2**40)\n");
            fprintf(stderr,"\t-scalelastline{num} scale last line by num (default=1024)\n");
            fprintf(stderr,"\t-maxnorm* ???? default=1\n");
            fprintf(stderr,"\t-delta_low{num} delta for first LLL reduction\n");
            fprintf(stderr,"\t-delta_med{num} delta for second reduction\n");
            fprintf(stderr,"\t-delta_high{num} delta for second reduction and for third reduction in case of -iterate\n");
            fprintf(stderr,"\t-delta_higher{num} delta for bkz reduction\n");

            fprintf(stderr,"\t-o{string} write solutions to file 'string' (default='solutions')\n");
            fprintf(stderr,"\t-time{num} stop program after num seconds\n");
            fprintf(stderr,"\t-silent do not write solutions to stdout and solution file\n");
            fprintf(stderr,"\t-printntl write (shortened) lattice after third reduction in NTL format to stdout\n");
            fprintf(stderr,"\t-restart{string} Read dumped lattice basis from file 'string' and jump to third reduction phase\n");
            fprintf(stderr,"\t-lds{num} Use LDS enumeration up to num discrepancies, otherwise use dfs (default=dfs)\n");
            fprintf(stderr,"Signals:\n");
            fprintf(stderr,"\t 10: print lattice, e.g. kill -10 PID\n");
            fprintf(stderr,"\t 12 Dump lattice to file 'dump_lattice.b', e.g. kill -12 PID\n");

            exit(EXIT_HELP);
        }
    }

    /**
     * Set default values
     */
    if (argc < 2 ||
        (strlen(argv[argc-1]) > 1 && argv[argc-1][0] == '-')) {
        fprintf(stderr,"The last parameter on the command line has to be the input file name or '-'.\n");
        exit(EXIT_ERR_INPUT);
    }
    if (lattice.LLL_params.type == -1) {
        fprintf(stderr,"No reduction was chosen.\n");
        fprintf(stderr,"It is set to iterate=1.\n");

        lattice.LLL_params.type = ITERATE;
        lattice.LLL_params.iterate_no = 1;
    }
    if (lattice.LLL_params.type > ITERATE &&
            (lattice.LLL_params.bkz.beta == -1 /*|| lattice.LLL_params.bkz.p == -1 */)) {
        fprintf(stderr,"You have chosen bkz or pbkz reduction. You also have to specify the parameters");
        fprintf(stderr," -beta* [-p*]\n");
        exit(EXIT_ERR_INPUT);
    }
    if (mpz_cmp_si(lattice.matrix_factor, 0) <= 0) {
        fprintf(stderr,"You did not supply -c* or -C*. ");
        fprintf(stderr,"It is set to 1099511627776=2^40.\n");        // 2**40
        mpz_set_str(lattice.matrix_factor, "1099511627776", 10);
    } else {
        fprintf(stderr,"Scale factor c=");
        mpz_out_str(stderr, 10, lattice.matrix_factor);
        fprintf(stderr, "\n");
    }
    if (mpz_cmp_si(lattice.max_norm, 0) <= 0) {
        fprintf(stderr,"You did not supply -maxnorm*. ");
        fprintf(stderr,"It is set to 1.\n");
        mpz_set_si(lattice.max_norm, 1);
    }
    if (mpz_cmp_si(lattice.LLL_params.scalelastlinefactor, 0) <= 0) {
        fprintf(stderr,"You did not supply -scalelastline*. ");
        fprintf(stderr,"It is set to %d.\n", LASTLINESFACTOR);
        mpz_set_si(lattice.LLL_params.scalelastlinefactor, LASTLINESFACTOR);
    }

    if (lattice.LLL_params.exhaustive_enum.lds < 0) {
        fprintf(stderr,"Enumeration type is 'dfs'\n");
    } else {
        fprintf(stderr,"Enumeration type is 'lds'. ");
        fprintf(stderr,"lds_k_max = %d\n", lattice.LLL_params.exhaustive_enum.lds_k_max);
    }

    fprintf(stderr,"LLL deltas:\n");
    fprintf(stderr,"\t delta_low   =%lf\n", lattice.LLL_params.lll.delta_low);
    fprintf(stderr,"\t delta_med   =%lf\n", lattice.LLL_params.lll.delta_med);
    fprintf(stderr,"\t delta_high  =%lf\n", lattice.LLL_params.lll.delta_high);
    fprintf(stderr,"\t delta_higher=%lf\n", lattice.LLL_params.lll.delta_higher);

    inputfile_name = argv[argc-1];

    /**
     * Start alarm for max run time
     */
    if (maxruntime > 0) {
        fprintf(stderr, "Alarm: stop program after %d seconds\n", maxruntime);
        signal(SIGALRM, stop_program_sig);
        alarm(maxruntime);
    }
    /**
     * Allow debug output by kill -10 PID
     */
    signal(SIGUSR1, print_lattice_sig);
    signal(SIGUSR2, dump_lattice_sig);

    /**
     * Open file or stdin
     */
    if (strlen(inputfile_name) == 1 && inputfile_name[0] == '-') {
        txt = stdin;
    } else {
        txt = fopen(inputfile_name, "r");
    }

    if (txt == NULL) {
        printf("Could not open file '%s'!\n", inputfile_name);
        exit(EXIT_ERR_INPUT);
    }

    /**
     * Read additional options and system size in input file
     */
    flag = 0;
    lattice.free_RHS = 0;
    lattice.cut_after = -1;
    lattice.LLL_params.stop_after_loops = 0;
    lattice.LLL_params.stop_after_solutions = 0;
    do {
        res = fgets(zeile, ZLENGTH, txt);
        if (res == NULL) {
            exit(EXIT_ERR_INPUT);
        }
        if (strstr(zeile,"% stopafter")!=NULL) {
            sscanf(zeile,"%% stopafter %ld",&(lattice.LLL_params.stop_after_solutions));
        }
        if (strstr(zeile,"% stoploops")!=NULL) {
            sscanf(zeile,"%% stoploops %ld",&(lattice.LLL_params.stop_after_loops));
        }
        if (strstr(zeile,"% cutafter")!=NULL) {
            sscanf(zeile,"%% cutafter %d",&(lattice.cut_after));
        }
        if (strstr(zeile,"% FREERHS")!=NULL) {
            lattice.free_RHS = 1;
        }
    }
    while (zeile[0]=='%');

    /**
     * Read problem size
     */
    sscanf(zeile, "%d%d%d", &(LGS.num_rows), &(LGS.num_cols), &flag);

    /**
     * Allocate memory and read problem LGS from file
     */
	alloc_lgs(&LGS);
    read_linear_system(txt, &LGS);
    fclose(txt);

    solfile = fopen(sol_filename, "w");
    time_0 = os_ticks();

    diophant(&LGS, &lattice, solfile, restart, restart_filename);

    time_1 = os_ticks();
    fclose(solfile);

    user_time = time_1 - time_0;
    timestring[0] = 0;
    print_delta_time(user_time, timestring);
    fprintf(stderr,"total enumeration time: %s\n", timestring);
    /* fflush(stdout); */

	free_lgs(&LGS);

    return 0;
}
