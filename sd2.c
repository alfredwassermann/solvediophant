#include <signal.h>
#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/times.h>  /* For run time measurements */
#include <unistd.h>
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

int SILENT;
int PRINT_REQUIRED;
int DUMP_REQUIRED;

/**
 * Main program
 * @param  argc number of command line arguments
 * @param  argv command line arguments
 * @return
 * -- 0: normal program flow (reduction plus exhaustive enumeration)
 * -- 1: Input error or internal error
 * -- 2: Solution not possible, system not solvable over the reals. This may also come from parameter -c being too small
 * -- 3: Program has been called with parameters -? or -h
 * -- 8: Stopped after finding a random solution in phase one (''\% stopafter: 1'' has been set in the problem file)
 * -- 9: Stopped after the maximum number of solutions (''\% stopafter: n'' has been set in the problem file)
 * -- 10: Stopped after reaching the maximum number of loops (''\% stoploops: n'' has been set in the problem file)
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

    FILE* solfile;
    char sol_filename[1024];
    char restart_filename[1024];
    int restart = 0;

    strcpy(sol_filename, "solutions");

    /**
     * Init structs lattice and LLL_params
     */
    lattice.LLL_params.iterate = lattice.LLL_params.bkz.beta = lattice.LLL_params.bkz.p = -1;
    lattice.LLL_params.exhaustive_enum.lds = -1;
    lattice.LLL_params.exhaustive_enum.lds_k_max = 10;

    mpz_init(lattice.LLL_params.scalelastlinefactor);
    mpz_init(lattice.matrix_factor);
    mpz_init(lattice.max_norm);

    mpz_set_si(lattice.LLL_params.scalelastlinefactor, -1);
    mpz_set_si(lattice.matrix_factor, -1);
    mpz_set_si(lattice.max_norm, -1);
    lattice.LLL_params.silent = SILENT = 0;
    lattice.LLL_params.print_ntl = 0;
    PRINT_REQUIRED = 0;
    DUMP_REQUIRED = 0;

    /**
     * Read CLI parameters
     */
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i],"-silent") == 0) {
            lattice.LLL_params.silent = SILENT = 1;
            fprintf(stderr,"No output of solutions, just counting.\n");

        } else if (strncmp(argv[i],"-printntl", 9) == 0) {
            lattice.LLL_params.print_ntl = 1;

        } else if (strncmp(argv[i],"-iterate", 8) == 0) {
            strcpy(suffix, argv[i] + 8);
            lattice.LLL_params.iterate_no  = atoi(suffix);
            lattice.LLL_params.iterate = 1;

        } else if (strncmp(argv[i],"-bkz", 4) == 0) {
            lattice.LLL_params.iterate = 0;

        } else if (strncmp(argv[i],"-beta", 5) == 0) {
            strcpy(suffix, argv[i] + 5);
            lattice.LLL_params.bkz.beta = atoi(suffix);

        } else if (strncmp(argv[i],"-p", 2) == 0) {
            strcpy(suffix, argv[i] + 2);
            lattice.LLL_params.bkz.p = atof(suffix);

        } else if (strncmp(argv[i],"-lds", 4) == 0) {
            strcpy(suffix, argv[i] + 4);
            lattice.LLL_params.exhaustive_enum.lds = 1;
            if (strlen(suffix) > 0) {
                lattice.LLL_params.exhaustive_enum.lds_k_max = atoi(suffix);
            }

        } else if (strncmp(argv[i],"-time", 5) == 0) {
            strcpy(suffix, argv[i] + 5);
            maxruntime = atoi(suffix);

        } else if (strncmp(argv[i],"-c", 2) == 0) {
            strcpy(suffix, argv[i] + 2);
    #if 1
            mpz_set_str(lattice.matrix_factor, suffix, 10);  /* Regular version */
    #else
            mpz_ui_pow_ui(lattice.matrix_factor, 10, atoi(suffix)); /* Version for the NTL output */
    #endif
        } else if (strncmp(argv[i],"-maxnorm", 8) == 0) {
            strcpy(suffix, argv[i] + 8);
            mpz_set_str(lattice.max_norm,suffix,10);

        } else if (strncmp(argv[i],"-scalelastline", 14) == 0) {
            strcpy(suffix, argv[i] + 14);
            mpz_set_str(lattice.LLL_params.scalelastlinefactor, suffix, 10);

        } else if (strncmp(argv[i],"-i", 2) == 0) {
            strcpy(suffix, argv[i] + 2);

        } else if (strncmp(argv[i],"-o", 2) == 0) {
            strcpy(sol_filename, argv[i] + 2);

        } else if (strncmp(argv[i],"-restart", 8) == 0) {
            strcpy(restart_filename, argv[i] + 8);
            restart = 1;

        } else if (strcmp(argv[i],"-?") == 0 || strcmp(argv[i],"-h") == 0) {
            fprintf(stderr,"\nsd2 --- multiple precision version --- \n");
            fprintf(stderr,"sd2");
            fprintf(stderr," -iterate*|(-bkz -beta* -p*) [-c*] [-maxnorm*] [-scalelastline*] [-lds*] [-time*(sec)] [-silent] [-o*] [-restart*] [-printntl]");
            fprintf(stderr," inputfile\n");
            fprintf(stderr," Print lattice: kill -10 PID\n");
            fprintf(stderr," Dump lattice: kill -12 PID\n");

            exit(3);
        }
    }

    /**
     * Set default values
     */
    if (argc < 2 || strncmp(argv[argc-1], "-",1) == 0) {
        fprintf(stderr,"The last parameter on the command line\n");
        fprintf(stderr,"has to be the input file name.\n");
        exit(1);
    }
    if (lattice.LLL_params.iterate == -1) {
        fprintf(stderr,"No reduction was chosen.\n");
        fprintf(stderr,"It is set to iterate=1.\n");

        lattice.LLL_params.iterate = 1;
        lattice.LLL_params.iterate_no = 1;
    }
    if (lattice.LLL_params.iterate == 0 &&
            (lattice.LLL_params.bkz.beta == -1 /*|| lattice.LLL_params.bkz.p == -1 */)) {
        fprintf(stderr,"You have chosen bkz reduction. You also have to specify the parameters");
        fprintf(stderr," -beta* [-p*]\n");
        exit(1);
    }
    if (mpz_cmp_si(lattice.matrix_factor, 0) <= 0) {
        fprintf(stderr,"You did not supply the options -c*. ");
        fprintf(stderr,"It is set to 10000000000000.\n");
        mpz_set_str(lattice.matrix_factor, "10000000000000", 10);
    }
    if (mpz_cmp_si(lattice.max_norm, 0) <= 0) {
        fprintf(stderr,"You did not supply the options -maxnorm*. ");
        fprintf(stderr,"It is set to 1.\n");
        mpz_set_si(lattice.max_norm, 1);
    }
    if (mpz_cmp_si(lattice.LLL_params.scalelastlinefactor, 0) <= 0) {
        fprintf(stderr,"You did not supply the options -scalelastline*. ");
        fprintf(stderr,"It is set to 10000.\n");
        mpz_set_si(lattice.LLL_params.scalelastlinefactor, 10000);
    }

    if (lattice.LLL_params.exhaustive_enum.lds < 0) {
        fprintf(stderr,"Enumeration type is 'dfs'\n");
    } else {
        fprintf(stderr,"Enumeration type is 'lds'. ");
        fprintf(stderr,"lds_k_max = %d\n", lattice.LLL_params.exhaustive_enum.lds_k_max);
    }

    inputfile_name = argv[argc-1];

    /**
     * Start alarm for max run time
     */
    if (maxruntime > 0) {
        signal(SIGALRM, stop_program_sig);
        alarm(maxruntime);
    }
    /**
     * Allow debug output by kill -10 PID
     */
    signal(SIGUSR1, print_lattice_sig);
    signal(SIGUSR2, dump_lattice_sig);

    /**
     * Read options and system size in input file
     */
    txt = fopen(inputfile_name, "r");
    if (txt==NULL) {
        printf("Could not open file '%s'!\n", inputfile_name);
        fflush(stdout);
        exit(1);
    }
    flag = 0;

    lattice.free_RHS = 0;
    lattice.cut_after = -1;
    lattice.LLL_params.stop_after_loops = 0;
    lattice.LLL_params.stop_after_solutions = 0;
    do {
        fgets(zeile, ZLENGTH, txt);
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
     * Allocate memory and read LGS
     */
	lgs_allocate_mem(&LGS);
    read_linear_system(txt, &LGS);
    fclose(txt);
    read_upper_bounds(inputfile_name, &LGS);
    read_selected_cols(inputfile_name, &LGS);

    solfile = fopen(sol_filename, "w");
    time_0 = os_ticks();

    diophant(&LGS, &lattice, solfile, restart, restart_filename);

    time_1 = os_ticks();
    fclose(solfile);

    user_time = time_1 - time_0;
    timestring[0] = 0;
    print_delta_time(user_time, timestring);
    fprintf(stderr,"total enumeration time: %s\n", timestring);
    fflush(stdout);

	lgs_free_mem(&LGS);

    return 0;
}
