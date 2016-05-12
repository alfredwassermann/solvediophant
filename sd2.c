#include <signal.h>
#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/times.h>  /* For run time measurements */
#include <unistd.h>
#include "gls.h"
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
	gls_t GLS;
    lll_params_t LLL_params;
    lattice_t lattice;

    int maxruntime = 0;

    FILE *txt;
    char *inputfile_name;

    char zeile[ZLENGTH];
    char suffix[1024];

    FILE* solfile;
    char solfilename[1024];

    strcpy(solfilename,"solutions");

    LLL_params.iterate = LLL_params.bkz.beta = LLL_params.bkz.p = -1;
    mpz_init(LLL_params.scalelastlinefactor);
    mpz_init(lattice.matrix_factor);
    mpz_init(lattice.max_norm);

    mpz_set_si(LLL_params.scalelastlinefactor, -1);
    mpz_set_si(lattice.matrix_factor, -1);
    mpz_set_si(lattice.max_norm, -1);
    LLL_params.silent = 0;

    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i],"-silent") == 0) {
            LLL_params.silent = 1;
            fprintf(stderr,"No output of solutions, just counting.\n");

        } else if (strncmp(argv[i],"-iterate",8) == 0) {
            strcpy(suffix,argv[i]+8);
            LLL_params.iterate_no  = atoi(suffix);
            LLL_params.iterate = 1;

        } else if (strncmp(argv[i],"-bkz",4) == 0) {
            LLL_params.iterate = 0;

        } else if (strncmp(argv[i],"-beta",5) == 0) {
            strcpy(suffix,argv[i]+5);
            LLL_params.bkz.beta = atoi(suffix);

        } else if (strncmp(argv[i],"-p",2) == 0) {
            strcpy(suffix,argv[i]+2);
            LLL_params.bkz.p = atoi(suffix);

        } else if (strncmp(argv[i],"-time",5) == 0) {
            strcpy(suffix,argv[i]+5);
            maxruntime = atoi(suffix);

        } else if (strncmp(argv[i],"-c",2) == 0) {
            strcpy(suffix,argv[i]+2);
    #if 1
            mpz_set_str(lattice.matrix_factor, suffix, 10);  /* Regular version */
    #else
            mpz_ui_pow_ui(lattice.matrix_factor, 10, atoi(suffix)); /* Version for the NTL output */
    #endif
        } else if (strncmp(argv[i],"-maxnorm",8) == 0) {
            strcpy(suffix,argv[i]+8);
            mpz_set_str(lattice.max_norm,suffix,10);

        } else if (strncmp(argv[i],"-scalelastline",14) == 0) {
            strcpy(suffix,argv[i]+14);
            mpz_set_str(LLL_params.scalelastlinefactor,suffix,10);

        } else if (strncmp(argv[i],"-i",2) == 0) {
            strcpy(suffix,argv[i]+2);

        } else if (strncmp(argv[i],"-o",2) == 0) {
            strcpy(solfilename,argv[i]+2);

        } else if (strcmp(argv[i],"-?") == 0 || strcmp(argv[i],"-h") == 0) {
            fprintf(stderr,"\nsd2 --- multiple precision version --- \n");
            fprintf(stderr,"sd2");
            fprintf(stderr," -iterate*|(-bkz -beta* -p*) [-c*] [-maxnorm*] [-scalelastline*] [-time*] [-silent] [-o*]");
            fprintf(stderr," inputfile\n\n");

            exit(3);
        }
    }
    if (argc < 2 || strncmp(argv[argc-1], "-",1) == 0) {
        fprintf(stderr,"The last parameter on the command line\n");
        fprintf(stderr,"has to be the input file name.\n");
        exit(1);
    }
    if (LLL_params.iterate == -1) {
        fprintf(stderr,"No reduction was chosen.\n");
        fprintf(stderr,"It is set to iterate=1.\n");

        LLL_params.iterate = 1;
        LLL_params.iterate_no = 1;
    }
    if (LLL_params.iterate == 0 && (LLL_params.bkz.beta == -1 ||  LLL_params.bkz.p == -1)) {
        fprintf(stderr,"You have chosen bkz reduction. You also have to specify the parameters");
        fprintf(stderr," -beta* -p*\n");
        exit(1);
    }
    if (mpz_cmp_si(lattice.matrix_factor, 0) <= 0) {
        fprintf(stderr,"You did not supply the options -c*. ");
        fprintf(stderr,"It is set to 10000000000000.\n");
        mpz_set_str(lattice.matrix_factor, "10000000000000", 10);
    }
    if (mpz_cmp_si(norm_input, 0) <= 0) {
        fprintf(stderr,"You did not supply the options -maxnorm*. ");
        fprintf(stderr,"It is set to 1.\n");
        mpz_set_si(LLL_params.max_norm, 1);
    }
    if (mpz_cmp_si(LLL_params.scalelastlinefactor, 0) <= 0) {
        fprintf(stderr,"You did not supply the options -scalelastline*. ");
        fprintf(stderr,"It is set to 1.\n");
        mpz_set_si(LLL_params.scalelastlinefactor, 1);
    }

    inputfile_name = argv[argc-1];

    /**
     * Start alarm
     */
    if (maxruntime > 0) {
        signal(SIGALRM, stop_program);
        alarm(maxruntime);
    }
    /**
     * Allow debug output by kill -10 PID
     */
    signal(SIGUSR1, show_lattice);

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

    LLL_params.free_RHS = 0;
    LLL_params.stop_after_loops = 0;
    LLL_params.stop_after_solutions = 0;
    LLL_params.cut_after = -1;
    do {
        fgets(zeile, ZLENGTH, txt);
        if (strstr(zeile,"% stopafter")!=NULL) {
            sscanf(zeile,"%% stopafter %ld",&(LLL_params.stop_after_solutions));
        }
        if (strstr(zeile,"% stoploops")!=NULL) {
            sscanf(zeile,"%% stoploops %ld",&(LLL_params.stop_after_loops));
        }
        if (strstr(zeile,"% cutafter")!=NULL) {
            sscanf(zeile,"%% cutafter %d",&(LLL_params.cut_after));
        }
        if (strstr(zeile,"% FREERHS")!=NULL) {
            LLL_params.free_RHS = 1;
        }
    }
    while (zeile[0]=='%');

    sscanf(zeile,"%d%d%d", &(GLS.num_rows), &(GLS.num_cols), &flag);

	gls_allocate_mem(&GLS);
    read_linear_system(txt, &GLS);
    fclose(txt);
    read_upper_bounds(inputfile_name, &GLS);
    read_selected_cols(inputfile_name, &GLS);

    solfile = fopen(solfilename, "w");
    time_0 = os_ticks();

    diophant(&GLS, &LLL_params, solfile);

    time_1 = os_ticks();
    fclose(solfile);

    user_time = time_1 - time_0;
    timestring[0] = 0;
    print_delta_time(user_time, timestring);
    fprintf(stderr,"total enumeration time: %s\n", timestring);
    fflush(stdout);

    /**
     * Free memory
     */
    mpz_clear(factor_input);
    mpz_clear(norm_input);
    mpz_clear(LLL_params.scalelastlinefactor);

	gls_free_mem(&GLS);

    return 0;
}
