CC=gcc
#-march=pentium4 -malign-double 
CFLAGS= -O3 -Wall -ffast-math \
        -funroll-all-loops -fomit-frame-pointer -finline-functions  \
        -march=core2 -msse2 
        #-mfpmath=sse+387
        #
#CFLAGS= -g -Wall  
#ASSEMBLERLIB=-L./blas -lblas1
ASSEMBLERLIB=
#GMPLIB=-L../gmp-4.2.1/bin/lib
#GMPLIB=-L../gmp-4.2.1/bin/lib
#GMPINC=-I../gmp-4.2.1/bin/include
#GMPLIB=
#GMPINC=

all: solvediophant.dvi solvediophant diophant.pdf 

solvediophant.tex: solvediophant.w
	cweave solvediophant.w 

solvediophant.dvi: solvediophant.tex
	tex solvediophant.tex

solvediophant.c: solvediophant.w
	ctangle solvediophant.w 

solvediophant: solvediophant.c diophant.o diophant.h 
	$(CC) $(CFLAGS) -o solvediophant solvediophant.c diophant.o \
	-lm -static -lgmp $(GMPLIB) $(ASSEMBLERLIB)  $(GMPINC)
	strip solvediophant
#	$(CC) -static $(CFLAGS) -o solvediophant solvediophant.c diophant.o \

diophant.tex: diophant.w
	#cweave diophant.w 20erVersion.ch 
	cweave diophant.w

diophant.dvi: diophant.tex
	tex diophant.tex

diophant.ps: diophant.dvi
	dvips diophant.dvi

diophant.pdf: diophant.tex
	pdftex diophant.tex

diophant.c: diophant.w 
	#ctangle diophant.w 20erVersion.ch
	ctangle diophant.w 

diophant.o: diophant.c 
	$(CC) $(CFLAGS) -c diophant.c $(GMPINC)
#	$(CC) $(CFLAGS)-DBLAS -c diophant.c $(GMPINC)






