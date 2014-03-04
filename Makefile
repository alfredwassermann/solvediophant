CC=gcc

CFLAGS= -O3 -Wall -ffast-math \
        -funroll-all-loops    \
        -march=native -mtune=native -msse -msse2 -msse3 -mssse3 -m64 \
        -pipe  \
        -fomit-frame-pointer -finline-functions  \
		-ftree-vectorize \
		-floop-interchange -floop-strip-mine -floop-block # gcc >=4.4
	#-ftree-vectorizer-verbose=2
	#-mveclibabi=svml 
	#-mfpmath=sse+387  # unstable
	#-ffast-math 
	#-fprofile-generate -fprofile-use\

#CFLAGS= -O3 -mcpu=i686 -march=i686 -fforce-addr -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -malign-functions=4

#CFLAGS= -g -Wall  
#ASSEMBLERLIB=./GotoBLAS2/libgoto2.a
#ASSEMBLERLIB=/usr/lib/openblas-base/libopenblas.a
#ASSEMBLERLIB=./OpenBLAS/libopenblas.a
ASSEMBLERLIB=./OpenBLAS/libopenblas_sandybridgep-r0.2.8.a
#GMPLIB=-L../gmp-4.2.1/bin/lib
#GMPINC=-I../gmp-4.2.1/bin/include
#GMPLIB=
#GMPINC=

all: solvediophant.dvi diophant.pdf solvediophant  

solvediophant.tex: solvediophant.w
	cweave solvediophant.w 

solvediophant.dvi: solvediophant.tex
	tex solvediophant.tex

solvediophant.c: solvediophant.w
	ctangle solvediophant.w 

solvediophant: solvediophant.c diophant.o diophant.h 
	$(CC) $(CFLAGS) -o solvediophant solvediophant.c diophant.o \
	$(ASSEMBLERLIB) \
	-lm -static -lgmp $(GMPLIB) $(GMPINC) -lpthread
	strip solvediophant
#	$(CC) -static $(CFLAGS) -o solvediophant solvediophant.c diophant.o \

diophant.tex: diophant.w
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






