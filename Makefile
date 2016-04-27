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
#ASSEMBLERLIB=./OpenBLAS/libopenblas.a
#ASSEMBLERLIB=/usr/lib/openblas-base/libblas.a
ASSEMBLERLIB=./OpenBLAS/libopenblas.a
#ASSEMBLERLIB=

#GMPLIB=-L../gmp-4.2.1/bin/lib
#GMPINC=-I../gmp-4.2.1/bin/include
#GMPLIB=
#GMPINC=

#all: solvediophant.dvi diophant.pdf solvediophant
all: sd2

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< $(GMPINC)

sd2: sd2.o dio2.o gls.o
	$(CC) $(CFLAGS) -o sd2 sd2.o dio2.o gls.o $(ASSEMBLERLIB) -lm -static -lgmp $(GMPLIB) $(GMPINC) -lpthread

clean:
	rm *.o
