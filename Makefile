CC=gcc

CFLAGS= -O3 -Wall \
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

VIMFLAGS=-c 'set printoptions=number:y' -c 'set printfont=Courier:h9'
#ASSEMBLERLIB=/usr/lib/openblas-base/libblas.a
ASSEMBLERLIB=./OpenBLASsub/libopenblas.a
#ASSEMBLERLIB=
BLAS=USEBLAS
#BLAS=NOBLAS

#GMPLIB=-L../gmp-4.2.1/bin/lib
#GMPINC=-I../gmp-4.2.1/bin/include
GMPLIB=
GMPINC=

#all: solvediophant.dvi diophant.pdf solvediophant
all: sd2

%.o: %.c %.h
	$(CC) $(CFLAGS) -D$(BLAS) -c $< $(GMPINC)

# With BLAS
sd2: sd2.o dio2.o lgs.o lattice.o lll.o
	$(CC) $(CFLAGS) -o sd2 sd2.o dio2.o lll.o lattice.o lgs.o $(ASSEMBLERLIB) -lm -static -lgmp $(GMPLIB) $(GMPINC) -lpthread

dio2.pdf: dio2.c
	#vim -c 'set printfont=DejaVu\ Sans\ Mono:h9' -c 'set printoptions=number:y' -c 'hardcopy > dio2.ps' -c quit dio2.c
	vim $(VIMFLAGS) -c 'hardcopy > dio2.ps' -c quit dio2.c
	ps2pdf dio2.ps
	rm dio2.ps

sd2.pdf: sd2.c
	vim $(VIMFLAGS) -c 'hardcopy > sd2.ps' -c quit sd2.c
	ps2pdf sd2.ps
	rm sd2.ps

gls.pdf: lgs.c
	vim $(VIMFLAGS) -c 'hardcopy > gls.ps' -c quit lgs.c
	ps2pdf lgs.ps
	rm lgs.ps

.PHONY: clean
clean:
	rm *.o
