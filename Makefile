CC=gcc

# CFLAGS= -O3 -Wall -g -pg

CFLAGS= -O3 -Wall \
    -funroll-all-loops    \
    --param max-unroll-times=2 \
	-flto \
    -march=native -mtune=native -msse -msse2 -msse3 -mssse3 -m64 \
    -pipe  \
	-finline-functions  \
    -fomit-frame-pointer \
	-ftree-vectorize \
	-floop-interchange -floop-strip-mine -floop-block # gcc >=4.4
# -O3 -ffast-math -fno-cx-limited-range -funroll-loops --param max-unroll-times=2 -march=native

#CFLAGS= -O3 -Wall -funroll-all-loops

	#	-fprofile-use \
	#-ffinite-math-only -fno-trapping-math -fno-signaling-nans -fno-signed-zeros \
	#-ftree-vectorizer-verbose=2
	#-mveclibabi=svml
	#-mfpmath=sse+387  # unstable
	#-ffast-math
	#-fprofile-generate -fprofile-use\

#CFLAGS= -O3 -Wall

#CFLAGS= -O3 -mcpu=i686 -march=i686 -fforce-addr -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -malign-functions=4
#CFLAGS= -g -Wall

VIMFLAGS=-c 'set printoptions=number:y,left:2pc,right:2pc' -c 'set printfont=Courier:h9'

###################################
# USE BLAS library
#
# First option: do not use BLAS at all.
# Uncomment these:
# BLAS=NOBLAS
# BLASINC=.
# BLASLIB=
# Second option: use OpenBLAS as installed in ubuntu.
# Uncomment BLAS, BLASINC and BLASLIB:
# Modern machine:
#BLAS=USE_BLAS
# Old machine, needed for old UBT compute cluster:
# Compile it at btmdxe
#BLAS=USE_BLAS_OLD
#BLASINC=/usr/lib/x86_64-linux-gnu/
#BLASLIB=-L/usr/lib/x86_64-linux-gnu/ -lopenblas -lpthread
# Third option: use OpenBLAS installed and compiled in a folder.
# Uncomment these:
BLAS=USE_BLAS_DEV
BLASINC=../OpenBLAS/
BLASLIB=-L../OpenBLAS/ -lopenblas -lpthread
###################################

###################################
# Use of the GMP long integer library
# First option: standard install in ubuntu
GMPLIB=
GMPINC=
# Second option: use a dedicated folder containing GMP
#GMPLIB=-L../gmp-4.2.1/bin/lib
#GMPINC=-I../gmp-4.2.1/bin/include
###################################

#all: solvediophant.dvi diophant.pdf solvediophant
all: sd3 tags

%.o: %.c %.h datastruct.h
	$(CC) $(CFLAGS) -D$(BLAS) -I$(BLASINC) -c $< $(GMPINC)

#sd3: sd2.o dio2.o lgs.o lattice.o lll.o bkz.o dualbkz.o
#	$(CC) $(CFLAGS) -o sd3 sd2.o dio2.o bkz.o dualbkz.o lll.o lattice.o lgs.o $(BLASLIB) -lm -static -lgmp $(GMPLIB) $(GMPINC) -lpthread
sd3: sd2.o dio2.o lgs.o lattice.o lll.o bkz.o enum.o
	$(CC) $(CFLAGS) -o sd3 sd2.o dio2.o bkz.o lll.o lattice.o lgs.o enum.o \
	-static $(BLASLIB) \
	-lgmp $(GMPLIB) $(GMPINC) -lm -lc

dio2.pdf: dio2.c
	vim $(VIMFLAGS) -c 'hardcopy > dio2.ps' -c quit dio2.c
	ps2pdf dio2.ps
	rm dio2.ps

enum.pdf: enum.c
	vim $(VIMFLAGS) -c 'hardcopy > enum.ps' -c quit enum.c
	ps2pdf enum.ps
	rm enum.ps

lattice.pdf: lattice.c
	vim $(VIMFLAGS) -c 'hardcopy > lattice.ps' -c quit lattice.c
	ps2pdf lattice.ps
	rm lattice.ps

lll.pdf: lll.c
	vim $(VIMFLAGS) -c 'hardcopy > lll.ps' -c quit lll.c
	ps2pdf lll.ps
	rm lll.ps

sd2.pdf: sd2.c
	vim $(VIMFLAGS) -c 'hardcopy > sd2.ps' -c quit sd2.c
	ps2pdf sd2.ps
	rm sd2.ps

lgs.pdf: lgs.c
	vim $(VIMFLAGS) -c 'hardcopy > lgs.ps' -c quit lgs.c
	ps2pdf lgs.ps
	rm lgs.ps

tags: $(SRCS) Makefile
	ctags *

.PHONY: clean
clean:
	rm *.o

