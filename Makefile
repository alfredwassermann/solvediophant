# CC=clang
CC=gcc
BIN=bin
SRC=src
PDF=pdf

CFLAGS= -O3                  \
  -DUSE_AVX                  \
  -fno-math-errno            \
  -ftree-vectorize           \
  -mavx -mavx2               \
  -msse3 -mssse3 -msse4.1    \
  -march=haswell             \
  -Wall


# Profile:   -pg \
# Find memory leaks
#  -fsanitize=address         \

VIMFLAGS=-c 'set printoptions=number:y,left:2pc,right:2pc' -c 'set printfont=Courier:h8'

###################################
# USE BLAS library
#
#
# ------ First option: do not use BLAS at all.
# Uncomment these:
# BLAS=NOBLAS
# BLASINC=.
# BLASLIB=
#
# ------  Second option: use OpenBLAS as installed in ubuntu.
# Uncomment BLAS, BLASINC and BLASLIB:
# Modern machine:
#BLAS=USE_BLAS
# Old machine, needed for old UBT compute cluster:
# Compile it at btmdxe
#BLAS=USE_BLAS_OLD
#BLASINC=/usr/lib/x86_64-linux-gnu/
#BLASLIB=-L/usr/lib/x86_64-linux-gnu/ -lopenblas -lpthread
#
#  ------ Third option: use OpenBLAS installed and compiled in a folder.
# Uncomment these:
BLAS=USE_BLAS_DEV
BLASINC=../OpenBLAS/
BLASLIB=-L../OpenBLAS/ -lopenblas -lpthread
###################################

###################################
# Use of the GMP long integer library
#  ------ First option: system install in linux
GMPLIB=
GMPINC=
#  ------ Second option: use a dedicated folder containing GMP
#GMPLIB=-L../gmp-4.2.1/bin/lib
#GMPINC=-I../gmp-4.2.1/bin/include
###################################

###################################
# Print files gsa.out and gsa1.out
# Otherwise comment out
GSA_OUT=FALSE
###################################

OBJFILES=$(SRC)/bkz.o   $(SRC)/dio2.o   $(SRC)/dualbkz.o   $(SRC)/enum.o   $(SRC)/lattice.o   $(SRC)/lgs.o   $(SRC)/lll.o   $(SRC)/sd2.o   $(SRC)/arith.o
PDFFILES=$(PDF)/bkz.pdf $(PDF)/dio2.pdf $(PDF)/dualbkz.pdf $(PDF)/enum.pdf $(PDF)/lattice.pdf $(PDF)/lgs.pdf $(PDF)/lll.pdf $(PDF)/sd2.pdf $(PDF)/arith.pdf

all: $(BIN)/sd2 tags $(PDFFILES) $(BIN)/test_la

$(SRC)/%.o: $(SRC)/%.c $(SRC)/%.h $(SRC)/datastruct.h $(SRC)/const.h
	$(CC) $(CFLAGS) $(PGO_CFLAGS) -D$(GSA_OUT) -D$(BLAS) -I$(BLASINC) -c $< $(GMPINC) -o $@
	@echo "Compiled "$<" successfully!"

$(PDF)/%.pdf: $(SRC)/%.c
	vim $(VIMFLAGS) -c 'hardcopy > $*.ps' -c quit $<
	ps2pdf $*.ps $@
	rm $*.ps
	@echo "Converted "$<" to LaTeX successfully!"

$(BIN)/sd2: $(OBJFILES)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PGO_LDFLAGS) \
	-o $(BIN)/sd2 $(OBJFILES) \
	-static $(BLASLIB) \
	-lgmp $(GMPLIB) $(GMPINC) -lm -lc

.PHONY: pgo
# Do profile generated optimization
# This is recommended for the final executable
pgo:
	./profiledmake.sh

$(BIN)/test_la: $(SRC)/test_la.o $(SRC)/arith.o
	$(CC) $(CFLAGS) -o $(BIN)/test_la $(SRC)/test_la.o $(SRC)/arith.o -lm -lc

tags: $(SRC) Makefile
	ctags $</*.c $</*.h

.PHONY: clean
clean:
	rm -f $(SRC)/*.o $(PDF)/*.pdf dump_lattice.b */dump_lattice.b

