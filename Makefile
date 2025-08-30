# CC=clang
CC=gcc
BIN=bin
SRC=src
PDF=pdf

CFLAGS= -O3                  \
  -DUSE_AVX                  \
  -flto                      \
  -fno-math-errno            \
  -ftree-vectorize           \
  -msse3 -mssse3 -msse4.1    \
  -Wall

HASWELL_CFLAGS= -mavx -mavx2 -march=skylake
# ICELAKE_CFLAGS= -mavx512f -mavx512dq -march=                       # Not yet supported
SAPPHIRERAPIDS_CFLAGS= -mavx512f -mavx512dq -march=sapphirerapids    # Not yet supported
ARCH_CFLAGS= $(HASWELL_CFLAGS)

# Profile:   -pg \
#  -g -pg \
# Find memory leaks
#  -fsanitize=address         \

VIMFLAGS=-c 'set printoptions=number:y,left:2pc,right:2pc' -c 'set printfont=Courier:h8'

###################################
# USE BLAS library
#
#
# ------ First option: do not use BLAS at all.
# Uncomment these:
BLAS=NOBLAS
BLASINC=.
BLASLIB=
#
# ------  Second option: use OpenBLAS as installed in e.g. ubuntu.
# Uncomment BLAS, BLASINC and BLASLIB:
# Modern machine:
#BLAS=USE_BLAS
#  ------ Third option: use OpenBLAS installed and compiled in a folder.
# Uncomment these:
#BLAS=USE_BLAS_DEV
#BLASINC=../OpenBLAS/
#BLASLIB=-L../OpenBLAS/ -lopenblas -lpthread
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
GSA_OUT=FALSE
###################################

OBJFILES=$(SRC)/bkz.o   $(SRC)/dio2.o   $(SRC)/dualbkz.o   $(SRC)/enum.o   $(SRC)/lattice.o   $(SRC)/lgs.o   $(SRC)/lll.o   $(SRC)/sd2.o   $(SRC)/arith.o
PDFFILES=$(PDF)/bkz.pdf $(PDF)/dio2.pdf $(PDF)/dualbkz.pdf $(PDF)/enum.pdf $(PDF)/lattice.pdf $(PDF)/lgs.pdf $(PDF)/lll.pdf $(PDF)/sd2.pdf $(PDF)/arith.pdf

all: $(BIN)/sd2 tags $(PDFFILES) $(BIN)/test_la

$(SRC)/%.o: $(SRC)/%.c $(SRC)/%.h $(SRC)/datastruct.h $(SRC)/const.h
	$(CC) $(CFLAGS) $(ARCH_CFLAGS) $(PGO_CFLAGS) -D$(GSA_OUT) -D$(BLAS) -I$(BLASINC) -c $< $(GMPINC) -o $@
	@echo "Compiled "$<" successfully!"

$(BIN)/sd2: $(OBJFILES)
	$(CC) $(CFLAGS) $(ARCH_CFLAGS) $(LDFLAGS) $(PGO_LDFLAGS) \
		-o $(BIN)/sd2 $(OBJFILES) \
		-static $(BLASLIB) \
		-lgmp $(GMPLIB) $(GMPINC) -lm -lc

$(PDF)/%.pdf: $(SRC)/%.c
	vim $(VIMFLAGS) -c 'hardcopy > $*.ps' -c quit $<
	ps2pdf $*.ps $@
	rm $*.ps
	@echo "Converted "$<" to LaTeX successfully!"

# Generate a binary with profile generated optimization
# This is recommended for the final executable
.PHONY: pgo
pgo:
	./profiledmake.sh

$(SRC)/test_la.o: $(SRC)/test_la.c
	$(CC) $(CFLAGS) $(ARCH_CFLAGS) $(PGO_CFLAGS) -D$(GSA_OUT) -D$(BLAS) -I$(BLASINC) -c $< $(GMPINC) -o $@
	@echo "Compiled "$<" successfully!"

$(BIN)/test_la: $(SRC)/test_la.o $(SRC)/arith.o
	$(CC) $(CFLAGS) $(ARCH_CFLAGS) -o $(BIN)/test_la $(SRC)/test_la.o $(SRC)/arith.o -static $(BLASLIB) -lm -lc

tags: $(SRC) Makefile
	ctags $</*.c $</*.h

.PHONY: clean
clean:
	rm -f $(SRC)/*.o $(PDF)/*.pdf dump_lattice.b */dump_lattice.b

