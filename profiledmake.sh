#!/bin/sh
PGO_CFLAGS="-g -fprofile-generate -fprofile-update=prefer-atomic" PGO_LDFLAGS="-lgcov --coverage" make clean bin/sd2
rm src/*.gcda
bin/sd2 -bkz -beta40 Tests/KM_PSL_2_23_7_8.in 
PGO_CFLAGS="-fprofile-use -fprofile-correction -fprofile-partial-training -static" make clean bin/sd2