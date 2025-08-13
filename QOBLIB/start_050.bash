#!/bin/bash

OUTFOLDER=Results

for file in 01-marketsplit/instances/ms_*_050_*.dat
do
	echo ""
	echo "INSTANCE ${file} -----------------------------------"
	instance=${file##*/}
	instance="${instance%.dat}"

	(time ../bin/sd2 -bkz -beta56 -tours20 -lds20 -t1 -o${OUTFOLDER}/050_${instance}.sol ${file} >/dev/null) 2>>${OUTFOLDER}/050_lds.log
done
