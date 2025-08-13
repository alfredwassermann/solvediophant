#!/bin/bash

OUTFOLDER=Results

for file in 01-marketsplit/instances/*.dat
do
	echo ""
	echo "INSTANCE ${file} -----------------------------------"
	instance=${file##*/}
	instance="${instance%.dat}"

	SECONDS=0
	../bin/sd2 -bkz -beta40 -o${OUTFOLDER}/all_${instance}.sol ${file} >/dev/null 2>>${OUTFOLDER}/all_dfs_full.log
	wc ${OUTFOLDER}/all_${instance}.sol
done	
