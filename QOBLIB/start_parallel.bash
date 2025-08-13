#!/bin/bash

OUTFOLDER=Results

for file in 01-marketsplit/instances/*.dat
do
	echo ""
	echo "INSTANCE ${file} -----------------------------------"
	instance=${file##*/}
	instance="${instance%.dat}"

	# lds and dfs enumeration in parallel

	SECONDS=0
	../bin/sd2 -bkz -beta40 -t1      -o${OUTFOLDER}/${instance}_dfs.sol ${file} >/dev/null 2>>${OUTFOLDER}/parallel_run_dfs_full.log & pids=( $! )
	../bin/sd2 -bkz -beta40  -lds20 -t1 -o${OUTFOLDER}/${instance}_lds.sol ${file} >/dev/null 2>>${OUTFOLDER}/parallel_run_lds_full.log & pids+=( $! )
	wait -n
	kill "${pids[@]}" 2>/dev/null
	duration=$SECONDS
	
	if [ -s ${instance}_dfs.sol ]
	then
		echo "SOL_DFS"
		cat ${OUTFOLDER}/${instance}_dfs.sol
	fi

	if [ -s ${instance}_lds.sol ]
	then
		echo "SOL_LDS"
		cat ${OUTFOLDER}/${instance}_lds.sol
	fi
	echo "time $((duration / 60))m$((duration % 60))s"
done	
