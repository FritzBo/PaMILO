#!/bin/sh
configs=$(ls instances/kirlik/ILP | grep dat$ | sed 's/_ins-[0-9]*.dat//' | uniq | sed 's/ILP_p-//' | sed 's/_n-/,/' | sed 's/_m-/,/')

for config in $configs; do
	nInts=$(expr $(echo $config  | sed 's/^[^,]*,//' | sed 's/,[^,]*$//') / 2)
	for i in $(seq 1 $1); do
		python3 scripts/kirlikToLP.py instances/kirlik/MILPgenerated/MILP $(echo $config | sed 's/,/ /g') $nInts
	done
done
