#!/bin/sh
make -j16


dir="results/"$(date +%m_%d_%Y_%H_%M_%S)
mkdir -p "results"
mkdir $dir

ulimit -v 12388608
fileList=$(ls ../molp-algo/instances/fritzins | grep "ap_[0-9]*_[0-9]*_[0-9]*$")
for file in $fileList
do
	#echo $file
	if [ "$RANDOM" -le 2000 ]
	then
		echo $file
		./scripts/runAP.sh ../molp-algo/instances/fritzins/$file $dir
	fi	
done
