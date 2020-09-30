#!/bin/sh
#
#  Author: Mirko H. Wagner 2020
#  This file is distributed for academics only
#  under the terms of an MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#

make -j16

dir="results/"$(date +%m_%d_%Y_%H_%M_%S)
mkdir -p "results"
mkdir $dir

ulimit -v 12388608
fileList=$(ls ../../fritsolve/molp-algo/instances/fritzins | grep "ap_3_[4-8][0-9]*_[0-9]*$\|ap_4_[1-4][0-9]_[0-9]*$\|ap_5_[1-2][0-9]_[0-9]*$\|ap_6_[0-9]_[0-9]*$\|ap_6_1[0-2]_[0-9]*$")

for file in $fileList
do
	#echo $file
	if [ "$RANDOM" -le 200 ]
	then
		echo $file
		./scripts/runAP.sh ../../fritsolve/molp-algo/instances/fritzins/$file $dir
		./scripts/runPolyScip.sh ../../fritsolve/molp-algo/instances/fritzins/$file.lp
	fi
done

