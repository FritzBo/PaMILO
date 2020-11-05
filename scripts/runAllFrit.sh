#!/bin/sh
#
#  Author: Mirko H. Wagner 2020
#  This file is distributed for academics only
#  under the terms of an MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#

#instancedir=../../fritsolve/molp-algo/instances/fritzins
instancedir=instances

make -j16

dir=$(date +%Y-%m-%d_%H-%M-%S)
resDir="results/$dir"
mkdir -p "results"
mkdir $resDir

ulimit -v 12388608
#fileList=$(ls $instancedir | grep "ap_3_[4-8][0-9]*_[0-9]*$\|ap_4_[1-4][0-9]_[0-9]*$\|ap_5_[1-2][0-9]_[0-9]*$\|ap_6_[0-9]_[0-9]*$\|ap_6_1[0-2]_[0-9]*$")
fileList=$(ls instances | grep "TestFil_n.*lp\|TestFil-.*lp" | rev | cut -c 4- |rev)

parallel -j8 'sh ./scripts/runAP.sh {2}/{1} results/{3}' ::: $fileList ::: $instancedir ::: $dir
#cd results && zip -rv {3}/{1}.zip {3}/{1}* && cd .. && \
	#rm results/{3}/{1}_*
