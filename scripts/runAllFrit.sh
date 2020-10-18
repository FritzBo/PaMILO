#!/bin/sh
#
#  Author: Mirko H. Wagner 2020
#  This file is distributed for academics only
#  under the terms of an MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#

instancedir=../molp-algo/instances/fritzins

make -j16

dir=$(date +%Y-%m-%d_%H-%M-%S)
resDir="results/$dir"
mkdir -p "results"
mkdir $resDir
touch $resDir/configsRun
cd results && zip $dir.zip $dir/configsRun && cd ..

ulimit -v 12388608
fileList=$(ls $instancedir | grep "ap_3_[4-8][0-9]*_[0-9]*$\|ap_4_[1-4][0-9]_[0-9]*$\|ap_5_[1-2][0-9]_[0-9]*$\|ap_6_[0-9]_[0-9]*$\|ap_6_1[0-2]_[0-9]*$")

for file in $fileList
do
	#echo $file
	if [ "$RANDOM" -le 50000 ]
	then
		echo $file
		echo $file >> $resDir/configsRun
		bash ./scripts/runAP.sh $instancedir/$file $resDir
		cd results && zip -rv $dir.zip $dir/$file* && cd ..
		rm $resDir/$file*
	fi
done

cd results && zip -rv $dir.zip $dir/configsRun && cd ..
rm -rf $resDir
