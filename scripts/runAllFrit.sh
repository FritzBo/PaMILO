#!/bin/sh
make -j16

ulimit -v 12388608
fileList=$(ls ../molp-algo/instances/fritzins | grep "ap_3_[4-8][0-9]*_[0-9]*$\|ap_4_[1-4][0-9]_[0-9]*$\|ap_5_[1-2][0-9]_[0-9]*$\|ap_6_[0-9]_[0-9]*$\|ap_6_1[0-2]_[0-9]*$")
for file in $fileList
do
	#echo $file
	if [ "$RANDOM" -le 2000 ]
	then
		echo $file
		./scripts/runAP.sh ../molp-algo/instances/fritzins/$file
	fi
done
