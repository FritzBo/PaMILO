#!/bin/sh
make -j16 -C ~/Downloads/cddlib 1>/dev/null 2>/dev/null
make install -C ~/Downloads/cddlib >/dev/null

make -j16

ulimit -v 12388608
fileList=$(ls ../molp-algo/instances/fritzins | grep "ap_[0-9]*_[0-9]*_[0-9]*$")
for file in $fileList
do
	#echo $file
	if [ "$RANDOM" -le 2000 ]
	then
		echo $file
		timeout 90s ./runAP.sh ../molp-algo/instances/fritzins/$file
	fi	
done
