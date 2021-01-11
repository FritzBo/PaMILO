#!/bin/sh
#
#  Author: Mirko H. Wagner 2020
#  This file is distributed for academics only
#  under the terms of an MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#

make -j16

fileList=""

resultDir=$(date +%Y-%m-%d_%H-%M-%S)
mkdir -p "results"
mkdir "results/$resultDir"

ulimit -v 12388608

# --- ap ---
#apDir="instances/fritz15"
apDir="../molp-algo/instances/fritzins/"
apfiles=$(find $apDir | grep "[0-9]$")
#apfiles=$(echo $apfiles | grep "ap_3_[4-8][0-9]*_[0-9]*$\|ap_4_[1-4][0-9]_[0-9]*$\|ap_5_[1-2][0-9]_[0-9]*$\|ap_6_[0-9]_[0-9]*$\|ap_6_1[0-2]_[0-9]*$")
fileList="$fileList $apfiles"

# --- kim ---
#kimfiles="$fileList $(ls instances | grep "TestFil_n.*lp\|TestFil-.*lp" | rev | cut -c 4- |rev)
#fileList="$fileList $kimfiles"

# --- kirlik ---
addKirlikInstances() {
	fileList="$(find $1 | grep "\.lp$" | rev | cut -c 4- |rev)\n$fileList"
}
addKirlikInstances "instances/kirlik/ILP"
addKirlikInstances "instances/kirlik/MILPgenerated"
addKirlikInstances "instances/kirlik/KP"


#echo $fileList
#echo $dirList


#parallel -j8 'sh ./scripts/runAP.sh {2}/{1} results/{3}' ::: $fileList ::: $instancedir ::: $resultDir
parallel -j8 'zsh ./scripts/runAP.sh {1} results/{2} && inst=$(basename {1}) && \
	cd results && zip -rv {2}/${inst}.zip {2}/${inst}* && rm {2}/${inst}_*' ::: $fileList ::: $resultDir