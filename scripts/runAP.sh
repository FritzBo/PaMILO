#!/bin/sh
#
#  Author: Mirko H. Wagner 2020
#  This file is distributed for academics only
#  under the terms of an MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#

python3 scripts/instanceFritzToVLP.py "../$1"
python3 ../../fritsolve/molp-algo/scripts/instanceFritzToVLP.py "$1"

#make -j16
mkdir -p $2

bensolve="bensolve"
cplex=~/cplex/cplex/bin/x86-64_linux/cplex
polyscip=~/Downloads/scip6.0.2-release/bin/applications/polyscip

pamiloOut="$2/$(basename $1)_pamilo"
rm pamiloOut >> /dev/null 2>&1
/usr/bin/time -p ./pamilo_cli $1.lp -o $pamiloOut > $pamiloOut

benOutDual="$2/$(basename $1)_ben_dual"
rm ${benOutDual}* >> /dev/null 2>&1
/usr/bin/time -p $bensolve $1.vlp -Adual -adual -o $benOutDual >> /dev/null

benOutPrimal="$2/$(basename $1)_ben_primal"
rm ${benOutPrimal}* >> /dev/null 2>&1
/usr/bin/time -p $bensolve $1.vlp -Aprimal -aprimal -o $benOutPrimal >> /dev/null

polyscipOut="$2/$(basename $1)_polyscip_out"
rm -f $1.mps
$cplex -c "read $1" "write $1.mps" > /dev/null
$polyscip -x -p params.set "$1.mps" > $polyscipOut


#python3 scripts/diffEps.py ${pamiloOut} ${benOut}_img_p.sol

