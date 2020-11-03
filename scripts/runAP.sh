#!/bin/sh
#
#  Author: Mirko H. Wagner 2020
#  This file is distributed for academics only
#  under the terms of an MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#
bensolve="bensolve"
cplex=~/cplex12_10/cplex/bin/x86-64_linux/cplex
polyscip=~/Downloads/scipoptsuite-6.0.2/build/bin/applications/polyscip

python3 scripts/instanceFritzToVLP.py "../$1"
python3 ../molp-algo/scripts/instanceFritzToVLP.py "$1"
rm -f $1.mps
$cplex -c "read $1.lp" "write $1.mps" > /dev/null

#make -j16
mkdir -p $2

pamiloOutCDD="$2/$(basename $1)_pamilo_cdd"
rm $pamiloOutCDD >> /dev/null 2>&1
/usr/bin/time -p ./pamilo_cli -E cdd $1.lp -o ${pamiloOutCDD} > ${pamiloOutCDD} 2> ${pamiloOutCDD}_time

pamiloOutGL="$2/$(basename $1)_pamilo_graphless"
rm $pamiloOutGL >> /dev/null 2>&1
/usr/bin/time -p ./pamilo_cli -E graphless $1.lp -o ${pamiloOutGL} > ${pamiloOutGL} 2> ${pamiloOutGL}_time

benOutPrimal="$2/$(basename $1)_ben_primal"
rm ${benOutPrimal}* >> /dev/null 2>&1
/usr/bin/time -p ${bensolve} $1.vlp -Aprimal -aprimal -o ${benOutPrimal} >> /dev/null 2> ${benOutPrimal}_time

benOutDual="$2/$(basename $1)_ben_dual"
rm ${benOutDual}* >> /dev/null 2>&1
/usr/bin/time -p ${bensolve} $1.vlp -Adual -adual -o ${benOutDual} >> /dev/null 2> ${benOutDual}_time

polyscipOut="$2/$(basename $1)_polyscip_out"
rm ${polyscipOut}* >> /dev/null 2>&1
/usr/bin/time -p ${polyscip} -x -p params.set "$1.mps" > ${polyscipOut} 2> ${polyscipOut}_time

echo ""


#python3 scripts/diffEps.py ${pamiloOut} ${benOut}_img_p.sol

