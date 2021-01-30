#!/bin/sh
#
#  Author: Mirko H. Wagner 2020
#  This file is distributed for academics only
#  under the terms of an MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#
#bensolve=~/bensolve-2.1.0/bensolve
#cplex=~/cplex12_10/cplex/bin/x86-64_linux/cplex
#polyscip=~/Downloads/scipoptsuite-6.0.2/build/bin/applications/polyscip

timelimit="1h"

bensolve="../bensolve-2.1.0/bensolve"
#bensolve="../../fritsolve/fritsolve/bensolve-2.1.0/bensolve"

cplex=~/cplex12_10/cplex/bin/x86-64_linux/cplex
#cplex=~/cplex/cplex/bin/x86-64_linux/cplex

polyscip=~/Downloads/scipoptsuite-6.0.2/build/bin/applications/polyscip
#polyscip=~/Downloads/scip6.0.2-release/bin/applications/polyscip

inst=$(basename $1)

if [[ ${inst} == "ap_"* ]]; then
	python3 scripts/instanceFritzToVLP.py "$1"

	python3 ../molp-algo/scripts/instanceFritzToVLP.py "$1"
	#python3 ../../fritsolve/molp-algo/scripts/instanceFritzToVLP.py "$1"
fi

rm -f $1.mps
$cplex -c "read $1.lp" "write $1.mps" > /dev/null

mkdir -p $2

eps="1e-9"

#pamiloOutCDDeps="$2/${inst}_pamilo_cdd_eps"
#timeout $timelimit /usr/bin/time -p ./pamilo_cli -E cdd -e $eps -v "1e-7" $1.lp -o ${pamiloOutCDDeps} > ${pamiloOutCDDeps} 2> ${pamiloOutCDDeps}_time

pamiloOutCDD="$2/${inst}_pamilo_cdd"
/usr/bin/time -p timeout $timelimit ./pamilo_cli -E cdd -e $eps $1.lp -o ${pamiloOutCDD} > ${pamiloOutCDD} 2> ${pamiloOutCDD}_time

pamiloOutGL="$2/${inst}_pamilo_graphless"
/usr/bin/time -p timeout $timelimit ./pamilo_cli -E graphless -e $eps $1.lp -o ${pamiloOutGL} > ${pamiloOutGL} 2> ${pamiloOutGL}_time

if [[ ${inst} =~ "ap_"* ]]; then
	benOutPrimal="$2/${inst}_ben_primal"
	/usr/bin/time -p timeout $timelimit ${bensolve} $1.vlp -Aprimal -aprimal -o ${benOutPrimal} >> /dev/null 2> ${benOutPrimal}_time

	benOutDual="$2/${inst}_ben_dual"
	/usr/bin/time -p timeout $timelimit ${bensolve} $1.vlp -Adual -adual -o ${benOutDual} >> /dev/null 2> ${benOutDual}_time
fi

#polyscipOut="$2/${inst}_polyscip_out"
#/usr/bin/time -p timeout $timelimit ${polyscip} -x -p params.set "$1.mps" > ${polyscipOut} 2> ${polyscipOut}_time

echo ""


cd $2
rm ${inst}_ben_*_adj*
rm ${inst}_ben_*_inc*
rm ${inst}_ben_*_img*
rm ${inst}_pamilo_*_cplex
rm ${inst}_pamilo_*_sol
zip -rv ${inst}.zip ${inst}*
rm ${inst}_*
