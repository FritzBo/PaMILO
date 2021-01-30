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

bensolve="bensolve"
#bensolve="../../fritsolve/fritsolve/bensolve-2.1.0/bensolve"

cplex=~/cplex12_10/cplex/bin/x86-64_linux/cplex
#cplex=~/cplex/cplex/bin/x86-64_linux/cplex

polyscip=~/Downloads/scipoptsuite-6.0.2/build/bin/applications/polyscip
#polyscip=~/Downloads/scip6.0.2-release/bin/applications/polyscip

if [[ $(basename $1) == "ap_"* ]]; then
	python3 scripts/instanceFritzToVLP.py "$1"

	python3 ../molp-algo/scripts/instanceFritzToVLP.py "$1"
	#python3 ../../fritsolve/molp-algo/scripts/instanceFritzToVLP.py "$1"
fi

rm -f $1.mps
$cplex -c "read $1.lp" "write $1.mps" > /dev/null

mkdir -p $2

pamiloOutCDD="$2/$(basename $1)_pamilo_cdd"
/usr/bin/time -p timeout $timelimit ./pamilo_cli -E cdd $1.lp -o ${pamiloOutCDD} > ${pamiloOutCDD} 2> ${pamiloOutCDD}_time

pamiloOutGL="$2/$(basename $1)_pamilo_graphless"
/usr/bin/time -p timeout $timelimit ./pamilo_cli -E graphless $1.lp -o ${pamiloOutGL} > ${pamiloOutGL} 2> ${pamiloOutGL}_time

if [[ $(basename $1) =~ "ap_"* ]]; then
	benOutPrimal="$2/$(basename $1)_ben_primal"
	/usr/bin/time -p timeout $timelimit ${bensolve} $1.vlp -Aprimal -aprimal -o ${benOutPrimal} >> /dev/null 2> ${benOutPrimal}_time

	benOutDual="$2/$(basename $1)_ben_dual"
	/usr/bin/time -p timeout $timelimit ${bensolve} $1.vlp -Adual -adual -o ${benOutDual} >> /dev/null 2> ${benOutDual}_time
fi

polyscipOut="$2/$(basename $1)_polyscip_out"
/usr/bin/time -p timeout $timelimit ${polyscip} -x -p params.set "$1.mps" > ${polyscipOut} 2> ${polyscipOut}_time

echo ""


#python3 scripts/diffEps.py ${pamiloOut} ${benOut}_img_p.sol

