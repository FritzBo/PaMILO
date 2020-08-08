#!/bin/sh
python3 scripts/instanceFritzToVLP.py "../$1"
python3 ../molp-algo/scripts/instanceFritzToVLP.py "$1"

#make -j16
mkdir -p $2

bensolve="bensolve"

pamiloOut="$2/$(basename $1)_pamilo"
rm pamiloOut >> /dev/null 2>&1
/usr/bin/time -p ./cli/pamilo_cli $1.lp > $pamiloOut

benOut="$2/$(basename $1)_ben"
rm ${benOut}* >> /dev/null 2>&1
/usr/bin/time -p $bensolve -Adual -adual $1.vlp -o $benOut >> /dev/null

python3 ../molp-algo/scripts/diffEps.py ${pamiloOut} ${benOut}_img_p.sol
