#!/bin/sh
python3 scripts/instanceFritzToVLP.py "../$1"
python3 ../molp-algo/scripts/instanceFritzToVLP.py "$1"

#make -j16
mkdir -p $2

pamiloOut="$2/$(basename $1)_pamilo"
rm pamiloOut >> /dev/null 2>&1
/usr/bin/time -p ./pamilo_cli $1.lp -o $pamiloOut > $pamiloOut

benOut="$2/$(basename $1)_ben"
rm ${benOut}* >> /dev/null 2>&1
/usr/bin/time -p ../fritsolve/bensolve-2.1.0/bensolve $1.vlp -o $benOut >> /dev/null

python3 scripts/diffEps.py ${pamiloOut} ${benOut}_img_p.sol
