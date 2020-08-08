#!/bin/sh
python3 scripts/instanceFritzToVLP.py "../$1"
python3 ../molp-algo/scripts/instanceFritzToVLP.py "$1"

#make -j16

rm .temp >> /dev/null 2>&1
/usr/bin/time -p ./cli/mco_cli $1.lp >> .temp
rm $1_* >> /dev/null 2>&1
/usr/bin/time -p $bensolve -Adual -adual $1.vlp >> /dev/null

python3 ../molp-algo/scripts/diffEps.py .temp $1_img_p.sol
