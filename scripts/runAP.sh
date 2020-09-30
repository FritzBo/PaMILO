#!/bin/sh
#
#  Author: Mirko H. Wagner 2020
#  This file is distributed under the terms of
#
#  the GNU General Public License v3,
#  a copy of which can be found in the file LICENCE-GPLv3.txt
#
#  OR
#
#  for academics, a MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#

python3 scripts/instanceFritzToVLP.py "../$1"
python3 ../molp-algo/scripts/instanceFritzToVLP.py "$1"

#make -j16
mkdir -p $2

bensolve="bensolve"

pamiloOut="$2/$(basename $1)_pamilo"
rm pamiloOut >> /dev/null 2>&1
/usr/bin/time -p ./pamilo_cli $1.lp -o $pamiloOut > $pamiloOut

benOut="$2/$(basename $1)_ben"
rm ${benOut}* >> /dev/null 2>&1
/usr/bin/time -p $bensolve $1.vlp -Adual -adual -o $benOut >> /dev/null

python3 scripts/diffEps.py ${pamiloOut} ${benOut}_img_p.sol
