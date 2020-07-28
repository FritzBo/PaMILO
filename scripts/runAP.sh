#!/bin/sh
python3 instanceFritzToVLP.py "$1"
#make -j16
/usr/bin/time -p ./cli/mco_cli $1.lp
