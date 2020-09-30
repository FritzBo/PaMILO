#!/bin/sh

mkdir -p $2

out="$2/$(basename $1)_polyscip_out"
cplex=~/cplex/cplex/bin/x86-64_linux/cplex
polyscip=~/Downloads/scip6.0.2-release/bin/applications/polyscip

rm -f $1.mps
$cplex -c "read $1" "write $1.mps" > /dev/null
$polyscip -x -p params.set "$1.mps" > $out
