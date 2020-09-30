#! /bin/sh

infile=$(pwd)
infile+="/$1"

out=$infile
out+=".zpl"

temp=".temp"

cat $infile | sed -r '1,2!b;1h;1!H;2!d;x;s/^([^\n]*)(.*\n)(.*)/\3\2\1/' > $temp

echo "param prob_file := \"$(pwd)/$temp\";
param no_objs := read prob_file as \"1n\" use 1;
param no_vars := read prob_file as \"1n\" use 1 skip 1;
set I := {1..no_vars};
set T := {1..no_objs*no_vars*no_vars};
param coeffs[T] := read prob_file as \"n+\" match \"[0-9]+\" skip 2;
param offset := no_vars*no_vars;
param Obj1[<i,j> in I*I] := coeffs[(i-1)*no_vars + j];
param Obj2[<i,j> in I*I] := coeffs[(i-1)*no_vars + j + offset];
param Obj3[<i,j> in I*I] := coeffs[(i-1)*no_vars + j + 2*offset];
param Obj4[<i,j> in I*I] := coeffs[(i-1)*no_vars + j + 3*offset];
param Obj5[<i,j> in I*I] := coeffs[(i-1)*no_vars + j + 4*offset];
param Obj6[<i,j> in I*I] := coeffs[(i-1)*no_vars + j + 5*offset];
var x[I*I] binary;
minimize Obj1: sum <i,j> in I*I: Obj1[i,j]*x[i,j];
Obj2: sum <i,j> in I*I: Obj2[i,j]*x[i,j];
Obj3: sum <i,j> in I*I: Obj3[i,j]*x[i,j];
Obj4: sum <i,j> in I*I: Obj4[i,j]*x[i,j];
Obj5: sum <i,j> in I*I: Obj5[i,j]*x[i,j];
Obj6: sum <i,j> in I*I: Obj6[i,j]*x[i,j];
subto row: forall <i> in I do
	sum <j> in I: x[i,j] == 1;
subto col: forall <i> in I do
	sum <j> in I: x[j,i] == 1;
" > $out

python3 ~/Downloads/scipoptsuite-6.0.2/scip/applications/PolySCIP/mult_zimpl/mult_zimpl_to_mop.py $out -p "$(dirname $infile)" -o "$(basename $infile)" --path_to_zimpl ~/Downloads/scip6.0.2-release/bin
