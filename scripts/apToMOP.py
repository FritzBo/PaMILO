#
#  Author: Mirko H. Wagner 2020
#  This file is distributed for academics only
#  under the terms of an MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#

import os
import sys

baseDir = os.path.dirname(os.path.realpath(__file__))

out = open(os.path.join(baseDir, sys.argv[1]) + ".zpl", "w+")

out.write("param prob_file := \""+str(os.path.join(baseDir, sys.argv[1])) + "\";\n")
out.write("param no_objs := read prob_file as \"1n\" use 1;\nparam no_vars := read prob_file as \"1n\" use 1 skip 1;\nset I := {1..no_vars};\nset T := {1..no_objs*no_vars*no_vars};\nparam coeffs[T] := read prob_file as \"n+\" match \"[0-9]+\" skip 2;\nparam offset := no_vars*no_vars;\nparam Obj1[<i,j> in I*I] := coeffs[(i-1)*no_vars + j];\nparam Obj2[<i,j> in I*I] := coeffs[(i-1)*no_vars + j + offset];\nparam Obj3[<i,j> in I*I] := coeffs[(i-1)*no_vars + j + 2*offset];\nvar x[I*I] binary;\nminimize Obj1: sum <i,j> in I*I: Obj1[i,j]*x[i,j];\nObj2: sum <i,j> in I*I: Obj2[i,j]*x[i,j];\nObj3: sum <i,j> in I*I: Obj3[i,j]*x[i,j];\nsubto row: forall <i> in I dosum <j> in I: x[i,j] == 1;\nsubto col: forall <i> in I dosum <j> in I: x[j,i] == 1;\n");


