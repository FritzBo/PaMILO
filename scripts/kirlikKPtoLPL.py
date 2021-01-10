#
#  Author: Mirko H. Wagner 2020
#  This file is distributed for academics only
#  under the terms of an MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#

import os
import sys
import random

baseDir = os.path.dirname(os.path.realpath(__file__))

generate = False
if len(sys.argv) > 1:
    generate = True
else:
    inputFile = open(os.path.join(baseDir, sys.argv[1]), "r+")

out = open(os.path.join(baseDir, sys.argv[1]) + ".lp", "w+")

# read number of obj functions
if generate:
    nObjFunc = int(sys.argv[2])
    nItems = int(sys.argv[3])
    cap = 0
else:
    nObjFunc = int(inputFile.readline())
    nItems = int(inputFile.readline())
    cap = int(inputFile.readline())

objCoeffs = []
for i in range(nObjFunc):
    objCoeffs.append([])
    if not generate:
        row = inputFile.readline().split(',')
    for j in range(nItems):
        if generate:
            val = random.randint(1,1000)
        else:
            val = int(row[j].strip('[]\n'))
        objCoeffs[i].append(val)

out.write("Minimize multi-objectives\n")
for i in range(nObjFunc):
    out.write("  obj"+str(i)+":\n   z"+str(i)+"\n")

weights = []
if not generate:
    row = inputFile.readline().split(',')
for i in range(nItems):
    if generate:
        val = random.randint(1,1000)
        cap += val/2.0
    else:
        val = int(row[j].strip('[]\n'))
    weights.append(val)

out.write("Subject To\n");
# write constraints (each item is used exactly once)
for i in range(nItems):
    out.write(" ")
    if i > 0 and weights[i] >= 0:
        out.write("+ ")
    elif weights[i] < 0:
        out.write("- ")
        weights[i] *= -1
    out.write(str(weights[i]) + " x"+str(i))
out.write(" <= " + str(cap) + "\n")

for i in range(nObjFunc):
    out.write(" z"+str(i))
    for j in range(nItems):
        out.write(" ")
        if objCoeffs[i][j] >= 0:
            out.write("+ ")
        else:
            out.write("- ")
            objCoeffs[i][j] *= -1
        out.write(str(objCoeffs[i][j]) + " x"+str(j))
    out.write(" = 0\n")

out.write("Bounds\n")
for i in range(nObjFunc):
    out.write(" z" + str(i) + " Free\n")

out.write("BINARY\n")
for i in range(nItems):
    out.write(" x" + str(i))
    if i%10 == 9 or i+1 == nItems:
        out.write("\n")

out.write("End\n")

out.close()
if not generate:
    inputFile.close()

