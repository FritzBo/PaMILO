#
#  Author: Mirko H. Wagner 2020
#  This file is distributed for academics only
#  under the terms of an MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#

import os
import sys
import random

baseDir = os.getcwd() #path.dirname(os.path.realpath(__file__))

generate = False
if len(sys.argv) > 1:
    generate = True
else:
    inputFile = open(os.path.join(baseDir, sys.argv[1]), "r+")


# read number of obj functions
if generate:
    nObjFunc = int(sys.argv[2])
    nCols = int(sys.argv[3])
    nRows = int(sys.argv[4])
    intVarCutoff = int(sys.argv[5])
    outName = ""
    counter = -1
    while counter == -1 or os.path.isfile(outName):
        counter += 1
        outName = os.path.join(baseDir, sys.argv[1] + "_p-" + str(nObjFunc) + "_n-" + str(nCols) + "_m-" + str(nRows) + "_i-" + str(intVarCutoff) + "_ins-" + str(counter) + ".lp")
else:
    nObjFunc = int(inputFile.readline())
    nCols = int(inputFile.readline())
    nRows = int(inputFile.readline())
    intVarCutoff = nCols
    outName = os.path.join(baseDir, sys.argv[1] + ".lp")

out = open(outName, "w+")

objCoeffs = []
for i in range(nObjFunc):
    if not generate:
        row = inputFile.readline().split(',')
    objCoeffs.append([])
    for j in range(nCols):
        if generate:
            if random.randint(1,5) == 1:
                val = -random.randint(1, 100)
            else:
                val = random.randint(0,100)
        else:
            val = float(row[j].strip('[]\n'))
        objCoeffs[i].append(val)

out.write("Minimize multi-objectives\n")
for i in range(nObjFunc):
    out.write("  obj"+str(i)+":\n   z"+str(i)+"\n")

coeffs = []
rhs = []
for i in range(nRows):
    if not generate:
        row = inputFile.readline().split(',')
    else:
        rhs.append(0)
    coeffs.append([])
    for j in range(nCols):
        if generate:
            randVar = random.randint(1,10)
            if randVar == 1:
                val = -random.randint(1, 100)
            elif randVar == 2:
                val = 0
            else:
                val = random.randint(1,100)
            rhs[i] += val
        else:
            val = float(row[j].strip('[]\n'))
        coeffs[i].append(val)

if not generate:
    row = inputFile.readline().split(',')
for i in range(nRows):
    if not generate:
        rhs.append(float(row[i].strip('[]\n')))
    else:
        if rhs[i] > 100:
            rhs[i] = random.randint(100, rhs[i])
        else:
            rhs[i] = random.randint(rhs[i], 100)

out.write("Subject To\n");
# write constraints (each item is used exactly once)
for i in range(nRows):
    nSuccCoefsInRow = 0
    for j in range(nCols):
        out.write(" ")
        if not coeffs[i][j] == 0:
            if nSuccCoefsInRow > 0 and coeffs[i][j] >= 0:
                out.write("+ ")
            elif coeffs[i][j] < 0:
                out.write("- ")
                coeffs[i][j] *= -1
            out.write(str(coeffs[i][j]) + " x"+str(j))
            nSuccCoefsInRow += 1
    if nSuccCoefsInRow > 0 :
        out.write(" <= " + str(rhs[i]) + "\n")

for i in range(nObjFunc):
    out.write(" z"+str(i))
    for j in range(nCols):
        out.write(" ")
        if objCoeffs[i][j] >= 0:
            out.write("+ ")
        else:
            out.write("- ")
            objCoeffs[i][j] *= -1
        out.write(str(objCoeffs[i][j]) + " x"+str(j))
    out.write(" = 0\n")

out.write("Bounds\n")
for i in range(nCols):
    out.write(" x" + str(i) + " >= 0\n")
for i in range(nObjFunc):
    out.write(" z" + str(i) + " Free\n")

out.write("GENERAL\n")
for i in range(intVarCutoff):
    out.write(" x" + str(i))
    if i%10 == 9 or i+1 == intVarCutoff:
        out.write("\n")

out.write("End\n")

out.close()
if not generate:
    inputFile.close()

