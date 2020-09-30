#
#  Author: Mirko H. Wagner 2020
#  This file is distributed for academics only
#  under the terms of an MIT license based license,
#  a copy of which can be found in the file LICENSE-academic.txt.
#

import os
import sys

baseDir = os.path.dirname(os.path.realpath(__file__))

inputFile = open(os.path.join(baseDir, sys.argv[1]), "r+")

out = open(os.path.join(baseDir, sys.argv[1]) + ".mop", "w+")

# read number of items (of one partion)
nItems = int(inputFile.readline())

# read number of obj functions
nObjFunc = int(inputFile.readline())

# write program line (information line)
#TODO: check whether min or max
out.write("NAME ASSIGNMENT\nOBJSENSE\n MAX\nROWS\n")
for i in range(nObjFunc):
    out.write(" N obj"+str(i)+"\n")

# write constraints (each item is used exactly once)
for i in range(nItems):
    out.write(" "+'{0: <5}'.format(" E cLeft#" + str(i) + "\n")
    out.write(" E cRight#" + str(i) + "\n")

out.write("COLUMNS\n")
for i in range(nItems):
    for j in range(nItems):
        out.write("    x#"+str(i)+"#"+str(j)+"    cLeft#"+str(i)+"    1\n")
        out.write("    x#"+str(i)+"#"+str(j)+"    cRight#"+str(j)+"    1\n")

# read and write obj functions
lineNo = 0
curObj = 0
for line in inputFile.readlines():
    if(lineNo == 0):
        nObjectiveCoeffs = 0
    lineNo += 1
    line = line.split(" ")
    if(len(line) != nItems):
        print("wrong line size! " + len(line) + " but should be " + nItems + " !\nline is: " + line)

    noInLine=0
    for coef in line:
        out.write("    x#"+str(noInLine)+"#"+str(lineNo)+"    obj#"+str(curObj)+"    "+str(int(coef))+"\n")
        noInLine += 1

    if(lineNo == nItems):
        lineNo = 0
        curObj += 1

out.write("RHS\n")
for i in range(nItems):
    out.write("    RHS       cLeft#"+str(i)+"    1\n")
    out.write("    RHS       cRight#"+str(i)+"    1\n")

out.write("BOUNDS\n")
for i in range(nItems):
    for j in range(nItems):
        out.write(" BV    BOUND    x#"+str(i)+"#"+str(j)+"\n")
out.write("ENDATA\n")

out.close()
inputFile.close()

