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

import os
import sys

baseDir = os.path.dirname(os.path.realpath(__file__))

inputFile = open(os.path.join(baseDir, sys.argv[1]), "r+")

out = open(os.path.join(baseDir, sys.argv[1]) + ".lp", "w+")

# read number of items (of one partion)
nItems = int(inputFile.readline())

# read number of obj functions
nObjFunc = int(inputFile.readline())

# write program line (information line)
#TODO: check whether min or max
out.write("Minimize multi-objectives\n")
for i in range(nObjFunc):
    out.write("  obj"+str(i)+":\n   z"+str(i)+"\n")
out.write("\nSubject To\n");

# write constraints (each item is used exactly once)
for i in range(nItems):
    out.write("x" + str(i*nItems+nObjFunc))
    for j in range(1,nItems):
        out.write(" + x" + str(i*nItems+j+nObjFunc))
    out.write(" = 1\n")

for i in range(nItems):
    out.write("x" + str(i+nObjFunc))
    for j in range(1,nItems):
        out.write(" + x" + str(j*nItems+i+nObjFunc))
    out.write(" = 1\n")

curObj = 0
lineNo = 0
# read and write obj functions
for line in inputFile.readlines():
    if(lineNo == 0):
        nObjectiveCoeffs = 0
        out.write("-z" + str(curObj))
    lineNo += 1
    line = line.split(" ")
    if(len(line) != nItems):
        print("wrong line size! " + len(line) + " but should be " + nItems + " !\nline is: " + line)
    for coef in line:
        index = nObjectiveCoeffs %(nItems*nItems)
        index += nObjFunc

        out.write(" + " + str(int(coef)) + " x" + str(index))

        nObjectiveCoeffs += 1
    if(lineNo == nItems):
        lineNo = 0
        out.write(" = 0\n")
        curObj += 1

out.write("Bounds\n")
for i in range(nObjFunc):
    out.write("z" + str(i) + " Free\n")

out.write("Binaries\n")
for i in range(nItems):
    for j in range(nItems):
        out.write("x" + str(i*nItems+j+nObjFunc) + " ")
    out.write("\n")

out.write("End\n")

out.close()
inputFile.close()
