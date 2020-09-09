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
import numpy as np
import pandas
import math
import time
import numpy_indexed as npi

np.set_printoptions(threshold=np.inf)

def eps(num1, num2):
    factor = 1.1
    diffEps = 10
    if num1 * factor >= num2 and num2 * factor >= num1:
        return True
    if num1 + diffEps >= num2 and num2 + diffEps >= num1:
        return True
    return False


baseDir = os.path.dirname(os.path.realpath(__file__))

origIn1 = pandas.read_csv(sys.argv[1], sep=' ');
nCols1 = len(origIn1.columns)
origIn1 = pandas.read_csv(sys.argv[1], sep=' ', dtype="float64", names=["type"] + ["dim"+str(i) for i in range(nCols1-1)]);
origIn1 = origIn1.round(0)
origIn1 = origIn1.drop_duplicates()

origIn2 = pandas.read_csv(sys.argv[2], sep=' ');
nCols2 = len(origIn2.columns)
origIn2 = pandas.read_csv(sys.argv[2], sep=' ', dtype="float64", names=["type"] + ["dim"+str(i) for i in range(nCols2-1)]);
origIn2 = origIn2.round(0)
origIn2 = origIn2.drop_duplicates()

if not nCols1 == nCols2 :
    print("Different number of dimensions!")
    exit(0)

inDf = origIn1.merge(origIn2, how='inner')

diffFactor = 1
if inDf.shape[0] == origIn1.shape[0] and inDf.shape[0] == origIn2.shape[0]:
    print("-\t-", end="",sep="")
else:
    print(origIn1.shape[0] - inDf.shape[0], "\t", origIn2.shape[0] - inDf.shape[0], end="",sep="")
if not (origIn1.shape[0] * diffFactor >= origIn2.shape[0] and origIn1.shape[0] <= origIn2.shape[0] * diffFactor) :
    print("\t",origIn1.shape[0],"\t",origIn2.shape[0], sep="")
else:
    print("\t-")

