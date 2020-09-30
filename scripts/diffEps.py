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
#origIn1 = origIn1.to_numpy()
origIn1 = origIn1.round(0)
origIn1 = origIn1.drop_duplicates()
#origIn1 = origIn1[np.lexsort(np.rot90(origIn1))]
#
#row_mask_next = np.append([True], np.any(np.diff(origIn1,axis=0),1))
#row_mask_prev = np.append([True], np.any(np.diff(np.flipud(origIn1),axis=0),1))

#origIn1 = origIn1[np.logical_and(row_mask_next, np.flipud(row_mask_prev))]
#print(origIn1.shape[0])

origIn2 = pandas.read_csv(sys.argv[2], sep=' ');
nCols2 = len(origIn2.columns)
origIn2 = pandas.read_csv(sys.argv[2], sep=' ', dtype="float64", names=["type"] + ["dim"+str(i) for i in range(nCols2-1)]);
#origIn2 = origIn2.to_numpy()
origIn2 = origIn2.round(0)
origIn2 = origIn2.drop_duplicates()
#origIn2 = origIn2[np.lexsort(np.rot90(origIn2))]
#
#row_mask_next = np.append([True], np.any(np.diff(origIn2,axis=0),1))
#row_mask_prev = np.append([True], np.any(np.diff(np.flipud(origIn2),axis=0),1))
#
#origIn2 = origIn2[np.logical_and(row_mask_next, np.flipud(row_mask_prev))]
#print(origIn2.shape[0])

if not nCols1 == nCols2 :
    print("Different number of dimensions!")
    exit(0)

inDf = origIn1.merge(origIn2, how='inner')
#startRoundWidth = 0
#roundWidth = startRoundWidth
#rwEqual = -1
#
#in1 = origIn1
#in2 = origIn2
#
#
#roundWidth = startRoundWidth
#
#while True:
#    inDf = inDf[np.lexsort(np.rot90(inDf))]
#    #print(inDf.shape[0])
#
#    row_mask_next = np.append([True], np.any(np.diff(inDf,axis=0),1))
#    row_mask_prev = np.append([True], np.any(np.diff(np.flipud(inDf),axis=0),1))
#
#    inDf = inDf[np.logical_and(row_mask_next, np.flipud(row_mask_prev))]
#
#    inDf = inDf.round(roundWidth)
#
#    if inDf.size == 0:
#        rwEqual = roundWidth
#        break
#
#    if roundWidth < 1:
#        break
#
#    roundWidth = int((roundwitdh-1)/2)

diffFactor = 1
if inDf.shape[0] == origIn1.shape[0] and inDf.shape[0] == origIn2.shape[0]:
    print("-\t-", end="",sep="")
else:
    print(origIn1.shape[0] - inDf.shape[0], "\t", origIn2.shape[0] - inDf.shape[0], end="",sep="")
if not (origIn1.shape[0] * diffFactor >= origIn2.shape[0] and origIn1.shape[0] <= origIn2.shape[0] * diffFactor) :
    print("\t",origIn1.shape[0],"\t",origIn2.shape[0], sep="")
else:
    print("\t-")
