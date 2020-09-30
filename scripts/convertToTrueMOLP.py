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

out = open(os.path.join(baseDir, sys.argv[1][:-3]) + "_converted.lp", "w+")

for line in inputFile.readlines():
    if "obj:" in line:
        i = 1
        line = line.split(" ")
        for word in line:
            if "z" in word:
                out.write(" obj" + str(i) + ":\r\n  z" + str(i) + "\r\n")
                i += 1
    elif "imize" in line or "min" in line or "max" in line:
        out.write(line[:-1] + " multi-objectives\r\n")
    else:
        out.write(line[:-1] + "\r\n")

inputFile.close()
out.close()
