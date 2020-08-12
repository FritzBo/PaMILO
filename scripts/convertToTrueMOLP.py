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
