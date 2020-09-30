import os
for fileDir in ["./include", "./src", "./cli", "./cmake", "./doc", "./scripts"]:
    for root, dirs, files in os.walk(fileDir, topdown=False):
        for name in files:
            if ".md" in name or ".py" in name or ".sh" in name:
                print(os.path.join(root, name))
                f = open(os.path.join(root, name), "r+")
                lines = f.readlines()
                f.close()
                f = open(os.path.join(root, name), "w+")
                i = 0
                while i < len(lines) and not "This file is distributed u" in lines[i]:
                    f.write(lines[i])
                    i += 1

                if i < len(lines):
                    f.write("#  This file is distributed for academics only\n#  under the terms of an MIT license based license,\n")
                    i += 8

                while i < len(lines):
                    f.write(lines[i])
                    i += 1
