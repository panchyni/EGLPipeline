# IMPORT
import sys

# MAIN
print'''

Reads a list of control files [1] and the path
to "0_ComparsionToInocWrapper_vX.py" [2]
and writes *.cc file for qsub_hpc

'''

files = [ln.strip() for ln in open(sys.argv[1],"r").readlines()]
path = sys.argv[2]

command_lines = []

for f in files:
    command_lines.append("python " + path + " " + f + "\n")

output = open(sys.argv[1] + ".cc","w")
output.write("".join(command_lines))
output.close()
