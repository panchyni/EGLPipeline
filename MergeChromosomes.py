# IMPORT
import sys

# MAIN
print'''

Takes the suffix of a set of "CombinedMatrix_chr##" files [1]
and the number of chromsomes to use [2]

Optionally, add chromsome #s to skip as a comma seperated list [3]

Merges filescross chromosomes and replaces the header
line with spacers for visulization

'''
# Read inputs
file_suffix = sys.argv[1]
chr_num = sys.argv[2]
skip_values = []
if len(sys.argv) > 3:
    skip_values = [int(v) for v in sys.argv[3].split(",")] 
    print "Skipping " + ",".join([str(v) for v in skip_values])

chrs = ["chr"+str(n+1) for n in range(int(chr_num)) if not n+1 in skip_values]

# Read and merge files
outlines = []
for chr in chrs:
    file = "CombinedMatrix_" + chr + "_" + file_suffix
    lines = [ln for ln in open(file,"r").readlines()]
    lines[0] = chr + "\t" + "\t".join([v.split("/")[-1].split(".")[0] for v in lines[0].split("\t")[1:]]) + "\n"
    outlines.extend(lines)

# Write output
output_file = "CombinedMatrix_" + file_suffix + ".cross_genome"
output = open(output_file,"w")
output.write("".join(outlines))
output.close()

print "Outfile: " + output_file
