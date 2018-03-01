# IMPORT
import sys

# MAIN
# Reads a "*.summed" or "*.read_summary" file [1] 
# and splits the SNPs in to seperate files
# by contigs (1st column)

# Read the input file
in_lines = open(sys.argv[1],"r").readlines()

# Get all contigs
contigs = list(set([l.strip().split("\t")[0] for l in in_lines[1:]]))

# For each contig select line and write an outputfile 
for c in contigs:
    outlines = [l for l in in_lines if l.strip().split("\t")[0] == c]
    output = open(sys.argv[1]+"."+c,"w")
    output.write("".join(outlines))
    output.close()
