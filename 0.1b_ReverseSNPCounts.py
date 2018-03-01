# IMPORT
import sys

# MAIN
#print'''
''' 
 Reads that output of "1a_FilterVCFSummaryFile.py" and reverse
 the base calls, counts, and allele frequency for ref and SNP

 Inputs:
   -filtered_SNPs = filtered SNP list from "1a_FilterVCFSummaryFile.py"

 This is used to correct SNP counts from an alternative SNP list
 relative to a chosen genomic background

 IMPORTANT: We assume here that the alternative SNP are unique
 to the alternative genome and, therefore, all reference SNPs
 must have come from ther other genome

'''

# Get Arguments
for i in range(len(sys.argv)):
    if sys.argv[i] == "-filtered_SNPs":
        filteredSNP_file = sys.argv[i+1]

# Read lines
filteredSNP_file_lines = [ln.strip() for ln in open(filteredSNP_file,"r").readlines()]

# Reverse values in lines
outlines = []
for ln in filteredSNP_file_lines:
    #print ln
    #print ln.split("\t")
    [chr,pos,ref_base,snp_base,total,lost,ref_count,snp_count,ref_freq,snp_freq] = ln.split("\t")
    outlines.append("\t".join([chr,pos,snp_base,ref_base,total,lost,snp_count,ref_count,snp_freq,ref_freq]) + "\n")

output = open(filteredSNP_file + ".rev","w")
output.write("".join(outlines))
output.close()
