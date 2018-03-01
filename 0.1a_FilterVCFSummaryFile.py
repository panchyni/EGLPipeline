# IMPORT
import sys

# MAIN
print'''

 Filters a condensed VCF ".reads_summary" file, filters
 its using a 2 column SNP list and annotated the filtered
 file with a SNP ID

 Input:
      -summary = the ".reads_summary" file derived from a VCF file
      -SNPlist = the 2 col SNP list file
      -SNPid = a string used to ID the SNPs

 DevNote:
  This script is derived from "0.2a_FilterDataAnnotSNP.py" in
  "/mnt/home/panchyni/1_SideProjects/17_ChlamyRIL_NEW/0.1_RevisedCode"

'''

# Get arguments
for i in range(len(sys.argv)):
    if sys.argv[i] == "-summary":
        summary_file = sys.argv[i+1]
    if sys.argv[i] == "-SNPlist":
        SNP_list_file = sys.argv[i+1]
    if sys.argv[i] == "-SNPid":
        SNP_id = sys.argv[i+1]

# Get read files
input_lines = open(summary_file,"r").readlines()
SNP_list = [ln.strip().split("\t")[0] + "-" + ln.strip().split("\t")[1] for ln in open(SNP_list_file,"r").readlines()]

# Filter input lines
filtered_lines = [ln for ln in input_lines if ln.strip().split("\t")[0] + "-" + ln.strip().split("\t")[1] in SNP_list]

# Write ouput
output = open(summary_file + ".filtered_" + SNP_id, "w")
output.write("".join(filtered_lines))
output.close()
