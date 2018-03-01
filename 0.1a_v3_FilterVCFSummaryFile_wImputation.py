# IMPORT
import sys

# MAIN
#print'''
'''
 Filters a condensed VCF ".reads_summary" file, filters
 its using a 2 column SNP list and imputes the frequency
 of missing SNPs sites as being from the reference genome.
 Annotates the filtered file with a SNP ID

 Input:
      -summary = the ".reads_summary" file derived from a VCF file
      -SNPlist = the 2 col SNP list file
      -SNPid = a string used to ID the SNPs

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

# Read intpt file into a dictionary and record positions
filter_dict = {}
for ln in input_lines:
    position = ln.strip().split("\t")[0] + "-" + ln.strip().split("\t")[1]
    filter_dict[position] = ln

# For every SNP postion, either retrive or impute
filtered_lines = []
for SNP_position in SNP_list:
    if SNP_position in filter_dict.keys():
        filtered_lines.append(filter_dict[SNP_position])
    else:
       imputed_site = SNP_position.split("-")[0] + "\t" + SNP_position.split("-")[1] + "\t" + "N\tN\tX\tX\tX\tX\t1.0\t0.0\n"
       filtered_lines.append(imputed_site)

# Write ouput
output = open(summary_file + ".filtered_imputed_" + SNP_id, "w")
output.write("".join(filtered_lines))
output.close()
