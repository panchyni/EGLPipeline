# IMPORT 
import sys

# MAIN
#print'''
'''
 Reads a FILTERED "*.reads_summary" and calculated the allele frequency over
 a window of a specified size

 IMPORTANT: The input file must either be prefiltered for SNPs for a certain
 background or all SNPs need to be corrected for a certain background. Each
 file should contain all of the SNPs for a single chromsome/contig
 
 Input:
      -SNP_file = input, filtered SNP file
      -wsize = window size 
      -witer = step size between windows
      -minSNP = the minimum number of SNPs per window not to be called "NA"

 DevNote:
  This script is derived from "1.4_WindowedInocAlleleFreq.py" in 
  "/mnt/home/panchyni/1_SideProjects/17_ChlamyRIL_NEW/0.1_RevisedCode"
 
'''

# Get arguments
for i in range(len(sys.argv)):
    if sys.argv[i] == "-SNP_file":
        SNP_file = sys.argv[i+1]
    if sys.argv[i] == "-wsize":
        wsize = int(sys.argv[i+1])
    if sys.argv[i] == "-witer":
        witer = int(sys.argv[i+1])
    if sys.argv[i] == "-minSNP":
        min_SNP = int(sys.argv[i+1])

# Read file
SNP_lines = [ln.strip() for ln in open(SNP_file,"r").readlines()]

# Make a positional dictioanry of the SNPs
position_dict = {}
for ln in SNP_lines:
    split_ln = ln.strip().split("\t")
    position_dict[int(split_ln[1])] = split_ln[4:]

position_keys = position_dict.keys()
max_key = max(position_keys)

# Make a list of windows based on wsize, witer and the length of contig (max key)
window_end_list = range(wsize,max_key+witer,witer)

outlines = []
for w in window_end_list:

    # Get all the SNP positions in the current window
    w_position_keys = [i for i in position_keys if i < w and i >= w-wsize]
    window_id = str(w-wsize) + "_" + str(w)
    
    # Check the SNP number and proceed if it passes the threshold
    SNP_number = len(w_position_keys)
    if SNP_number > min_SNP:
        
        # Calculate the SNP frequency across the window
        SNP_Freq = [float(position_dict[k][-1]) for k in w_position_keys]
        Avg_SNP_Freq = sum(SNP_Freq)/len(SNP_Freq)
        outlines.append(window_id + "\t" + str(SNP_number) + "\t" + str(Avg_SNP_Freq) + "\n")

    else:
        outlines.append(window_id + "\t" + str(SNP_number) + "\tNA\n")

header = "Window\tSNP\tFreq\n"
output = open(SNP_file + ".Window_" + str(wsize) + "_" + str(witer) + ".AlleleFreq","w")
output.write(header)
output.write("".join(outlines))
output.close()
