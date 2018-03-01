# IMPORT
import sys

# MAIN
# Reads a standard *.vcf file [1] and write a simplified ".reads.summary" file

# Reads vcf file 
vcf_lines = [line for line in open(sys.argv[1],"r") if not line.startswith("#")]

# Record outliens
outlines = []
header = "Chrom" + "\t" +"Position" + "\t" + "RefBase" + "\t" + "AltBase" + "\t" + "TotalReads" + "\t" + "DroppedReads" + "\t" + "RefReads" + "\t" + "AltReads" + "\t" + "RefFreq" + "\t" + "AltSeq" + "\n"
for l in vcf_lines:

    split_l = l.strip().split("\t")
    chrom = split_l[0]		# Chromosome
    pos = split_l[1]		# Position
    ref_base = split_l[3] 	# Refenrence Base
    alt_base = split_l[4]	# Alternate/SNP Base
    info_string = split_l[7].split(";")		# Info string about this SNP site
    total_reads = int([val.split("=")[1] for val in info_string if val.startswith("DP=")][0])        # Parse total reads, [0] use 0 b/c we want a value not list
    reads_split = [val.split("=")[1].split(",") for val in info_string if val.startswith("DP4=")][0] # Get read distrib , ditto
    af1 = [val.split("=")[1] for val in info_string if val.startswith("AF1=")]
    ref_reads = int(reads_split[0]) + int(reads_split[1])	# Reads with Reference Base
    alt_reads = int(reads_split[2]) + int(reads_split[3])	# Reads with Alternate/SNP Base
    dropped_reads = total_reads - (ref_reads + alt_reads)	# Dropped Reads = Total Reads - Sum(Ref Reads, Alt Reads)
    ref_freq = float(ref_reads)/float(ref_reads+alt_reads)	# Frequecny of Ref Reads = Ref/(Ref+Alt)
    alt_freq = float(alt_reads)/float(ref_reads+alt_reads)	# Frequecny of Alt Reads = Alt/(Ref+Alt)
    outlines.append(chrom + "\t" + pos + "\t" + ref_base + "\t" + alt_base + "\t" + str(total_reads) + "\t" + str(dropped_reads) + "\t" + str(ref_reads) + "\t" + str(alt_reads) + "\t" + str(ref_freq) + "\t" + str(alt_freq) + "\n")  

# Write output
output = open(sys.argv[1].split("/")[-1]+".reads_summary","w")
output.write(header)
output.write("".join(outlines))
output.close()
