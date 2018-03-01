# IMPORT
import sys
import bisect

# MAIN 
print'''

Reads a condensed region file (*.condensed.select_#####')[1] 
and a list of summarized VCF files [2]

Note, both of these file should be for the same chromsome/contig.

Returns a list SNPs for each region 

'''

# Read Inputs
region_lines = [ln.strip() for ln in open(sys.argv[1],"r").readlines()]
path = "/".join(sys.argv[2].split("/")[:-1])
vcf_files = [ln.strip() for ln in open(sys.argv[2],"r").readlines()]

# Get a lsit of regions
regions = []
for r_ln in region_lines:
    [region,value] = r_ln.split("\t")
    regions.append(region)

# Make SNP dict indexed by location
snp_dict = {}
for vf in vcf_files:
    vcf_lines = [ln for ln in open(path + "/" + vf,"r").readlines()]
    for v_ln in vcf_lines:
        snp_index = v_ln.split("\t")[1]
        snp_dict[int(snp_index)] = v_ln

# Make a list of snp_positions
snp_indexes = snp_dict.keys()
snp_indexes.sort()

# Bisect regions to find SNPs 
region_SNP_dict = {}
for r in regions:
    [start,stop] = r.split("_")
    start = int(start)
    stop = int(stop)

    start_index = bisect.bisect_left(snp_indexes,start)
    stop_index = bisect.bisect_right(snp_indexes,stop)

    # If the stop index is beyond the list
    region_snps = []
    if stop_index > len(snp_indexes) - 1:
        region_snps = [snp_dict[index] for index in snp_indexes[start_index:]] 
    else:
        region_snps = [snp_dict[index] for index in snp_indexes[start_index:stop_index]]

    region_SNP_dict[r] = region_snps

# Write output
outlines = []
for r in region_SNP_dict.keys():
    r_snps = region_SNP_dict[r]
    outlines.append(">" + r + "\n")
    outlines.extend(r_snps)

output = open(sys.argv[1].split("/")[-1] + "_BY_" + sys.argv[2].split("/")[-1] + ".snps","w")
output.write("".join(outlines))
output.close()
