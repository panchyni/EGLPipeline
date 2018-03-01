# IMPORT 
import sys
import time
import bisect
import OverlapFunctions

# FUNCTIONS

def FormatTime(seconds):
    # Return a formatted string of the current time from the number of seconds that have passed since starting the program
    # Input
    #     seconds := number of seconds (use time.time() - start_imter)
    # Output
    #     current_time := formated string of time since program start
    m,s = divmod(seconds,60)
    h,m = divmod(m,60)
    current_time = "%d:%02d:%02d" % (h, m, s)
    return current_time

def GetFromGFFString(string,tag):
    # Takes a text object (string) representing the string field from
    # a GFF files and another text object (tag) indicating which value
    # to extract from the string
    try:
        substring = gene_string.split(",")[-1].split(tag+"=")[-1].split(";")[0]
    except:
        print "GFF string does not contain an " + tag + " tag"
        substring = "Missing"
    return substring

# ERROR CLASSES
class Error(Exception):
      """ Base error class """
      pass

class MissingParameter(Error):
      """ Exception raised for missing paramters

      Attributes:
         -message: The error message to be dispalyed

      """
      def __init__(self,message):
          self.message = message
      def __str__(self):
          return repr(self.message)

class UnknownValue(Error):
      """ Exception raised for missing paramters

      Attributes:
         -message: The error message to be dispalyed
         -value: The value of the variable which caused the error

      """
      def __init__(self,message):
          self.message = message
      def __str__(self):
          return repr(self.message + ": " + str(value))

# MAIN
print'''

Wrapper for mapping SNPs to region of significant allele frequency
change and to specific coding sequences within that region.

Takes a control file [1, see "SNPMappingTeplcate.ctl"] which
specifies all the parameters need for this pipeline

This pipeline involves the following steps:
 -(1) Mapping SNPs to a given region of the genome 
 -(2) Finding overlap between SNPs positio nand features (Genes,CDS,Exons) in that region

Future extension involve:
 -(3) Mapping SNPs changes to CDS sequences
 -(4) Determine putative change to protein sequence cause by SNP

'''

# Start timer
start_timer = time.time()
print "Read Inputs [" + FormatTime(time.time() - start_timer) + "]"


# Read Input
control_file = [ln.strip() for ln in open(sys.argv[1],"r").readlines()]

# Read parameters
parameter_dictionary = {}
for index in range(len(control_file)):
    param_ln = control_file[index]
    try:
        [name,value] = param_ln.split("\t")
        parameter_dictionary[name] = value
    except ValueError as e:
       e.args += ("Control file: line " + str(index+1),)
       raise e

# Check required parameters
required_parameters = ["RegionsFile","SNPFiles","SNPFilesPath","GFFfile","CDSFile"]

for param in required_parameters:
    if not param in parameter_dictionary.keys():
        raise MissingParameter("Control file is missing a line for parameter: " + param)
    else:
        print "Parameter " + param + " is present and valued"

print "All Inputs present and accounted for [" + FormatTime(time.time() - start_timer) + "]"

## Run Pipeline ##

# Step 1: Obtain SNPs for region #
#

print  "Start Step1: Obtain SNPs per region [" + FormatTime(time.time() - start_timer) + "]"

# Read input from files
region_lines = [ln.strip() for ln in open(parameter_dictionary["RegionsFile"],"r").readlines()]
snp_files_list = [ln.strip() for ln in open(parameter_dictionary["SNPFiles"],"r").readlines()]
snp_files_path = parameter_dictionary["SNPFilesPath"]
if not snp_files_path.endswith("/"):
   snp_files_path = snp_files_path + "/"

# Get a list of regions 
regions = []
for ln in region_lines:
    [region,value] = ln.split("\t")
    regions.append(region)

# Make SNP dict indexed by location
snp_dict = {}
for file in snp_files_list:
    snp_lines = [ln for ln in open(snp_files_path + file,"r").readlines()]
    for ln in snp_lines:
        snp_index = ln.split("\t")[1]
        snp_dict[int(snp_index)] = ln

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

print "Step1 Done [" + FormatTime(time.time() - start_timer) + "]"

# Step 2: Map SNPs to features in region #
#

print "Start Step2: Map SNPs to features in region [" + FormatTime(time.time() - start_timer) + "]"

# Read GFF file
gff_lines = [ln.strip() for ln in open(parameter_dictionary["GFFfile"],"r").readlines()]

# Filter GFF file by chromosome
snp_chr = region_SNP_dict[regions[0]][0].split("\t")[0]
gff_filtered = [ln for ln in gff_lines if ln.startswith(snp_chr)]

# Make a list of regions 
regions_list = [snp_chr + "\t" + "\t".join(key.split("_")) for key in region_SNP_dict.keys()]

# Subset GFF file
gene_lines = [ln.strip() for ln in gff_filtered if ln.split("\t")[2] == "gene"]
CDS_lines = [ln.strip() for ln in gff_filtered if ln.split("\t")[2] == "CDS"]
exon_lines = [ln.strip() for ln in gff_filtered if ln.split("\t")[2] == "exon"]

# Get Features in each regions
[GeneOverlaps,region_keys] = OverlapFunctions.FindOverlapsByIndexing([regions_list,0,1,2],[gene_lines,0,3,4],10000)
[CDSOverlaps,region_keys] = OverlapFunctions.FindOverlapsByIndexing([regions_list,0,1,2],[CDS_lines,0,3,4],10000)
[ExonOverlaps,region_keys] = OverlapFunctions.FindOverlapsByIndexing([regions_list,0,1,2],[exon_lines,0,3,4],10000)

# Get SNPs for each feature
gene_summary_out = []
cds_summary_out = []
exon_summary_out = []
gene_snps_out = []
cds_snps_out = []
exon_snps_out = []
for key in region_keys:
    region_genes = [snp_chr + "\t" + gene[0] + "\t" + gene[1] + "\t" + gene[-1] for gene in GeneOverlaps[key][1:-1]]
    region_CDS = [snp_chr + "\t" + cds[0] + "\t" + cds[1] + "\t" + cds[-1] for cds in CDSOverlaps[key][1:-1]]
    region_exons = [snp_chr + "\t" + exon[0] + "\t" + exon[1] + "\t" + exon[-1] for exon in ExonOverlaps[key][1:-1]]
    region_snps = region_SNP_dict[str(key[1]) + "_" + str(key[2])]
    
    # Sort features
    region_genes.sort()
    region_CDS.sort()
    region_exons.sort()

    # Find SNPs
    [GeneSNPs,GeneSNPKeys] = OverlapFunctions.FindOverlapByBisection([region_genes,0,1,2],[region_snps,0,1])
    [CDSSNPs,CDSSNPKeys] = OverlapFunctions.FindOverlapByBisection([region_CDS,0,1,2],[region_snps,0,1])
    [ExonSNPs,ExonSNPKeys] = OverlapFunctions.FindOverlapByBisection([region_exons,0,1,2],[region_snps,0,1])

    # Summarize each feature for output
    region_name = str(key[1]) + "_" + str(key[2])
    for gene_key in GeneSNPKeys:
        gene_string = GeneSNPs[gene_key][0]
        gene_name = GetFromGFFString(gene_string,"ID")
 
        # Write summary   
        gene_summary_out.append(region_name + "\t" + "\t".join([str(v) for v in gene_key]) + "\t" + gene_name + "\t" + str(GeneSNPs[gene_key][-1]) + "\n")
           
        # Write snps
        gene_chr = str(gene_key[0])
        gene_snps = GeneSNPs[gene_key][1:-1]
        gene_snps_out.append(gene_name+"\n") 
        for s in gene_snps:
            gene_snps_out.append(gene_chr + "\t" + s[0] + "\t" + "\t".join(s[1].split(",")))

    for cds_key in CDSSNPKeys:
        cds_string = CDSSNPs[cds_key][0]
        cds_name = GetFromGFFString(cds_string,"ID")

        # Write summary        
        cds_summary_out.append(region_name + "\t" + "\t".join([str(v) for v in cds_key]) + "\t" +  cds_name + "\t" + str(CDSSNPs[cds_key][-1]) + "\n")

        # Write snps
        cds_chr = str(cds_key[0])
        cds_snps = CDSSNPs[cds_key][1:-1]
        cds_snps_out.append(cds_name+"\n")
        for s in cds_snps:
            cds_snps_out.append(cds_chr + "\t" + s[0] + "\t" + "\t".join(s[1].split(",")))

    for exon_key in ExonSNPKeys:
        exon_string = ExonSNPs[exon_key][0]
        exon_name = GetFromGFFString(exon_string,"ID")
        exon_summary_out.append(region_name + "\t" + "\t".join([str(v) for v in exon_key]) + "\t" +  exon_name + "\t" + str(ExonSNPs[exon_key][-1]) + "\n")

        # Write summary
        exon_chr = str(exon_key[0])
        exon_snps = ExonSNPs[exon_key][1:-1]
        exon_snps_out.append(exon_name+"\n")
        for s in exon_snps:
            exon_snps_out.append(exon_chr + "\t" + s[0] + "\t" + "\t".join(s[1].split(",")))

# Write summary of each feature
header = "Region\tChromosome\tStart\tStop\tFeature#\tFeatureID\tSNPCount\n" 
output_gene = open(sys.argv[1] + ".gene_summary","w")
output_gene.write(header)
output_gene.write("".join(gene_summary_out))
output_gene.close()

output_cds = open(sys.argv[1] + ".cds_summary","w")
output_cds.write(header)
output_cds.write("".join(cds_summary_out))
output_cds.close()

output_exon = open(sys.argv[1] + ".exon_summary","w")
output_exon.write(header)
output_exon.write("".join(exon_summary_out))
output_exon.close()

# Write list of snps for each feature
output_gene_snps = open(sys.argv[1] + ".gene_snps","w")
output_gene_snps.write("".join(gene_snps_out))
output_gene_snps.close()

output_cds_snps = open(sys.argv[1] + ".cds_snps","w")
output_cds_snps.write("".join(cds_snps_out))
output_cds_snps.close()

output_exon_snps = open(sys.argv[1] + ".exon_snps","w")
output_exon_snps.write("".join(exon_snps_out))
output_exon_snps.close()

# TEST CODE
#print GeneSNPs[GeneSNPs.keys()[5]]
#print CDSSNPs[CDSSNPs.keys()[5]]
#print ExonSNPs[ExonSNPs.keys()[5]]

print "Step2 Done [" + FormatTime(time.time() - start_timer) + "]"

# Step 3: Map SNPs to CDS Sequence#
#

print "Start Step3: Map SNPS to CDS Sequence [" + FormatTime(time.time() - start_timer) + "]"

# Read CDS Sequence file
cds_lines = [ln.strip() for ln in open(parameter_dictionary["CDSFile"],"r").readlines()]

cds_dict = {}
gene = "start"
seq = "seq"
for ln in cds_lines:

    if ln.startswith(">"):
        cds_dict[gene] = seq
        gene = ln.split(" ")[0].strip(">")
        seq = ""
    else:
        seq = seq + ln

del cds_dict["start"]

# Test CODE
#print cds_dict[cds_dict.keys()[0]]
