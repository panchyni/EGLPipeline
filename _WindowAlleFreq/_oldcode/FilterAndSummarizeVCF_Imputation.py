# IMPORT
import sys

# FUNCTIONS

def VCFtoDict(vcf_file):
    '''

    Read a VCF file into a nested dictionary keyed first by chromsoeome and 
    then by position. Dictionary values are a list of ReferenceBase and SNPBase

    '''

    vcf_dict= {}
    count = 0
    for ln in open(vcf_file,"r").readlines():
        
        # If the line is not a header line
        if not ln.startswith("#"):

            # Split ln
            split_ln = ln.strip().split("\t")

            # Parse site information
            [chr,pos,id,ref,alt] = split_ln[0:5]
            
            info_string = split_ln[7].split(";")
            reads_split = [val.split("=")[1].split(",") for val in info_string if val.startswith("DP4=")][0] 
            ref_reads = int(reads_split[0]) + int(reads_split[1])       # Reads with Reference Base
            alt_reads = int(reads_split[2]) + int(reads_split[3])       # Reads with Alternate/SNP Base

            # If chromosome already in dictionary, add a new position entry
            if chr in vcf_dict.keys():
                vcf_dict[chr][pos] = [ref,alt,ref_reads,alt_reads]
            else:
                vcf_dict[chr] = {}
                vcf_dict[chr][pos] = [ref,alt,ref_reads,alt_reads]

            count = count + 1
            if count % 100000 == 0:
                print str(count/1000) + "k VCF lines read"

    return vcf_dict

def FilterVCF(vcf_dict,SNP_list,rel_strain):
    '''

    Takes a vcf_dict and retrives information based on a list of SNP positions
    Reads and allele frequency are reversed bsaed on the selected strain.
    Returns a dictionary of SNP data based on chromosome and position (for
    sorting) and a list of missing SNPs

    '''

    filtered_SNPs = []
    count = 0
    for snp in SNP_list:
        [chr,pos,strain] = snp[0:]
        
        try:

            # If strain is rel strain:
            if strain == rel_strain:
                [ref,alt,ref_reads,alt_reads] = vcf_dict[chr][pos]
            else:
                [alt,ref,alt_reads,ref_reads] = vcf_dict[chr][pos]
            
            filtered_SNPs.append([chr,pos,strain,ref,alt,ref_reads,alt_reads])
        except KeyError:
            
            # Infer the reference allele if no SNP
            
            if strain == rel_strain:
                [ref,alt,ref_reads,alt_reads] = ["*","*",100,0]
            else:
                [alt,ref,alt_reads,ref_reads] = ["*","*",100,0]
              
            filtered_SNPs.append([chr,pos,strain+"*",ref,alt,ref_reads,alt_reads])

        count = count + 1
        if count % 100000 == 0:
            print str(count/1000) + "k SNP position processed"

    # Make a dictionry of the filtered SNPs
    filtered_SNP_dict = {}
    for snp in filtered_SNPs:
        [chr,pos,strain,ref,alt,ref_reads,alt_reads] = snp
         
        if chr in filtered_SNP_dict.keys():
            filtered_SNP_dict[chr][pos] = [strain,ref,alt,ref_reads,alt_reads]
        else:
            filtered_SNP_dict[chr] = {}
            filtered_SNP_dict[chr][pos] = [strain,ref,alt,ref_reads,alt_reads]

    return filtered_SNP_dict

def SortSNPDictionary(snp_dict):
    '''

    Takes a nested dictionary of SNP information where the first level keys are
    chromosome/contig names and the second level of keys are the position of
    the SNP on the chromosome. Chromosome/contigs should have the same based name
    followed by a number seperated by "_". Return as sorted list of SNP information.

    '''

    # Sort contigs
    contigs = snp_dict.keys()
    contig_pairs = [(v.split("_")[0],int(v.split("_")[1])) for v in contigs]
    contig_pairs.sort(key=lambda x: x[1])
    contigs_sorted = [v[0] + "_" + str(v[1]) for v in contig_pairs]

    # Sort position within each contig
    snp_list = []
    for contig in contigs_sorted:

        snp_positions = [int(v) for v in snp_dict[contig].keys()]
        snp_positions.sort()
        snp_positions = [str(v) for v in snp_positions]

        for position in snp_positions:
            snp_info = [contig,position]
            snp_info.extend(snp_dict[contig][position])
            snp_list.append(snp_info)

    return snp_list

# MAIN
def main():
    print'''

    Converts a full VCF file to filtered summary using a list of SNP positions
    annotated with strain of origin and a strain interest such that all other
    strains are treated as the background.

    Inputs:
        -vcf 	= VCF file to be filtered
        -snp 	= SNP list for filtering
        -strain = strain to calculate SNP frequency relaitve to


    Ouputs:
        A file ending in *.VCF_summary that contains the position, base, and allele
        frequency for each filtered SNP

    '''
 
    # Read inputs arguments
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-vcf":
            vcf_file = sys.argv[i+1]
        if sys.argv[i] == "-snp":
            SNP_list_file = sys.argv[i+1]
        if sys.argv[i] == "-strain":
            rel_strain = sys.argv[i+1]

    # Read VCF file into dictionary
    vcf_dict = VCFtoDict(vcf_file)

    # Filter dictionary by SNP list
    SNP_list = [ln.strip().split("\t") for ln in open(SNP_list_file,"r").readlines()] 
    filtered_SNP_dict = FilterVCF(vcf_dict,SNP_list,rel_strain)
    
    # Sort filtered SNPs
    filtered_SNPs = SortSNPDictionary(filtered_SNP_dict)

    # Calculate Allele Frequency 
    filter_outlines = []
    for snp in filtered_SNPs:
        [ref_reads,alt_reads] = snp[-2:]
        snp_allele_freq = float(alt_reads)/float(ref_reads+alt_reads)
        filter_outlines.append("\t".join([str(v) for v in snp]) + "\t" + str(snp_allele_freq) + "\n")

    # Write output
    filter_outfile = open(vcf_file + ".VCF_summary","w")
    filter_outfile.write("".join(filter_outlines))
    filter_outfile.close()

# If not loaded as a module
if __name__ == '__main__':
    main()
