# IMPORT
import sys
import OverlapFunctions

# MAIN
print'''

Takes a *snps file of region snps [1] and a gff file
of genes [2].

The script first finds the genes that overlap with
each region. Then, for the SNPs in that region, 
finds overlaps with a gene.

Returns a FASTA style file of SNPs in each gene


'''
snp_lines = [ln.strip() for ln in open(sys.argv[1],"r").readlines()]
gff_lines = [ln.strip() for ln in open(sys.argv[2],"r").readlines()]
snp_chr = snp_lines[1].split("\t")[0]
gene_lines = [ln for ln in gff_lines if ln.startswith(snp_chr)]

# Filter genes by SNP chromsome
filtered_genes = [ln for ln in gff_lines if ln.startswith(snp_chr)]

# Make a dictioanry of SNPs by regions
snp_dict = {}
region = "start"
snps = []
for ln in snp_lines:
    if ln.startswith(">"):
        snp_dict[region] = snps
        region = snp_chr + "\t" + "\t".join(ln.strip(">").split("_"))
        snps = []
    else:
        snps.append(ln)

del snp_dict["start"]
snp_region_lines = snp_dict.keys()

# Get genes in each reigon
[GFFRegionOverlaps,sorted_keys] = OverlapFunctions.FindOverlapsByIndexing([snp_region_lines,0,1,2],[gene_lines,0,3,4],10000)

# Get SNPS for each gene
outlines = []
for key in sorted_keys:
    region_genes = [snp_chr + "\t" + gene[0] + "\t" + gene[1] + "\t" + gene[-1] for gene in  GFFRegionOverlaps[key][1:-1]]
    region_snps = snp_dict["\t".join(key[0:3])] 

    # Sort Genes
    region_genes.sort()

    [GeneSNPOverlap,GeneSNPSortedKeys] = OverlapFunctions.FindOverlapByBisection([region_genes,0,1,2],[region_snps,0,1])
    
    for gene in GeneSNPSortedKeys:
        gene_info = GeneSNPOverlap[gene][0]
        snps = GeneSNPOverlap[gene][1:-1]
       
        outlines.append(">" + gene_info.split("=")[-1] + "\n") 
        
        for s in snps:
            outlines.append(snp_chr + "\t" + s[0] + "\t" + "\t".join(s[1].split(",")) + "\n")

# Write output
     
output = open(sys.argv[1] + ".genes","w")
output.write("".join(outlines))
output.close()
