# IMPORT
import sys

# FUNCTIONS

def SummarytoDict(snp_file):
    '''
     
    Reads a summary snp file into a dictioanry keyed 
    by chromosome then position

    '''

    snp_dict= {}
    count = 0
    for ln in open(snp_file,"r").readlines():

        # Split ln
        split_ln = ln.strip().split("\t")

        # Parse site information
        [chr,pos] = split_ln[0:2]
        site_info = split_ln[2:]


        # If chromosome already in dictionary, add a new position entry
        if chr in snp_dict.keys():
            snp_dict[chr][pos] = site_info
        else:
            snp_dict[chr] = {}
            snp_dict[chr][pos] = site_info

        count = count + 1
        if count % 100000 == 0:
            print str(count/1000) + "k VCF lines read"

    return snp_dict

def SortContigs(contigs):
    '''

    Read a list of contigs and sorts them in numeric order .Chromosome/contigs 
    should have the same based name followed by a number seperated by "_".

    '''

    contig_pairs = [(v.split("_")[0],int(v.split("_")[1])) for v in contigs]
    contig_pairs.sort(key=lambda x: x[1])
    contigs_sorted = [v[0] + "_" + str(v[1]) for v in contig_pairs]

    return contigs_sorted

def DefineWindows(snp_positions,window_size,window_iter):
    '''

    Divide a chromosome into windows of given size with a given step
    between windows

    '''

    # Find the end of the chromosome
    max_position = max(snp_positions)

    # Define windows by their ends
    window_end_list = range(window_size,max_position+window_iter,window_iter)

    return window_end_list


def GetWindowedAlleleFreq(snp_dict,window_size,window_iter,min_SNP):
    '''

    Calculate the average allele frequency of SNPs across the a chromsome
    over windows of defined size (window_size) and step (window_iter)

    '''

    snp_positions = [int(v) for v in snp_dict.keys()]
    window_end_list = DefineWindows(snp_positions,window_size,window_iter)
    
    # For each window
    windows = []
    count = 0
    for w in window_end_list:
        w_position_keys = [i for i in snp_positions if i < w and i >= w-window_size]
        window_id = str(w-window_size) + "_" + str(w)
        SNP_number = len(w_position_keys)

        # Check if there is the necessary number of SNP per window
        if SNP_number >= min_SNP:

            SNP_Freq = []
            for k in w_position_keys:
                SNP_Freq.append(float(snp_dict[str(k)][-1]))

            Avg_SNP_Freq = sum(SNP_Freq)/len(SNP_Freq)

            windows.append([window_id, str(SNP_number), str(Avg_SNP_Freq)])

        else:
            windows.append([window_id, str(SNP_number), "NA"])

        count = count + 1
        if count % 1000 == 0:
            print str(count/1000) + "k windows processed for current contig"

    return windows

# MAIN
def main():
    print'''

    Reads a summarized SNP file with chromosome in the 
    first column, position in the second, and allele 
    frequency in the last column.

    Calculates the windowed allele frequency of windows
    of a specified size and iteration.

    Input:
        -snp = name of the summarized SNP file
        -wsize = window size
        -witer = window size iteration
        -minSNP = the minimum number of SNPs per window

    '''

    for i in range(len(sys.argv)):
        if sys.argv[i] == "-snp":
            snp_file = sys.argv[i+1]
        if sys.argv[i] == "-wsize":
            window_size = int(sys.argv[i+1])
        if sys.argv[i] == "-witer":
            window_iter = int(sys.argv[i+1])
        if sys.argv[i] == "-minSNP":
            min_SNP = int(sys.argv[i+1])

    # Read SNP data into a dictionary
    snp_dict = SummarytoDict(snp_file)
    
    # Sort chromosomes
    sorted_chromosomes = SortContigs(snp_dict.keys())

    # For each chromosome
    outlines = []
    for chr in sorted_chromosomes:
        
        window_allele_freq = GetWindowedAlleleFreq(snp_dict[chr],window_size,window_iter,min_SNP)

        for w in window_allele_freq:
            [id,count,freq] = w
            [start,stop] = id.split("_")
            outlines.append(chr + "\t" + start + "\t" + stop + "\t" + count + "\t" + freq + "\n")

    # Write output
    output = open(snp_file + ".WindowedAlleFreq","w")
    header = "Chr\tWStart\tWStop\tSNPCount\tAlleFreq\n"
    output.write(header)
    output.write("".join(outlines))
    output.close()

# If not loaded as a module
if __name__ == '__main__':
    main()


