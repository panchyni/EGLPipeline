# IMPORT
import sys

# MAIN
#print'''
'''
 Read a set of *.AlleleFreq files and combine into a common matrix

 Input: 
      -files_list = a file that is a list of files
      -outfile = root name of the outputfile

'''

# Get arguments
for i in range(len(sys.argv)):
    if sys.argv[i] == "-files_list":
        files_list = sys.argv[i+1]
    if sys.argv[i] == "-outfile":
        outfile = sys.argv[i+1]

# Get files to combine
files = [ln.strip() for ln in open(files_list,"r").readlines()]

# Read files and get windows
windows = []
files_dict = {}
for f in files:
    # Define a dictionary entry and window list for each file
    files_dict[f] = {}
    f_windows = []
    for ln in open(f,"r").readlines()[1:]:
        # Update dictionary and window list
        split_ln = ln.strip().split("\t")
        window = split_ln[0]
        f_windows.append(window)
        files_dict[f][window] = split_ln[1:]
    windows.extend(f_windows)

windows = list(set(windows))

# Sort our list of windows
pairs = [(int(w.split("_")[0]),int(w.split("_")[1])) for w in windows]
sorted_start = [A for (A,B) in sorted(pairs, key=lambda pair: pair[0])]
sorted_stop = [B for (A,B) in sorted(pairs, key=lambda pair: pair[0])]
sorted_windows = [str(sorted_start[i]) + "_" + str(sorted_stop[i]) for i in range(len(sorted_stop))]

# Get merge values across windows
header = "Window\t" + "\t".join(files) + "\n"
allele_lines = []
count_lines = []
# For each window
for w in sorted_windows:
    alleles = []
    counts = []
    # For each file
    for f in files:
        # If the given file has an entry for the current window
        if w in files_dict[f]:
            counts.append(files_dict[f][w][0])           
            alleles.append(files_dict[f][w][1])
        else:
            alleles.append("NA")
            counts.append("NA")
    allele_lines.append(w + "\t" + "\t".join(alleles) + "\n")
    count_lines.append(w + "\t" + "\t".join(counts) + "\n")

# Write output 
output = open(outfile + ".AlleleFreq","w")
output.write(header)
output.write("".join(allele_lines))
output.close()

output = open(outfile + ".SNPCounts","w")
output.write(header)
output.write("".join(count_lines))
output.close()
    
