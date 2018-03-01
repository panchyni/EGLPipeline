# IMPORT
import sys

# MAIN
#print'''
'''
 Takes two filtered SNP files and combines them into a
 single sorted file

 Input:
      -fileA = 1st SNP file
      -fileB = 2nd SNP file

 Optionally, SNP labes can be specified as arguments
 using the following format, but if they are absent,
 this script will default to the SNP name in the file
 name (assuming the files are from "1a_FilterVCFSummaryFile.py"
 or "1b_ReverseSNPCounts.py") or to "A" and "B"

 Optional Input:
      -labelA = 1st SNP label
      -labelB = 2nd SNP label

 DevNote: 
  This script is based off of "2.0_MergeSNPFiles.py" in
  "/mnt/home/panchyni/1_SideProjects/17_ChlamyRIL_NEW/0.1_RevisedCode"

'''

labelA = ""
labelB = ""

# Get arguments
for i in range(len(sys.argv)):
    if sys.argv[i] == "-fileA":
        fileA = sys.argv[i+1]
    if sys.argv[i] == "-fileB":
        fileB = sys.argv[i+1]
    if sys.argv[i] == "-labelA":
        lableA = sys.argv[i+1]
    if sys.argv[i] == "-labelB":
        labelB = sys.argv[i+1]

# Check labels
if labelA == "":
   try:
      label_div = [div for div in fileA.split(".") if div.startswith("filtered")]
      labelA = label_div[0].split("_")[1]
   except IndexError:
      labelA = "A"

if labelB == "":
   try:
      label_div = [div for div in fileB.split(".") if div.startswith("filtered")]
      labelB = label_div[0].split("_")[1]
   except IndexError:
      labelB = "B"

# Read SNP data set
input_A_lines = [l.strip() for l in open(fileA,"r")]
input_B_lines = [l.strip() for l in open(fileB,"r")]

position_SNP_dict = {}

for ln in input_A_lines:
    if ln:
        split_ln = ln.strip().split("\t")
        position_SNP_dict[int(split_ln[1])] = split_ln

for ln in input_B_lines:
    if ln:
        split_ln = ln.strip().split("\t")
        position_SNP_dict[int(split_ln[1])] = split_ln

# Sort position and records outlines
keys = position_SNP_dict.keys()
keys.sort()

outlines = ["\t".join(position_SNP_dict[k])+"\n" for k in keys]

# Write ouptut
if fileA.endswith(".rev"):
    root = ".".join(fileA.split(".")[0:-2]) 
else:
    root = ".".join(fileA.split(".")[0:-1])
output = open(root + ".Combined","w")
output.write("".join(outlines))
output.close()


