# IMPORT
import sys

# MAIN
#print'''
'''
Read a Combined*AlleleFreq or a Combined*SNPCount file and reports
the percentage of covered (i.e. non-NA) window per columns

'''

input_lines = [ln.strip() for ln in open(sys.argv[1],"r").readlines()]
columns = input_lines[0].split("\t")[1:]

for i in range(len(columns)):
    values = [ln.split("\t")[i+1] for ln in input_lines]
    percentage = float(len(values) - values.count("NA"))/float(len(values))
    print columns[i] + ": " +str(percentage*100) + "% Coverage"
