# IMPORT 
import sys
import os
import time

# FUNCTIONS

# Utility
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

# Step 1 #
def RunSNPFilter(read_summary,SNP,Scripts):
    ''' Runs "0.1a_v2_FilterVCFSummaryFile.py" with the specified files and parameters 

	Returns the name of the outputfile
    '''
    
    [SNP_file,SNP_id] = SNP.split(",")
    command = "python " + Scripts + "0.1a_v2_FilterVCFSummaryFile.py -summary " + read_summary + " -SNPlist " + SNP_file + " -SNPid " + SNP_id + "\n"
    os.system(command)
    return read_summary + ".filtered_" + SNP_id 

def RunReverseCounts(filtered_file,Scripts):
    ''' Runs "0.1b_ReverseSNPCounts.py" on the speicified files

    	Return the name of the outputfile
    '''
    
    command = "python " + Scripts + "0.1b_ReverseSNPCounts.py -filtered_SNPs " + filtered_file + "\n"
    os.system(command)
    return filtered_file + ".rev"

def MergeFiles(fileA,fileB,SNP_A,SNP_B,Scripts):
    ''' Runs the "0.1c_MergeFilteredSNPs.py" with the specified files and paramters

    '''
    [SNP_fileA,SNP_idA] = SNP_A.split(",")
    [SNP_fileB,SNP_idB] = SNP_B.split(",")
    if fileA.endswith(".rev"):
        root = ".".join(fileA.split(".")[0:-2])
    else:
        root = ".".join(fileA.split(".")[0:-1])

    command = "python " + Scripts + "0.1c_MergeFilteredSNPs.py -fileA " + fileA + " -fileB " + fileB + " -labelA " + SNP_idA + " -labelB " + SNP_idB + "\n"
    os.system(command)
    return root + ".Combined"

def SummarizeStep1(summary_file,BckSNP,AltSNP,Scripts):
    ''' Wrapper functions for all functions in Step1

        Return the final file generate by this step
    '''

    BckFilter = RunSNPFilter(summary_file,BckSNP,Scripts)
    AltFilter = RunSNPFilter(summary_file,AltSNP,Scripts)
    RevAltFilter = RunReverseCounts(AltFilter,Scripts)
    MergedFile = MergeFiles(BckFilter,RevAltFilter,BckSNP,AltSNP,Scripts)
    os.system("mv " + MergedFile + " " + MergedFile + ".Rel_" + BckSNP.split(",")[1])
 
    return MergedFile + ".Rel_" + BckSNP.split(",")[1]

# Step 2 #
def RunWindowedAlleleFreq(file,wsize,witer,minSNP,Scripts):
    ''' Runs the "0.2a_CalculateWindowedAlleleFrequency.py" function for all files from Step1

        Returns the files generate by this step
    '''

    command = "python " + Scripts + "0.2a_CalculateWindowedAlleleFrequency.py -SNP_file " + file + " -wsize " + wsize + " -witer " + witer + " -minSNP " + minSNP + "\n"
    os.system(command)

    outfile = file + ".Window_" + str(wsize) + "_" + str(witer) + ".AlleleFreq"

    if os.path.isfile(outfile):
        return outfile
    else:
         print "Error processing " + file + "\nCotinuing without this file."
    
def RunCombineAlleFreq(files,outfile,Scripts):
    '''Run the "0.2b_CombineAlleleFreqFile.py" function for all ".AlleleFreq" files
     
       Returns the name of the combined file
    '''

    command = "python " + Scripts + "0.2b_CombineAlleleFreqFile.py -files_list " + files + " -outfile " + outfile + "\n"
    os.system(command)
    
    return [outfile + ".AlleleFreq",outfile + ".SNPCounts"]

def RunCoverageCheck(allele_file,Scripts):
    ''' Run "0.2c_CheckCoverage.py" function on the combined AlleleFreq file from "2b_CombineAlleleFreqFile.py"

        Print converage statistic for each file in the matrix
    '''

    command = "python " + Scripts + "0.2c_CheckCoverage.py " + allele_file
    os.system(command)

# Step 3#
def ListMeanStdv(values):
    ''' Calculate the mean of a list of float values

        Return the mean value
    '''
    mean = sum(values)/len(values)
    divs = [(mean - i)**2 for i in values]
    stdv = (sum(divs)/len(divs))**(0.5)
    return mean, stdv
    
def InocDifference(inoc1,inoc2,AlleleMatrix):
    ''' Calculates the allele frequency difference between two Inoculum file

        Inputs:
              -inoc1 = name of the first Inoculum file 
              -inoc2 = name of the second Inoculum file
              -AlleleMatrix = name of the Allele Matrix file

        Returns the mean and standard deviation of the distribution of differences
    '''

    # Read Inoculum data for the AlleleMatrix file
    AlleleMatrixLines = [ln.strip() for ln in open(AlleleMatrix,"r").readlines()]
    in_header = AlleleMatrixLines[0].split("\t")
    inoc1_index = in_header.index(inoc1)
    inoc2_index = in_header.index(inoc2)
    inoc1_values = [ln.split("\t")[inoc1_index] for ln in AlleleMatrixLines]
    inoc2_values = [ln.split("\t")[inoc2_index] for ln in AlleleMatrixLines]
    
    # Take difference between Inoclumn
    diffs = []
    for index in range(len(inoc1_values[1:])):
        [inoc1,inco2] = [inoc1_values[index+1],inoc2_values[index+1]]
    
        if not inoc1 == "NA" and not inco2 == "NA":
            diffs.append(float(inoc1) - float(inco2))

    [mean,stdv] = ListMeanStdv(diffs)

    return abs(mean),stdv

def FoldedNormalPvalue(diff,mean,stdv):
    # Determine the p-value for derving the
    # value diff from a folded normal distribution
    # with parameters mean and stdv

    import numpy
    from scipy import stats
    shape = abs(mean)/stdv
    pvalue = 1 - (stats.foldnorm.cdf(diff,shape,loc=0,scale=stdv))

    # If p-value = 0, resent to 2.22e-16, which seems to be the 
    # smallest p-value you can generate before it goes to zero
    if pvalue == 0:
        pvalue = 2.22e-16
    return pvalue


def ExprDifference(expr_names,inoc_names,AlleleMatrix,mean,stdv,outfile):
    ''' Calculate the significance of allele frequency differences between Inoculum and Experimental data

        Inputs:
              -expr_names = names of the Experiment files
              -inoc_names = names of the Inoculum files
              -AlleleMatrix = name of the Allele Matrix file
              -mean = mean of the Inoculum difference distribution
              -stdv = stdv of the Inoculum difference distribution
              -outfile = root of the output file
        
        Return the name of the output matrix
    '''

    # Read Inoculum data for the AlleleMatrix file
    AlleleMatrixLines = [ln.strip() for ln in open(AlleleMatrix,"r").readlines()]
    in_header = AlleleMatrixLines[0].split("\t")
    windows = [ln.split("\t")[0] for ln in AlleleMatrixLines[1:]]

    # Get indexes
    expr_indexes = [in_header.index(expr) for expr in expr_names]
    inoc_indexes = [in_header.index(inoc) for inoc in inoc_names]

    # Make output dictionaries
    diff_dict = {}
    sig_dict = {}
    for w in windows:
        diff_dict[w] = []
        sig_dict[w] = []

    # For each combination of Experiment and Inoculum
    out_header = []
    for expr_i in expr_indexes:
        for inoc_i in inoc_indexes:
            out_header.append(in_header[expr_i] + "_vs_" + in_header[inoc_i])
            expr_values = [ln.split("\t")[expr_i] for ln in AlleleMatrixLines[1:]]
            inoc_values = [ln.split("\t")[inoc_i] for ln in AlleleMatrixLines[1:]]
                    
            # For each window
            for i in range(len(windows)):
                e_v = expr_values[i]
                i_v = inoc_values[i]

                if not e_v == "NA" and not i_v == "NA":
                    diff = float(e_v) - float(i_v)
                    pvalue = FoldedNormalPvalue(abs(diff),mean,stdv)
                    diff_dict[windows[i]].append(diff)
                    sig_dict[windows[i]].append(pvalue)
                else:
                    diff_dict[windows[i]].append("NA")
                    sig_dict[windows[i]].append("NA")

    # Make and write outlines
    diff_outlines = []
    sig_outlines = []
    diff_outlines.append("Window\t" + "\t".join(out_header) + "\n")
    sig_outlines.append("Window\t" + "\t".join(out_header) + "\n")
    for w in windows:
        diff_outlines.append(w + "\t" + "\t".join([str(v) for v in diff_dict[w]]) + "\n")
        sig_outlines.append(w + "\t" + "\t".join([str(v) for v in sig_dict[w]]) + "\n")

    freq_output = open(outfile + ".FreqDiff","w")
    freq_output.write("".join(diff_outlines))
    freq_output.close()

    diff_output = open(outfile + ".SigDiff","w")
    diff_output.write("".join(sig_outlines))
    diff_output.close()

    return outfile + ".FreqDiff", outfile + ".SigDiff"

# Step 4
def SummarizeDifference(AlleleDiffMatrix,ExprConditions):
    ''' Average the differece values across expression conditions

    '''

    # Read the difference matrix
    inlines = [ln.strip() for ln in open(AlleleDiffMatrix,"r").readlines()]
    conditions = inlines[0].split("\t")
    
    # Make a dictionary of indexes for each ExprCondition
    condiction_index_dict = {}
    for exp_cond in ExprConditions:
        indexes = [idx for idx in range(len(conditions)) if exp_cond == conditions[idx].split("_vs_")[0]]
        condiction_index_dict[exp_cond] = indexes

    # For each line, average condictions accordng to condiction_index_dict
    outheader = "Window\t" + "\t".join(ExprConditions) + "\n"
    outlines = []
    for ln in inlines[1:]:
        split_ln = ln.split("\t")
        current_line = split_ln[0]
        for exp_cond in ExprConditions:
            values = [float(split_ln[i]) for i in condiction_index_dict[exp_cond] if not split_ln[i] == "NA"]
            if len(values) > 0:
                average = sum(values)/float(len(values))
            else:
                average = "NA"
            current_line = current_line + "\t" + str(average)
        outlines.append(current_line+"\n")

    # Write output
    output = open(AlleleDiffMatrix+".summary","w")
    output.write(outheader)
    output.write("".join(outlines))
    output.close()
    
    # Return file name
    return AlleleDiffMatrix+".summary"

def SummarizeDiffPvalue(AlleleSigMatrix,ExprConditions):
    ''' Average the log p-values across expression conditions

    '''
    
    import math

    # Read the difference matrix
    inlines = [ln.strip() for ln in open(AlleleSigMatrix,"r").readlines()]
    conditions = inlines[0].split("\t")

    # Make a dictionary of indexes for each ExprCondition
    condiction_index_dict = {}
    for exp_cond in ExprConditions:
        indexes = [idx for idx in range(len(conditions)) if exp_cond == conditions[idx].split("_vs_")[0]]
        condiction_index_dict[exp_cond] = indexes

    # For each line, average condictions accordng to condiction_index_dict
    outheader = "Window\t" + "\t".join(ExprConditions) + "\n"
    outlines = []
    for ln in inlines[1:]:
        split_ln = ln.split("\t")
        current_line = split_ln[0]
        for exp_cond in ExprConditions:
            values = [-math.log(float(split_ln[i]),10) for i in condiction_index_dict[exp_cond] if not split_ln[i] == "NA"]
            if len(values) > 0:
                average = sum(values)/float(len(values))
            else:
                average = "NA"
            current_line = current_line + "\t" + str(average)
        outlines.append(current_line+"\n")

    # Write output
    output = open(AlleleSigMatrix+".summary","w")
    output.write(outheader)
    output.write("".join(outlines))
    output.close()

    # Return file name
    return AlleleSigMatrix+".summary"

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

# MAIN

print'''

Wrapper script for time point analysis pipeline.

Takes a control file [1, see "ControlFileTemplate.ctl"] which
specifies all paramters

This pipeline involves the following steps, with scripts in </path/to/scripts>:
    -(1) Filter SNPs
    -(2) Calculte Windowed AlleleFreq
    -(3) Calcualte Difference in AlleFreq and Significance 

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
required_parameters = ["Inc","Ctl","Trt","BckSNP","AltSNP","Wsize","Witer","minSNP","Scripts","MatrixFile"]

for param in required_parameters:
    if not param in parameter_dictionary.keys():
        raise MissingParameter("Control file is missing a line for parameter: " + param) 
    else:
        print "Parameter " + param + " is present and valued"

print "All Inputs present and accounted for [" + FormatTime(time.time() - start_timer) + "]"

# Make sure the path to scripts is properly formatted
if not parameter_dictionary["Scripts"].endswith("/"):
    parameter_dictionary["Scripts"] = parameter_dictionary["Scripts"] + "/"

## Run Pipeline ##

# Step 1: Filter SNPs
#
# Functions for this section:
# - RunSNPFilter(read_summary,SNP,Scripts)
# - RunReverseCounts(filtered_file,Scripts)
# - MergeFiles(fileA,fileB,SNP_A,SNP_B,Scripts)
# - SummarizeStep1(summary_file,BckSNP,AltSNP,Scripts)

print "Start Step1: Filtering SNPS [" + FormatTime(time.time() - start_timer) + "]"

IncMergedFiles = []
for file in parameter_dictionary['Inc'].split(","):
    try:
        IncMergedFiles.append(SummarizeStep1(file,parameter_dictionary['BckSNP'],parameter_dictionary['AltSNP'],parameter_dictionary['Scripts']))
    except:
        print "Error processing " + file + "\nCotinuing without this file."

CtlMergedFiles = []
for file in parameter_dictionary['Ctl'].split(","):
    try:
        CtlMergedFiles.append(SummarizeStep1(file,parameter_dictionary['BckSNP'],parameter_dictionary['AltSNP'],parameter_dictionary['Scripts']))
    except:
        print "Error processing " + file + "\nCotinuing without this file."

TrtMergedFiles = []
for file in parameter_dictionary['Trt'].split(","):
    try:
        TrtMergedFiles.append(SummarizeStep1(file,parameter_dictionary['BckSNP'],parameter_dictionary['AltSNP'],parameter_dictionary['Scripts']))
    except:
        print "Error processing " + file + "\nCotinuing without this file."

print "Step1 Finsihed [" + FormatTime(time.time() - start_timer) + "]"

# Test Code
#print IncMergedFiles
#print CtlMergedFiles
#print TrtMergedFiles

# Step 2: Make a combined AlleleFreq Matrix
#
# Functions for this section:
# -RunWindowedAlleleFreq(file,wsize,witer,minSNP,Scripts)
# -RunCombineAlleFreq(files,outfile,Scripts)
# -RunCoverageCheck(allele_file,Scripts)

print "Start Step2: Caculate Windowed Allele Freq [" + FormatTime(time.time() - start_timer) + "]"

# Caculate the windowed AlleleFreq
IncWindowedFiles = [RunWindowedAlleleFreq(f,parameter_dictionary['Wsize'],parameter_dictionary['Witer'],parameter_dictionary['minSNP'],parameter_dictionary['Scripts']) for f in IncMergedFiles]
CtlWindowedFiles = [RunWindowedAlleleFreq(f,parameter_dictionary['Wsize'],parameter_dictionary['Witer'],parameter_dictionary['minSNP'],parameter_dictionary['Scripts']) for f in CtlMergedFiles]
TrtWindowedFiles = [RunWindowedAlleleFreq(f,parameter_dictionary['Wsize'],parameter_dictionary['Witer'],parameter_dictionary['minSNP'],parameter_dictionary['Scripts']) for f in TrtMergedFiles]

# Purge missing files from file list
IncWindowedFiles = [i for i in IncWindowedFiles if i]
CtlWindowedFiles = [c for c in CtlWindowedFiles if c]
TrtWindowedFiles = [t for t in TrtWindowedFiles if t]

# Test Code
#print IncWindowedFiles
#print CtlWindowedFiles
#print TrtWindowedFiles

# Write files to file
files = []
files.extend(IncWindowedFiles)
files.extend(CtlWindowedFiles)
files.extend(TrtWindowedFiles)
output = open("files.tmp","w")
output.write("\n".join(files))
output.close()

# Combine AlleFreq files
[CombinedAlleleFile,CombinedSNPFile] = RunCombineAlleFreq("files.tmp",parameter_dictionary['MatrixFile'],parameter_dictionary['Scripts'])

# Check coverage of windows
RunCoverageCheck(CombinedAlleleFile,parameter_dictionary['Scripts'])

print "Step2 Finished [" + FormatTime(time.time() - start_timer) + "]"

# Step3: Determine significance of change in AlleleFreq
#
# Functions for this section:
# -InocDifference(inoc1,inoc2,AlleleMatrix)
# -ExprDifference(expr_names,inoc_names,AlleleMatrix,mean,stdv,outfile)

print "Step3 Calculate Significance of change in Allele Frequency [" + FormatTime(time.time() - start_timer) + "]"

# Calculate difference between Inoculum files
try:
   [mean,stdv] = InocDifference(IncWindowedFiles[0],IncWindowedFiles[1],CombinedAlleleFile)
except IndexError as e:
   print "Insufficient Inoculum files for contrast"
   raise e

print "Inoculum Difference: " + str(mean) + "( +/- " + str(stdv) + " )"

# Get all of the Inoculum conditions
ExprConditions = []
ExprConditions.extend(CtlWindowedFiles)
ExprConditions.extend(TrtWindowedFiles)

# Calculate the significane of the difference between Inoclum and Experimental conditions
[AlleleDiffMatrix,AlleleSigMatrix] = ExprDifference(ExprConditions,IncWindowedFiles,CombinedAlleleFile,mean,stdv,parameter_dictionary['MatrixFile'])

print "Allele Difference File: " + AlleleDiffMatrix
print "Allele Difference P-value File: " + AlleleSigMatrix

print "Step3 Finished [" + FormatTime(time.time() - start_timer) + "]"

# Step 4: Summarize Results
#
# Functions for this section
# -SummarizeDifference(AlleleDiffMatrix,ExprConditions)
# -SummarizeDiffPvalue(AlleleSigMatrix,ExprConditions)

print "Step4 Summarize Results"

SummaryDiff = SummarizeDifference(AlleleDiffMatrix,ExprConditions)
SummarySigDiff = SummarizeDiffPvalue(AlleleSigMatrix,ExprConditions)

print "Summarized Allele Difference File: " + SummaryDiff
print "Summarized Allele SigDifference File: " + SummarySigDiff
