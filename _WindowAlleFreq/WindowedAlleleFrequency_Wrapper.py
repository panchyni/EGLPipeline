# NOTES:
# Considering adding a optional file path to append to being of summary and SNP files
#

# IMPORT
import sys
import os

# FUNCTIONS
def RunSNPFilter(read_summary,SNP_file,SNP_id,Scripts):
    ''' Runs "0.1a_v2_FilterVCFSummaryFile.py" with the specified files and parameters

        Returns the name of the outputfile
    '''

    command = "python " + Scripts + "0.1a_v2_FilterVCFSummaryFile.py -summary " + read_summary + " -SNPlist " + SNP_file + " -SNPid " + SNP_id + "\n"
    os.system(command)
    return read_summary + ".filtered_" + SNP_id


def RunSNPFilter_Imputation(read_summary,SNP_file,SNP_id,Scripts):
    ''' Runs "0.1a_v3_FilterVCFSummaryFile_wImputation.py" with the specificed files and parameters

        Return the name of the outputfile
    '''

    command = "python " + Scripts + "0.1a_v3_FilterVCFSummaryFile_wImputation.py -summary " + read_summary + " -SNPlist " + SNP_file + " -SNPid " + SNP_id + "\n"
    os.system(command)
    return read_summary + ".filtered_imputed_" + SNP_id

def RunReverseCounts(filtered_file,Scripts):
    ''' Runs "0.1b_ReverseSNPCounts.py" on the speicified files

        Return the name of the outputfile
    '''

    command = "python " + Scripts + "0.1b_ReverseSNPCounts.py -filtered_SNPs " + filtered_file + "\n"
    os.system(command)
    return filtered_file + ".rev"

def RunReverseCounts(filtered_file,Scripts):
    ''' Runs "0.1b_ReverseSNPCounts.py" on the speicified files

        Return the name of the outputfile
    '''

    command = "python " + Scripts + "0.1b_ReverseSNPCounts.py -filtered_SNPs " + filtered_file + "\n"
    os.system(command)
    return filtered_file + ".rev"

def MergeFiles(fileA,fileB,SNP_idA,SNP_idB,Scripts):
    ''' Runs the "0.1c_MergeFilteredSNPs.py" with the specified files and paramters

    '''

    if fileA.endswith(".rev"):
        root = ".".join(fileA.split(".")[0:-2])
    else:
        root = ".".join(fileA.split(".")[0:-1])

    command = "python " + Scripts + "0.1c_MergeFilteredSNPs.py -fileA " + fileA + " -fileB " + fileB + " -labelA " + SNP_idA + " -labelB " + SNP_idB + "\n"
    os.system(command)
    return root + ".Combined"

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
def main():
    print'''

    Wrapper script for generating windowed allele frequency files

    Takes a control file (See "WindowAllele_CtlFile.txt")


    '''

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
    required_parameters = ["SumFiles","BckSNP","AltSNP","BckLabel","AltLabel","Wsize","Witer","minSNP","Scripts"]

    for param in required_parameters:
        if not param in parameter_dictionary.keys():
            raise MissingParameter("Control file is missing a line for parameter: " + param)

    # Check optional parameters
    imputation = "FALSE"
    if 'Imputation' in parameter_dictionary.keys():
        imputation = parameter_dictionary['Imputation']

    # Order SNP by Chromosome
    BckSNP_dict = {}
    for bck_file in parameter_dictionary['BckSNP'].split(","):
        bck_file_chr = bck_file.split(".")[-1]
        BckSNP_dict[bck_file_chr] = bck_file

    AltSNP_dict = {}
    for alt_file in parameter_dictionary['AltSNP'].split(","):
        alt_file_chr = alt_file.split(".")[-1]
        AltSNP_dict[alt_file_chr] = alt_file

    # Simplify some variables
    BckLabel = parameter_dictionary['BckLabel']
    AltLabel = parameter_dictionary['AltLabel']
    Scripts =  parameter_dictionary['Scripts']


    # Go through each summary file
    processed_files = []
    for file in parameter_dictionary['SumFiles'].split(","):

        # Get all contigs
        contigs = list(set([l.strip().split("\t")[0] for l in open(file,"r").readlines()]))
        chromosomes = [c for c in contigs if c.startswith("chromo")]

        # For each contig select line and write an outputfile
        for c in chromosomes:
            
            # Seperate out chromosomes 
            chr_file = file + "." + c
            command = 'grep -P "' + c + '\\t" ' + file + ' > ' + chr_file
            os.system(command)

            # Filter SNPs
            if imputation == 'TRUE':
                BckFilter = RunSNPFilter_Imputation(chr_file,BckSNP_dict[c],BckLabel,Scripts)
                AltFilter = RunSNPFilter_Imputation(chr_file,AltSNP_dict[c],AltLabel,Scripts)
            else:
                BckFilter = RunSNPFilter(chr_file,BckSNP_dict[c],BckLabel,Scripts)
                AltFilter = RunSNPFilter(chr_file,AltSNP_dict[c],AltLabel,Scripts)

            # Reverse freq for AltSNP file
            RevAltFilter = RunReverseCounts(AltFilter,Scripts)

            # Merge Files
            MergedFile = MergeFiles(BckFilter,RevAltFilter,BckLabel,AltLabel,Scripts)
            os.system("mv " + MergedFile + " " + MergedFile + ".Rel_" + BckLabel)

            MergedFile = MergedFile + ".Rel_" + BckLabel

            # Calculate Windowed AlleleFreq
            AlleleFreqFile = RunWindowedAlleleFreq(MergedFile,parameter_dictionary["Wsize"],parameter_dictionary["Witer"],parameter_dictionary["minSNP"],Scripts)
            processed_files.append(AlleleFreqFile)

# If not loaded as a module
if __name__ == '__main__':
    main()

