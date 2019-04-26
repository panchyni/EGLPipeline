# IMPORT
import sys
import time

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


def UpdateFileDictionary(file_list,file_dict):
    """ Takes file information form a file_list and assigns 
        each file to a group(s) in the file_dict
    """

    # Process files, and assign files to group in the dictionary
    for split_ln in file_list:
        file = split_ln[0]
        contig = split_ln[1]
        time = split_ln[2]
        ctl_idxs = [int(i)-1 for i in split_ln[3].split(",")]
        trt_idxs = [int(i)-1 for i in split_ln[4].split(",")]

        # If we already have a record for the current contig
        if contig in file_dict.keys():
            file_dict[contig]["all"].append((file,ctl_idxs,trt_idxs))
            # If we already have a record for the current time 
            if time in file_dict[contig].keys():
                file_dict[contig][time].append((file,ctl_idxs,trt_idxs))
            else:
                file_dict[contig][time] = [(file,ctl_idxs,trt_idxs)]
        else:
            file_dict[contig] = {}
            file_dict[contig]["all"] = [(file,ctl_idxs,trt_idxs)]
            file_dict[contig][time] = [(file,ctl_idxs,trt_idxs)]


def CheckValues(values,cutoff,sign):
    """ Compares a list of values to a cutoff with
        the sign of the comparison determined by sign
    """

    # Determine type of comparison
    comp_list = []
    if sign == "Plus":
        comp_list = [1 for v in values if v >= cutoff]
    elif sign == "Minus":
        comp_list = [1 for v in values if v <= cutoff]
    else:
        raise UnknownValue("Unkown values for comparsion sign",sign)

    # Evalute comparison list
    if sum(comp_list) == len(values) and len(values) > 0:
        return True
    else:
        return False

def Avg(values):
    """ Calculates average of a list of values """
 
    if len(values) > 0:
        return sum(values)/len(values)
    else:
        return "NA" 

def ContrastContigs(contig_files,cutoff,sign):
    """ Contrasts windows using data in contig_files
        the specified cutoff 

        Returns a list of regions which pass for all
        values, control values, treatment values,
        as well as regions unique to control and
        treatment
    """

    # Make dictionary for all, control, and treatment data
    all_dict = {}
    ctl_dict = {}
    trt_dict = {}

    sorted_windows = []

    # For each file
    for f in contig_files:

        #print f

        # Get file name and location indexes for ctl and treatment columns
        file = f[0]
        ctl_idxs = f[1]
        trt_idxs = f[2]
        
        file_lines = [ln.strip().split("\t") for ln in open(file,"r").readlines()[1:]] # Skip file header
        
        for split_ln in file_lines:
            #print split_ln
            window = split_ln[0]
            all_values = [float(v) for v in split_ln[1:] if not v == "NA"]
            ctl_values = [float(split_ln[i]) for i in ctl_idxs if not split_ln[i] == "NA"]
            trt_values = [float(split_ln[i]) for i in trt_idxs if not split_ln[i] == "NA"]
            
            if window in all_dict.keys():
                all_dict[window].extend(all_values)
                ctl_dict[window].extend(ctl_values)
                trt_dict[window].extend(trt_values)
            else:
                all_dict[window] = all_values
                ctl_dict[window] = ctl_values
                trt_dict[window] = trt_values
                sorted_windows.append(window)

    # For each window
    window_lines = []
    out_dict = {"all_out": [], "clt_out": [], "trt_out": [], "ctl_unique_out": [], "trt_unique_out": []}

    for window in sorted_windows:
        all_pass = CheckValues(all_dict[window],cutoff,sign)
        ctl_pass = CheckValues(ctl_dict[window],cutoff,sign)
        trt_pass = CheckValues(trt_dict[window],cutoff,sign)
 
        #window_values = ['0','0','0','0','0','0']
        window_values = ['0','0','0','0','0']

        if all_pass == True:
            out_dict["all_out"].append(window + "\t" + str(Avg(all_dict[window])) + "\n")
            window_values[0] = "1"
        if ctl_pass == True:
            out_dict["clt_out"].append(window + "\t" + str(Avg(ctl_dict[window])) + "\n")
            window_values[1] = "1"
        if trt_pass == True:
            out_dict["trt_out"].append(window + "\t" + str(Avg(trt_dict[window])) + "\n")
            window_values[2] = "1"
        if ctl_pass == True and trt_pass == False:
	    out_dict["ctl_unique_out"].append(window + "\t" + str(Avg(ctl_dict[window])) + "\n")
            window_values[3] = "1"
        if trt_pass == True and ctl_pass == False:
            out_dict["trt_unique_out"].append(window + "\t" + str(Avg(trt_dict[window])) + "\n")
            window_values[4] = "1"
        #if trt_pass == True and ctl_pass == True:
        #    out_dict["overlapping_out"].append(window + "\t" + str(Avg(all_dict[window])) + "\n")
        #    window_values[5] = "1"

        window_lines.append(window + "\t" + "\t".join(window_values) + "\n")

    return out_dict, window_lines

def CondenseRegions(file):
    """ For a list of files with annotated region data
        and condense regions which overlap into a single
        region

        Return a list of condensed regions

    """
 
    # Get regions from the contrast data file
    file_regions = [ln.strip().split("\t") for ln in open(file,"r").readlines()]

    # Iniatie files for condensing regions
    condense_regions = []
    first_region = file_regions[0]
    current_region = {"start": int(first_region[0].split("_")[0]), "stop": int(first_region[0].split("_")[1]), "vals": [float(first_region[1])]}

    # Condense region and update contig 
    for region in file_regions[1:]:

        # Read values from region
        [start,stop] = [int(v) for v in region[0].split("_")]
        value = float(region[1])

        # Combine overlapping regions
        if start > current_region["start"] and start <= current_region["stop"]:
            current_region["stop"] = stop
            current_region["vals"].append(value)
        else:
            condense_regions.append(str(current_region["start"]) + "_" + str(current_region["stop"]) + "\t" + str(sum(current_region["vals"])/len(current_region["vals"])) + "\n")
            current_region = {"start": start, "stop": stop, "vals": [value]}

    condense_regions.append(str(current_region["start"]) + "_" + str(current_region["stop"]) + "\t" + str(sum(current_region["vals"])/len(current_region["vals"])) + "\n")

    return condense_regions

def Filterregions(regions,size):
    """ For a list of regions from "CondenseRegions"
        filter them based on region size

    """
    
    select_regions = []    
    for r in regions:
        [start,stop] = r.split("\t")[0].split("_")
        if int(stop) - int(start) >= size:
            select_regions.append(r)
    
    return select_regions

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

Wrapper for post-processing and analysis 
of windowed regions of the genome

Takes a control file [1, see "ProcessingFileTemplate.ctl"] which
specifies all paramters

This pipeline involves the following steps:
 -(1) Building contrast regions
 -(2) Comparing contrast regions to specific threshold
 -(3) Condensing overlapping contrast regions into single regions

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
required_parameters = ["FreqDiffFiles","SigDiffFiles","FreqDiffCutoff","SigDiffCutoff","FreqDiffSign","SigDiffSign","RegionSize"]

for param in required_parameters:
    if not param in parameter_dictionary.keys():
        raise MissingParameter("Control file is missing a line for parameter: " + param)
    else:
        print "Parameter " + param + " is present and valued"

print "All Inputs present and accounted for [" + FormatTime(time.time() - start_timer) + "]"


## Run Pipeline ##

# Step 1: Build Contrasts #
#
# Functions for this section
# -UpdateFileDictionary(file_list,file_dict)

print "Start Step1: Build Contrast Groups [" + FormatTime(time.time() - start_timer) + "]"

# Frequency Difference

freq_diff_dict = {}
freq_diff_files = [ln.strip().split("\t") for ln in open(parameter_dictionary["FreqDiffFiles"],"r").readlines() if not ln.startswith("#")]
UpdateFileDictionary(freq_diff_files,freq_diff_dict)

# Significance of Difference

sig_diff_dict = {}
sig_diff_files = [ln.strip().split("\t") for ln in open(parameter_dictionary["SigDiffFiles"],"r").readlines() if not ln.startswith("#")]
UpdateFileDictionary(sig_diff_files,sig_diff_dict)

print "Step1 Done [" + FormatTime(time.time() - start_timer) + "]"

# Step 2: Compare Contrasts #
#
# Functions for this section
# -ContrastContigs(contig_files,cutoff,sign)

print "Start Step2: Do Contrast and Write [" + FormatTime(time.time() - start_timer) + "]"

# Frequency Difference
freq_cutoff = float(parameter_dictionary["FreqDiffCutoff"])
freq_sign = parameter_dictionary["FreqDiffSign"]
freq_files = []
for contig in freq_diff_dict.keys():
    
    # For each group in the contig
    for file_key in freq_diff_dict[contig].keys():

        out_dict, window_lines  = ContrastContigs(freq_diff_dict[contig][file_key],freq_cutoff,freq_sign)    
    
        base_out = "Contrast_FreqDiff_" + contig + "_" + file_key + "_"

        # Write region lists
        for out_key in out_dict.keys():

            output = open(base_out + out_key + ".passed_" + str(freq_cutoff) + "_" + freq_sign, "w")
            freq_files.append(base_out + out_key + ".passed_" + str(freq_cutoff) + "_" + freq_sign)
            output.write("".join(out_dict[out_key]))
            output.close()

        # Write window
        output = open(base_out + ".passed_" + str(freq_cutoff) + "_" + freq_sign + ".windows_matrix","w")
        output.write("".join(window_lines))
        output.close()
       

# Significant Difference
sig_cutoff = float(parameter_dictionary["SigDiffCutoff"])
sig_sign = parameter_dictionary["SigDiffSign"]
sig_files = []
for contig in sig_diff_dict.keys():

    # For each group in the contig
    for file_key in sig_diff_dict[contig].keys():
        out_dict, window_lines  = ContrastContigs(sig_diff_dict[contig][file_key],sig_cutoff,sig_sign)

        # Write region lists
        base_out = "Contrast_SigDiff_" + contig + "_" + file_key + "_"
        for out_key in out_dict.keys():
            output = open(base_out + out_key + ".passed_" + str(sig_cutoff) + "_" + sig_sign, "w")
            sig_files.append(base_out + out_key + ".passed_" + str(sig_cutoff) + "_" + sig_sign)
            output.write("".join(out_dict[out_key]))
            output.close()

        # Write window 
        output = open(base_out + ".passed_" + str(sig_cutoff) + "_" + sig_sign + ".windows_matrix", "w")
        output.write("".join(window_lines))
        output.close()

print "Step2 Done [" + FormatTime(time.time() - start_timer) + "]"

# Step 3: Condense Contrast Regions and Annotate Contig #
#
# Functions for this section
# -CondenseRegions(file)
# -Filterregions(regions,size)

print "Start Step3: Condense Regions [" + FormatTime(time.time() - start_timer) + "]"

# Frequency Difference

for file in freq_files:

    # Check that the file is not empty
    if open(file,"r").readline():
        condensed_regions = CondenseRegions(file)
        select_regions = Filterregions(condensed_regions,int(parameter_dictionary["RegionSize"]))

        # Write output
        output = open(file + ".condensed","w")
        output.write("".join(condensed_regions))
        output.close()

        output = open(file + ".condensed.select_" + parameter_dictionary["RegionSize"],"w")
        output.write("".join(select_regions))
        output.close()
  

# Significant Difference
for file in sig_files:
 
    # Check that the file is not empty
    if open(file,"r").readline():
        condensed_regions = CondenseRegions(file)
        select_regions = Filterregions(condensed_regions,int(parameter_dictionary["RegionSize"]))

        # Write output
        output = open(file + ".condensed","w")
        output.write("".join(condensed_regions))
        output.close()

        output = open(file + ".condensed.select_" + parameter_dictionary["RegionSize"],"w")
        output.write("".join(select_regions))
        output.close()

print "Step3 Done [" + FormatTime(time.time() - start_timer) + "]"
