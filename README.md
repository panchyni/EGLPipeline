# The EGL annotation pipeline is executed using the '0_ComparsionToInocWrapper_vx.py' (where x is the version number) script
# and a control file which will be passed to the wrapper script in the following format

python 0_ComparsionToInocWrapper_v5.py ControlFile.txt

# Everything else required for the pipeline is described in the control file, which has the following format:

Inc     {List of processed VCF files for the target chromosome from the inoculum or simulated inoculum, comma seperated}
Ctl     {List of processed VCF files for the target chromosome from the control experiments, comma seperated}
Trt     {List of processed VCF files for the target chromosome from the treatment experiments, comma seperated}
BckSNP  {Background SNP locations,SNP annotation Name; example: bwa_msdr_MR_ih_lc_nr503_f.vcf.1009vs2343.unique.snps.CR2343.PASS.filter_5.positions.chromosome_15,CR2343]
AltSNP  {Alternative SNP locations,SNP annotation Name; example: bwa_msdr_MR_ih_lc_nr503_f.vcf.1009vs2343.unique.snps.CR1009.PASS.filter_5.positions.chromosome_15,CR1009}
Wsize   {size of window in bases, interger value: example: 10000}
Witer   {size of window steps in bases, interger value: example: 2000}
minSNP  {minimum number of SNPs per window for analysis: example: 5}
Scripts {pather to pipelien scritps; example: /mnt/home/panchyni/4_GitLocal/EGLPipline}
MatrixFile      {base name for matrix outputfile; example: CombinedMatrix_runLR_time1213_chr15}

# The wrapper script uses the data described in the control file to execute the pipeline, which has the following steps (note that the wrapper script
#is broken up into the same sections

Step 1: Filter SNPs
 - Filter SNPs by location, reverse counts for background SNPs and combine

 Functions for this section:
 - RunSNPFilter(read_summary,SNP,Scripts)
 - RunReverseCounts(filtered_file,Scripts)
 - MergeFiles(fileA,fileB,SNP_A,SNP_B,Scripts)
 - SummarizeStep1(summary_file,BckSNP,AltSNP,Scripts)

Step 2: Make a combined AlleleFreq Matrix
 - Calculate SNP AlleleFreq using the described window size, step and MinSNPs and merge the results
 - Do this for each set of files (Inc, Ctrl, Trt)

 Functions for this section:
 -RunWindowedAlleleFreq(file,wsize,witer,minSNP,Scripts)
 -RunCombineAlleFreq(files,outfile,Scripts)
 -RunCoverageCheck(allele_file,Scripts)

Step3: Determine significance of change in AlleleFreq
 -Model distribution between incolumn samples (Inoc vs Inoc)
 -Calculate change in AlleleFreq (Ctrl vs Inoc and Trt vs Inoc)
 -Determine significance of difference by comparing it to the Inoc vs Inoc distribution

 Functions for this section:
 -InocDifference(inoc1,inoc2,AlleleMatrix)
 -ExprDifference(expr_names,inoc_names,AlleleMatrix,mean,stdv,outfile)

Step 4: Summarize Results
 -Summerize results

 Functions for this section
 -SummarizeDifference(AlleleDiffMatrix,ExprConditions)
 -SummarizeDiffPvalue(AlleleSigMatrix,ExprConditions)
