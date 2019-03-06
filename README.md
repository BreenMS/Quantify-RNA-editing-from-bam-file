# Identify and quantify RNA editing from mapped bam files

Written by: Michael S. Breen, PhD <br />
Contact: michael.breen@mssm.edu <br />

PART 1. Identify and quantify RNA editing in each mapped bam file <br />
Our code for quantifying RNA editing events uses a series of inter-linked perl scripts. There are four scripts that initate the pipeline. The most effecient way to run these scripts is to set up a shell script so as to run each script sequentially (e.g. see Rscript_1_callEditing.R). <br /> <br />

Step 1. [Step1_Query_Editing_Level.pl].<br /> 
This script will call RNA editing sites and compute editing rates across the genome for each site listed in a large RNA editing database. The user must specifiy three input files: <br /> 
1) A list of editing sites (e.g. RADAR database) <br /> 
2) INDEXED BAM alignment file <br /> 
3) Output file name <br /> <br /> 

Usage: perl Step1_Query_Editing_Level.pl Radar_Database.txt SampleName.bam SampleName.bam.editing.txt <br /> <br /> 

Step 2. [Step2_OverallEditing.pl].<br /> 
This script will compute a) the total number of detected sites, b) the total number of edited reads, c) the total coverage (i.e. reads overall), d) the overall editing rates (i.e. total numnber edited reads / overall editing rates). The user must specifiy one input file: <br /> 
1) The input file should be the output file from Step 1 containing all sites detected and editing levels <br /> <br /> 

Usage: perl Step2_OverallEditing.pl SampleName.bam.editing.txt > SampleName.bam.Overallrates.txt <br /> <br /> 

Step 3. [Step3_OverallEditingbyRegion.pl].<br /> 
This script will compute a) the total number of detected sites, b) the total number of edited reads, c) the total coverage (i.e. reads overall), d) the overall editing rates based on genic regions (e.g. 3UTRs, 5UTRs, exonic etc..) The user must specifiy one input file: <br /> 
1) The input file should be the output file from Step 1 containing all sites detected and editing levels <br /> <br /> 

Usage: perl Step3_OverallEditingbyRegion.pl SampleName.bam.editing.txt > SampleName.bam.OverallratesbyRegion.txt <br /> <br /> 

Step 4. [Step4_OverallEditingbyGeneSet.pl].<br /> 
This script will compute a) the total number of detected sites, b) the total number of edited reads, c) the total coverage (i.e. reads overall), d) the overall editing rates based on curated genesets (e.g. glutamatergic and serotonergic receptors etc..) The user must specifiy one input file: <br /> 
1) The input file should be the output file from Step 1 containing all sites detected and editing levels <br /> <br /> 

Usage: perl Step4_OverallEditingbyGeneSet.pl Curated_GeneSets.txt SampleName.bam.editing.txt > SampleName.bam.OverallratesbyGeneSets.txt <br /> <br /> 


