setwd("/sc/orga/projects/buxbaj01a/Breen/NIMH/")

file_list=read.table("BamFiles.txt", stringsAsFactors=FALSE) #a text file listing all mapped bam files
files=basename(file_list$V1)
full_run_input <- NULL
full_run_input <- paste('#!/bin/bash\n', sep="")
for(item in 1:length(files)){
        file_name1=(files[item])
        sh_file <- paste(file_name1,".sh", sep="")
        full_run_input <- paste(full_run_input,'bsub < ', sh_file, ';\n sleep 1;\n', sep="")
        sh_runfile <- "run_all_samples.sh"
        fileConn <- file(sh_file)

writeLines(paste('
#!/bin/bash
#BSUB -W 08:00
#BSUB -n 5
#BSUB -q premium
#BSUB -P acc_buxbaj01a
#BSUB -cwd /sc/orga/projects/buxbaj01a/Breen/NIMH/
#BSUB -o /sc/orga/projects/buxbaj01a/Breen/NIMH/', file_name1, '.edit.log.out
#BSUB -e /sc/orga/projects/buxbaj01a/Breen/NIMH/', file_name1, '.edit.log.err
#BSUB -u michael.breen@mssm.edu
#BSUB -R rusage[mem=14000]
#BSUB -R span[hosts=1]
#BSUB -L /bin/bash

module load samtools
cd /sc/orga/projects/buxbaj01a/Breen/NIMH 

perl Step1_Query_Editing_Level.pl Radar_extended.txt ', file_name1, ' ', file_name1,'.editing.txt                                         #Identify RNA editing sites
perl Step2_OverallEditing.pl ',file_name1,'.editing.txt > ', file_name1,'.overallEditing.txt                                              #Quantify overall RNA editing per sample
perl Step3_OverallEditingbyRegion.pl ',file_name1,'.editing.txt >> ', file_name1,'.overallEditing.txt                              #Quantify overall RNA editing within defined genic regions per sample
perl Step4_OverallEditingbyGeneSet.pl Curated_GeneSets.txt ',file_name1,'.editing.txt >> ', file_name1,'.overallEditing.txt       #Quantify overall RNA editing within specific genesets per sample

echo Finished running ', file_name1, '', sep=""), fileConn)
        close(fileConn)

}
fileConn <- file(sh_runfile)
writeLines(paste(full_run_input, sep=""), fileConn)
close(fileConn)
