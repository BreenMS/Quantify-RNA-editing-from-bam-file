# Identify and quantify RNA editing from mapped bam files

Written by: Michael S. Breen, PhD <br />
Contact: michael.breen@mssm.edu <br />

File descriptions: <br />
RapidPreservation.R = Rscript that will execute the preservation analysis. This script has several dependencies: 1) R libraries WGCNA, GeneOverlap and clusterRepro 2) two user defined matrices (e.g. case matrix and control matrix), which as a rule of thumb should contain greater than 15 samples per matrix, and 3) a curated .gmt file of pathways to test for preservation.   <br />
Pathways4Preservation.gmt = .gmt file of pathways to be tested <br />
