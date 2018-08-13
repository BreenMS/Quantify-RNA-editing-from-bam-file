library(GeneOverlap)
library(clusterRepro)
library(WGCNA)

#set minimum overlap of genes within a pathway for preservation analysis
minimumoverlap = 10; #there must be 10 genes in a pathway to test for preservation

case <- read.delim("Disease_samples_matrix.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(case) =names(case)
dataExpr0<-as.data.frame(t(case))
head(dataExpr0)
dim(dataExpr0)
gsg=goodSamplesGenes(dataExpr0,verbose=3)
gsg$allOK

control <- read.delim("Controls_samples_matrix.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(control) =names(control)
dataExpr1<-as.data.frame(t(control))
head(dataExpr1)
dim(dataExpr1)
gsg=goodSamplesGenes(dataExpr1,verbose=3)
gsg$allOK

Genes=colnames(dataExpr0)
nGenes=length(Genes)

list1 <- scan("Pathways4Preservation.gmt", what="", sep="\n")
list1a <- strsplit(list1, "[[:space:]]+")
names(list1a) <- sapply(list1a, `[[`, 1)
list1a <- lapply(list1a, `[`, -1)

for(i in 1:length(list1a))
{
PathwayName=names(list1a[i])
printFlush("=====================================");
printFlush(paste("Now working on pathway =>", PathwayName, "\n"));

GeneSet <- lapply(list1a[i], `[`, )
NumberGenes=length(GeneSet[[1]])
Convergent=  GeneSet[[1]][GeneSet[[1]] %in% Genes]
Divergent= NumberGenes-length(Convergent);
    
printFlush("MATCHING GENES:",length(Convergent), "      FALSE GENES:", Divergent);
printFlush("=====================================\n");
#doModulePreservation = TRUE;
#doClusterRepro = TRUE;
common=length(Convergent)

#Set
if (common > minimumoverlap) {
presGenes = c(1:nGenes)[GeneSet[[1]] %in% Genes];
preserved = "PathwayGene";
pathLabels0 = rep(preserved, length(Convergent));

pathwayLabels = rep(0, nGenes);
pathwayLabels[sample(presGenes, length(Convergent))] = pathLabels0;

multiColor = list();
multiColor = list(Control = pathwayLabels)
#Set multiExpr so that controls are the reference network
multiExpr  = list(Control=list(data=(dataExpr1)),PTSD=list(data=(dataExpr0)))

mp1=modulePreservation(multiExpr,multiColor,
referenceNetworks=1,verbose=3,networkType="unsigned", corFnc="bicor",
nPermutations=200,maxGoldModuleSize=common,maxModuleSize=common)

Zsummary1 = mp1$preservation$Z$ref.Control$inColumnsAlsoPresentIn.PTSD[,"Zsummary.pres"]
Zdensity1 = mp1$preservation$Z$ref.Control$inColumnsAlsoPresentIn.PTSD[,"Zdensity.pres"]
Zconnectivity1 = mp1$preservation$Z$ref.Control$inColumnsAlsoPresentIn.PTSD[,"Zconnectivity.pres"]
log.p1 = mp1$preservation$log.p$ref.Control$inColumnsAlsoPresentIn.PTSD[,"log.psummary.pres"]
p.value1 = exp(log.p1)
moduleSize1=mp1$preservation$Z$ref.Control$inColumnsAlsoPresentIn.PTSD$moduleSize
Zsummary.ControlRef.vPTSD = cbind(Zsummary1, Zdensity1, Zconnectivity1,log.p1, p.value1,moduleSize1)
rownames(Zsummary.ControlRef.vPTSD) = rownames(mp1$preservation$Z$ref.Control$inColumnsAlsoPresentIn.PTSD)
Digit=which(rownames(Zsummary.ControlRef.vPTSD)=="PathwayGene")
rownames(Zsummary.ControlRef.vPTSD)[Digit] <- PathwayName
write.table(Zsummary.ControlRef.vPTSD, file=paste("Zsummary.ControlRef.vPTSD.",PathwayName,".txt", sep=""), sep="\t")
printFlush(paste("\n Finished (Controls=reference network) for pathway", PathwayName, "\n"));

}


}

