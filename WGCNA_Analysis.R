bsub -q alloc -P acc_buxbaj01a -n 48 -R rusage[mem=12000] -R span[ptile=4] -W 05:00 -m manda -Ip /bin/bash

library(WGCNA)
library(cluster)
options(stringsAsFactors  =  FALSE)
allowWGCNAThreads(n=34)

##################################################################################################################################################
######### HERE WE CREATE ONE LARGE NETWORK FOR EACH DISEASED CONDITION  ##################
###################################################################################################################################################
#INPUT PRE
PRE <- read.delim("Norm_median_clean.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(PRE) =names(PRE)
dataExpr0<-as.data.frame(t(PRE))
head(dataExpr0)
dim(dataExpr0)
gsg=goodSamplesGenes(dataExpr0,verbose=3)
gsg$allOK

nGenes = ncol(dataExpr0)
nSamples = nrow(dataExpr0)

traitData <- read.delim("some.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData)

#MATCH TRAITS TO PD-DEPLOYMENT
rowsExpr <- rownames(dataExpr0)
traitRows <- match(rowsExpr,traitData$SID.1)
datTraits = traitData[traitRows, -1];
rownames (datTraits) = traitData[traitRows, 1];
table(rownames(datTraits)==rownames(dataExpr0)) # EVERYTHING OK?
names(datTraits)

powers=c(1:15) # in practice this should include powers up to 20.
sft0=pickSoftThreshold(dataExpr0,powerVector=powers, networkType="signed")

pdf("Soft_Threshold_q.pdf")
par(mfrow=c(1,2))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],labels=powers,col="red")
abline(h=0.80,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

#M1: 8
#M2: 10
#M3: test 5 then 7
#M4: 6
adjacencyPre = adjacency((dataExpr0),power=7 ,type="signed") #FOR ONE LARGE NETWORK WE CHOOSE A POWER OF 10
diag(adjacencyPre)=0
dissTOMPre   = 1-TOMsimilarity(adjacencyPre, TOMType="signed")
geneTreePre  = hclust(as.dist(dissTOMPre), method="average")

pdf("GeneTree.pdf",height=6,width=12)
par(mfrow=c(1,1))
plot(geneTreePre,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity", labels=FALSE,hang=0.04);
dev.off()

#MODULE ASSIGNMENTS
mColorh=NULL
for (ds in 0:4){
 tree = cutreeHybrid(dendro = geneTreePre, pamStage=FALSE,
   minClusterSize = (50), cutHeight = 0.99, 
   deepSplit = ds, distM = dissTOMPre)
 mColorh=cbind(mColorh,labels2colors(tree$labels));
}

pdf("DeepSplit_Choices_q.pdf", height=10,width=25); 
plotDendroAndColors(geneTreePre, mColorh, paste("dpSplt =",0:4), main = "Co-Expression Network",dendroLabels=FALSE);
dev.off()

#SET DEEP SPLIT CHOICE AND NAME OUR COLORS
modulesPRE =  mColorh[,3]
table(modulesPRE)

#Check to see if disease network modules can be cut and merged...
MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method ="average")

pdf("Module_Relationshipsq.pdf")
plot(METree, main ="Clustering of Module Eigengenes",xlab ="",sub="")
abline(h=0.1,col="red")
abline(h=0.15,col="red")
abline(h=0.2,col="red")
abline(h=0.25,col="red")
dev.off()
MEDissThres = 0.25

####################################### SHOULD WE MERGE MODULES BASED ON THE ABOVE? #######################################
merge = mergeCloseModules(dataExpr0, modulesPRE, cutHeight = MEDissThres, verbose=3)
mergedColors = merge$colors
table(mergedColors)
table(modulesPRE)
mergedMEs = merge$newME

#COMPARE UNMERGED MODULES TO MERGED MODULES
pdf("PreDeploy_Network_Unmerged_Merged.pdf", w=9)
plotDendroAndColors(geneTreeControlPre, cbind(modulesPRE, mergedColors),c("Dynamic Tree Cut","Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

modulesPRE = merge$colors
table(modulesPRE)
MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method ="average")

##########################################################################################################################
##Our Options
# deep split = 4
# Min Module Size = 60
# Merge past = 0.2
# cutHeight = 0.99

#CALCULATE PC FOR VISUALIZATION FOR CASE PD-DEPLOYMENT
PCsPD    = moduleEigengenes((dataExpr0),  colors=modulesPRE) 
ME_PD    = PCsPD$eigengenes
distPCPD = 1-abs(cor(ME_PD,use="p"))
distPCPD = ifelse(is.na(distPCPD), 0, distPCPD)
pcTreePD = hclust(as.dist(distPCPD),method="average") 
MDS_PD   = cmdscale(as.dist(distPCPD),2)
colorsPD = names(table(modulesPRE))
names = row.names((dataExpr0))

pdf("Module_Visualizationq1.pdf",height=8,width=8)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 3) + 0.1, cex=1)
plot(pcTreePD, xlab="",ylab="",main="",sub="")
plot(MDS_PD, col= colorsPD,  main="MDS plot", cex=2, pch=19)

for (which.module in names(table(modulesPRE)))
{
 par(mfrow=c(2,1), mar=c(4, 4.1, 4.1, 2))
 plotMat(t(scale(dataExpr0[,modulesPRE==which.module])),
,cex.axis=2,nrgcols=100,rlabels=F,tck=0, rcols=which.module,main=paste("Heatmap",which.module,"Module"))

  ME = ME_PD[, paste("ME",which.module, sep="")] 
  n<- barplot(ME, col=which.module, cex.main=1, ylab="Eigengene Expression",xlab="")
  axis(1,at=n, labels=row.names(dataExpr0), las=2, cex.axis=0.4, font=2)
};
dev.off();


#DETERMINE GENE-SIG. ACROSS MODULES
names(datTraits)

Depression = as.data.frame(datTraits$Depression)
names(Depression)="Depression"
GS.Depression=as.numeric(cor(dataExpr0,Depression,use="p"))
GS.DepressionColor=numbers2colors(GS.Depression,signed=T)

PTSD.Dep = as.data.frame(datTraits$PTSD.Dep)
names(PTSD.Dep)="PTSD.Dep"
GS.PTSD.Dep=as.numeric(cor(dataExpr0,PTSD.Dep,use="p"))
GS.PTSD.DepColor=numbers2colors(GS.PTSD.Dep,signed=T)

PTSD = as.data.frame(datTraits$PTSD)
names(PTSD)="PTSD"
GS.PTSD=as.numeric(cor(dataExpr0,PTSD,use="p"))
GS.PTSDColor=numbers2colors(GS.PTSD,signed=T)

TE = as.data.frame(datTraits$TE)
names(TE)="TE"
GS.TE=as.numeric(cor(dataExpr0,TE,use="p"))
GS.TEColor=numbers2colors(GS.TE,signed=T)

Control = as.data.frame(datTraits$Control)
names(Control)="Control"
GS.Control=as.numeric(cor(dataExpr0,Control,use="p"))
GS.ControlColor=numbers2colors(GS.Control,signed=T)

Nicotine = as.data.frame(datTraits$Nicotine)
names(Nicotine)="Nicotine"
GS.Nicotine=as.numeric(cor(dataExpr0,Nicotine,use="p"))
GS.NicotineColor=numbers2colors(GS.Nicotine,signed=T)

Alcohol = as.data.frame(datTraits$Alcohol)
names(Alcohol)="Alcohol"
GS.Alcohol=as.numeric(cor(dataExpr0,Alcohol,use="p"))
GS.AlcoholColor=numbers2colors(GS.Alcohol,signed=T)

IPV = as.data.frame(datTraits$IPV)
names(IPV)="IPV"
GS.IPV=as.numeric(cor(dataExpr0,IPV,use="p"))
GS.IPVColor=numbers2colors(GS.IPV,signed=T)

Delivery = as.data.frame(datTraits$Delivery.natural.)
names(Delivery)="Delivery"
GS.Delivery=as.numeric(cor(dataExpr0,Delivery,use="p"))
GS.DeliveryColor=numbers2colors(GS.Delivery,signed=T)

Gender = as.data.frame(datTraits$Gender)
names(Gender)="Gender"
GS.Gender=as.numeric(cor(dataExpr0,Gender,use="p"))
GS.GenderColor=numbers2colors(GS.Gender,signed=T)

RIN = as.data.frame(datTraits$RIN)
names(RIN)="RIN"
GS.RIN=as.numeric(cor(dataExpr0,RIN,use="p"))
GS.RINColor=numbers2colors(GS.RIN,signed=T)

Batch = as.data.frame(datTraits$Batch)
names(Batch)="Batch"
GS.Batch=as.numeric(cor(dataExpr0,Batch,use="p"))
GS.BatchColor=numbers2colors(GS.Batch,signed=T)

Chip = as.data.frame(datTraits$Chip)
names(Chip)="Chip"
GS.Chip=as.numeric(cor(dataExpr0,Chip,use="p"))
GS.ChipColor=numbers2colors(GS.Chip,signed=T)



datColors0=data.frame(modulesPRE, GS.DepressionColor,GS.PTSD.DepColor,GS.PTSDColor, GS.TEColor,GS.ControlColor)


pdf("Drakenstein_ColoredModules_q2.pdf",height=8,width=14)
plotDendroAndColors(geneTreePre, colors=datColors0, main="Gene Dendrogram and Module Colors (Beta=8)", groupLabels=c("Module colors", "Depression", "PTSD.Dep","PTSD", "TE", "Control"), dendroLabels=FALSE, hang=0.03, addGuide=FALSE, guideHang=0.05) 
dev.off()


# PLOT MODULE-TRAIT RELATIONSHIP
MEs0  =  moduleEigengenes(dataExpr0,  modulesPRE)$eigengenes
MEs = orderMEs(MEs0)
write.table(MEs, "Drakenstein_MEs_q1.txt", sep="\t")

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
data = cbind(moduleTraitCor,moduleTraitPvalue)
pdf("Drakenstein_Correlationsq3.pdf")
textMatrix = paste(signif(moduleTraitCor, 3), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.4,
zlim = c(-.3,.3),
main = paste("ME-Trait Relationships"))
dev.off()


ySymbols = paste ("M", (unique(modulesPRE)), ": ", table (modulesPRE), sep="")

labeledHeatmap(Matrix=numMatlarge,
               xLabels=xLabels, #xSymbols=xSymbols,
               yLabels=yLabels, ySymbols=ySymbols,


pdf("Module_GS_q.pdf")
par(mfrow=c(3,3))
plotModuleSignificance(GS.Depression, modulesPRE, boxplot = FALSE, main = "", ylab = "MS (cor, Dep)", las=2, cex.axis=0.7)
plotModuleSignificance(GS.PTSD, modulesPRE, boxplot = FALSE, main = "", ylab = "MS (cor, PTSD)", las=2, cex.axis=0.7)
plotModuleSignificance(GS.TE, modulesPRE, boxplot = FALSE, main = "", ylab = "MS (cor, TE)", las=2, cex.axis=0.7)
plotModuleSignificance(GS.PTSD.Dep, modulesPRE, boxplot = FALSE, main = "", ylab = "MS (cor, PTSD.Dep)", las=2, cex.axis=0.7)
dev.off()

pdf("MDS_q1.pdf")
par(mfrow=c(2,2))
plot(MDS_PD, col= colorsPD,  main="MDS plot", cex=2, pch=16, las=1, cex.axis=0.8)
abline(h=0,lty=2)
abline(v=0,lty=2)
dev.off() 

#DETERMINE MODULE SIG. USING -LOG PVALUE
datSummary <- read.delim("Master.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
GS_Depression=-log10(datSummary$Depression)
GS_PTSD=-log10(datSummary$PTSD_c)
GS_TE=-log10(datSummary$TE)
GS_PTSD.Dep=-log10(datSummary$PTSD.Dep_c)

GS_Depression_Color=numbers2colors(GS_Depression,signed=T)
GS_PTSD_Color=numbers2colors(GS_PTSD,signed=T)
GS_TE_Color=numbers2colors(GS_TE,signed=T)
GS_PTSD.Dep_Color=numbers2colors(GS_PTSD.Dep,signed=T)

colfunc <- colorRampPalette(c("steelblue4","white"))

pdf("ModuleSignificance_q8.pdf")
par(mfrow=c(3,3))
verboseBarplot(GS_Depression,modulesPRE,color=colorsPD, border=colorsPD,main="Depression" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9, ylim=c(0,1.2))
abline(h=0.6,col="black",lty = 2)   
verboseBarplot(GS_PTSD,modulesPRE,color=colorsPD,border=colorsPD,main="PTSD" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = T,ylim=c(0,1.2))
abline(h=0.6,col="black",lty = 2)   
verboseBarplot(GS_TE,modulesPRE,color=colorsPD,border=colorsPD,main="TE" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = T,ylim=c(0,1.2))
abline(h=0.6,col="black",lty = 2)    
verboseBarplot(GS_PTSD.Dep,modulesPRE,color=colorsPD,border=colorsPD,main="PTSD.Dep" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = T,ylim=c(0,1.2))
abline(h=0.6,col="black",lty = 2)   
dev.off()


##CELL-Type Enrichment Per Module
Gene  = colnames(dataExpr0)
enrichments = userListEnrichment(Gene, modulesPRE,
fnIn = NULL,
catNmIn = NULL,
#fnIn = c("GeneList","ModuleColors"),
#catNmIn =  c("Genes","Modules"),
nameOut = "ModuleEnrichment2.csv", 
useBrainLists = TRUE,
useBloodAtlases = TRUE,
omitCategories = "grey", 
outputCorrectedPvalues = TRUE,
useStemCellLists = TRUE, 
outputGenes = TRUE, 
minGenesInCategory = 2, 
useBrainRegionMarkers = TRUE,
useImmunePathwayLists = TRUE,
usePalazzoloWang = TRUE)


#PLOT CONNECTIVITY VS GS
CM_Pre=intramodularConnectivity(adjacencyPre,colors=modulesPRE)
names(CM_Pre) 
whichmodule="red";
restrict1=modulesPRE==whichmodule
pdf("Red_Connectivity.pdf")
verboseScatterplot (CM_Pre$kWithin[restrict1],GS.Eventual_PTSD[restrict1],col=modulesPRE[restrict1],
xlab="Connectivity (k) ",ylab="G.S. Eventual PTSD(p)") 
dev.off()


# CALCULATE MODULE MEMBERSHIP VALUES (aka. module eigengene based connectivity kME)
datKME=signedKME(dataExpr0, MEs)
colorOfColumn=substring(names(datKME),4)

pdf("Regress_Modules_q9.pdf",w=13,h=5)
par(mfrow = c(2,5))
 for (module in names(table(modulesPRE))) 
  {
  column = match(module,colorOfColumn) 
  restModule=modulesPRE==module
  
  verboseScatterplot(datKME[restModule,column],GS.Depression[restModule],
  xlab=paste("MM ",module,""),ylab="GS.Depression",
  main=paste("kME vs GS"),col=module,abline=T,las=1, cex.axis=0.9, pch=16)
  
      verboseScatterplot(datKME[restModule,column],GS.PTSD.Dep[restModule],
  xlab=paste("MM ",module,""),ylab="GS.PTSD.Dep",
  main=paste("kME vs GS"),col=module,abline=T,las=1, cex.axis=0.9, pch=16)

    verboseScatterplot(datKME[restModule,column],GS.PTSD[restModule],
  xlab=paste("MM ",module,""),ylab="GS.PTSD",
  main=paste("kME vs GS"),col=module, abline=T,las=1, cex.axis=0.9, pch=16)

    verboseScatterplot(datKME[restModule,column],GS.TE[restModule],
  xlab=paste("MM ",module,""),ylab="GS.TE",
  main=paste("kME vs GS"),col=module,abline=T,las=1, cex.axis=0.9, pch=16)
    
        verboseScatterplot(datKME[restModule,column],GS.Control[restModule],
  xlab=paste("MM ",module,""),ylab="GS.Control",
  main=paste("kME vs GS"),col=module,abline=T,las=1, cex.axis=0.9, pch=16)
  }
dev.off()

#MODULE MEMBERSHIP (kME) KME is defined as correlation between expression and modules
#USED TO MEASURE CORRELATIONS BETWEEN EACH GENE AND EACH MODULE EIGENGENE
geneModuleMembership1 = signedKME((dataExpr0), MEs)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(dataExpr0)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");

Gene       = rownames(t(dataExpr0))
kMEtable1  = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))
write.csv(kMEtable1,"Drakenstein_kMEtable.csv",row.names=FALSE)

topGenesKME = NULL
for (c in 1:length(colorsPD)){
 kMErank1    = rank(-geneModuleMembership1[,c])
 maxKMErank  = rank(apply(cbind(kMErank1+.00001),1,max))
 topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=20])
}; colnames(topGenesKME) = colorsPD
topGenesKME

#WRITE OUTPUT
Pre<-as.data.frame(datTraits)
names(Pre)<-"Pre"
modNames = substring(names(MEs), 3)
nGenes = ncol(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
nSamples = nrow(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
geneModuleMembership<-as.data.frame(cor(dataExpr0,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr0, Pre, use = "p")) #  Generate correlations and p-value for each gene against the trait.  
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Pre), sep="")
names(GSPvalue) = paste("p.GS.", names(Pre), sep="")

Colors=modulesPRE
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="Network_WGCNA_OUTPUT_q.csv")

#HERE WE EXPORT 20 TOP INTRAMODULAR HUB GENES PER MODULE OF INTEREST

modules = c("blue") #"red", "yellow", "lightgreen") #,"yellow", "lightgreen", "black"
probes=names(dataExpr0)
inModule=is.finite(match(modulesPRE, modules));
modProbes=probes[inModule] ##
modTom=dissTOMPre[inModule, inModule]
dimnames(modTom)=list(modProbes, modProbes)
nTopHubs = 50
kIN = softConnectivity(dataExpr0[, modProbes])
selectHubs = (rank (-kIN) <= nTopHubs)

#EXPORT ENTIRE NETWORK
cyt=exportNetworkToCytoscape(modTom,
edgeFile=paste("CytoscapeInput_", paste(modules,collapose="_") , "edges.txt", sep=""),
nodeFile=paste("CytoscapeInput_", paste(modules,collapose="_") , "nodes.txt", sep=""),
weighted=TRUE,
threshold=0.50,
nodeNames=modProbes,
nodeAttr=modulesPRE[inModule])


geneModuleMembership1 = signedKME((dataExpr0), ME_PD)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(dataExpr1)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");
Gene= rownames(t(dataExpr1))
kMEtable1 = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
             kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))
write.csv(kMEtable1,"Control_Pre_kMEtable_for_ControlModules.csv",row.names=FALSE)

topGenesKME = NULL
for (c in 1:length(colorsPD)){
 kMErank1    = rank(geneModuleMembership[,c])
 maxKMErank  = rank(apply(cbind(kMErank1),1,max))
 topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; colnames(topGenesKME) = colorsPD
#PRINT TOP 10 GENES PER MODULE BASED ON KME IN BOTH NETWORKS
topGenesKME




matrix <- read.delim("list3.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
head(matrix)

par(mfrow=c(3,4))
boxplot(matrix$MEblue[1:31],matrix$MEblue[32:44],matrix$MEblue[45:64],matrix$MEblue[65:87], matrix$MEblue[88:149],medcol="blue", lwd=1.3, las=2, cex.axis=0.6, cex=.2, notch=F,xaxt='n')
abline(h=0,lty=2,col="black",lwd=0.8)
axis(side=1,at=c(1,2,3,4,5), cex.axis=0.6,labels=c("Depression","PTSD.Dep","PTSD","TE", "Control"),las=2)

boxplot(matrix$MEpink[1:31],matrix$MEpink[32:44],matrix$MEpink[45:64],matrix$MEpink[65:87], matrix$MEpink[88:149],medcol="pink", lwd=1.3, las=2, cex.axis=0.6, cex=.2, notch=F,xaxt='n')
abline(h=0,lty=2,col="black",lwd=0.8)
axis(side=1,at=c(1,2,3,4,5), cex.axis=0.6,labels=c("Depression","PTSD.Dep","PTSD","TE", "Control"),las=2)

boxplot(matrix$MEturquoise[1:31],matrix$MEturquoise[32:44],matrix$MEturquoise[45:64],matrix$MEturquoise[65:87], matrix$MEturquoise[88:149],medcol="turquoise", lwd=1.3, las=2, cex.axis=0.6, cex=.2, notch=F,xaxt='n')
abline(h=0,lty=2,col="black",lwd=0.8)
axis(side=1,at=c(1,2,3,4,5), cex.axis=0.6,labels=c("Depression","PTSD.Dep","PTSD","TE", "Control"),las=2)

boxplot(matrix$MEbrown[1:31],matrix$MEbrown[32:44],matrix$MEbrown[45:64],matrix$MEbrown[65:87], matrix$MEbrown[88:149],medcol="brown", lwd=1.3, las=2, cex.axis=0.6, cex=.2, notch=F,xaxt='n')
abline(h=0,lty=2,col="black",lwd=0.8)
axis(side=1,at=c(1,2,3,4,5), cex.axis=0.6,labels=c("Depression","PTSD.Dep","PTSD","TE", "Control"),las=2)

boxplot(matrix$MEcyan[1:31],matrix$MEcyan[32:44],matrix$MEcyan[45:64],matrix$MEcyan[65:87], matrix$MEcyan[88:149],medcol="cyan", lwd=1.3, las=2, cex.axis=0.6, cex=.2, notch=F,xaxt='n')
abline(h=0,lty=2,col="black",lwd=0.8)
axis(side=1,at=c(1,2,3,4,5), cex.axis=0.6,labels=c("Depression","PTSD.Dep","PTSD","TE", "Control"),las=2)

boxplot(matrix$MEmagenta[1:31],matrix$MEmagenta[32:44],matrix$MEmagenta[45:64],matrix$MEmagenta[65:87], matrix$MEmagenta[88:149],medcol="magenta", lwd=1.3, las=2, cex.axis=0.6, cex=.2, notch=F,xaxt='n')
abline(h=0,lty=2,col="black",lwd=0.8)
axis(side=1,at=c(1,2,3,4,5), cex.axis=0.6,labels=c("Depression","PTSD.Dep","PTSD","TE", "Control"),las=2)

boxplot(matrix$MEyellow[1:31],matrix$MEyellow[32:44],matrix$MEyellow[45:64],matrix$MEyellow[65:87], matrix$MEyellow[88:149],medcol="yellow", lwd=1.3, las=2, cex.axis=0.6, cex=.2, notch=F,xaxt='n')
abline(h=0,lty=2,col="black",lwd=0.8)
axis(side=1,at=c(1,2,3,4,5), cex.axis=0.6,labels=c("Depression","PTSD.Dep","PTSD","TE", "Control"),las=2)




stripchart(list(matrix$MEblue[1:31],matrix$MEblue[32:44],matrix$MEblue[45:64],matrix$MEblue[65:87], matrix$MEblue[88:149]),
            vertical = TRUE, method = "jitter", cex=0.2,
            pch = 21, col = "grey30", bg = "grey30",
            add = TRUE)

























#CREATE BOXPLOTS OF EXPDSSION VALUES FOR CANIDATE BIOMARKERS
ContrastBoxPlot <- file.path("Targets.txt")
cc1 <- as.matrix(read.table(ContrastBoxPlot, header=TRUE, sep="", row.names=1, as.is=TRUE))
head(cc1)

ContrastBoxPlot <- file.path("Targets.2.txt")
cc2 <- as.matrix(read.table(ContrastBoxPlot, header=TRUE, sep="", row.names=1, as.is=TRUE))
head(cc2)

pdf("BioMarkers_Post.pdf")
par(mar=c(8, 4, 4, 4) + 0.1)
boxplot(cc1, las=2, cex.axis=.9, main="PTSD Biomarkers (Cohort 1)", ylab="Average Expression",  col=(c("grey80","grey35")), outline=F, xlab="")
 legend("topleft", inset=0, title="",
   c("Healthy Control","PTSD"), fill=c("grey80","grey35"))
boxplot(cc2, las=2, cex.axis=.9, main="PTSD Biomarkers (Cohort 2)", ylab="Average Expression",  col=(c("grey80","grey35")), outline=F, xlab="")
 legend("topleft", inset=0, title="",
   c("Healthy Control","PTSD"), fill=c("grey80","grey35"))
dev.off()

mydata = read.table("PIN.txt", sep="\t")
par(mfrow = c(3,3))
par(oma=c(1,25,1,1)) 
barplot((mydata$V2), main="",horiz=TRUE, xlab="", log="x",
names.arg=mydata$V1,col=c("lightblue", "steelblue1","steelblue2","steelblue3","steelblue4"),
space=0.15, cex.axis=0.5, las=1,cex=1, xaxt="n", lend=2, xlim=c(1,1e-41),cex.main=0.5, font=8)
axis(3,at=c(1,1e-20,1e-41) ) # draw y axis with required labels
abline(v=5.00e-02,col="black", lty=2)

par(mfrow = c(3,2))
par(oma=c(1,20,1,1)) 
barplot((mydata$V2), main="Enrichment P-Value\n",horiz=TRUE, xlab="", log="x",
names.arg=mydata$V1,col=c("grey95","grey70","grey70","grey70","grey30","grey30","grey30"), space=0.15, cex.axis=0.9, las=1,cex=1, xaxt="n", lend=2, xlim=c(1,1e-30),cex.main=0.9, font=1)
axis(3,at=c(1,0 ,1e-15,1e-30) ) # draw y axis with required labels


par(oma=c(1,20,1,1)) 
barplot((mydata$V2), main="-log10 FDR P-Value\n",horiz=TRUE, xlab="", log="x",
names.arg=mydata$V1,col=c("grey95","grey95","grey95","grey70","grey70","grey70","grey70","grey70","grey70","grey30","grey30","grey30"), space=0.15, cex.axis=1, las=1,cex=1, xaxt="n", lend=2, xlim=c(1,3.22E-30),cex.main=0.9, font=3)
axis(3,at=c(1,1.0E-5,1.0E-10,1.0E-15 ,1.0E-20,1.0E-25) ) # draw y axis with required labels



par(oma=c(1,20,1,1)) 
barplot((mydata$V2), main="-log10 FDR P-Value\n",horiz=TRUE, xlab="", log="x",
names.arg=mydata$V1,col=c("white","white","white","white","white","white","black","black","black"), space=0.15, cex.axis=1, las=1,cex=1, xaxt="n", lend=2, xlim=c(1,3.22E-30),cex.main=0.9, font=3)
axis(3,at=c(1,1.0E-5,1.0E-10,1.0E-15 ,1.0E-20,1.0E-25) ) # draw y axis with required labels


barplot((mydata$V2), main="-log10 FDR P-Value\n",horiz=TRUE, xlab="", log="x",
names.arg=mydata$V1,col=c("grey30","grey30","grey30","grey30","grey30","grey30", "grey20","grey20","grey20","grey20","grey20","grey20"), space=0.1, cex.axis=1, las=1,cex=1, xaxt="n", lend=2, xlim=c(1,3.22E-30),cex.main=0.9, font=3)
axis(3,at=c(1,1.0E-05,1.0E-10,1.0E-15,1.0E-20 ,1.0E-25) ) # draw y axis with required labels

par(mfrow = c(4,4))
par(oma=c(1,6,1,1)) 
barplot((mydata2$V2), main="Cell-Type Enrichment\n", names.arg=mydata2$V1,horiz=TRUE, las=2,col="orange", space=0.09, cex.axis=1,cex=1, lend=2, cex.main=0.9,font=1, ylab="",lend=2,xlim=c(0,20),xaxt="n")
axis(3,at=c(0,5,10,15,20) ) # draw y axis with required label)
barplot((mydata2$V2), main="Cell-Type Enrichment\n", names.arg=mydata2$V1,horiz=TRUE, las=2,col="grey50", space=0.09, cex.axis=1,cex=1, lend=2, cex.main=0.9,font=1, ylab="",lend=2,xlim=c(0,20),xaxt="n")
axis(3,at=c(0,5,10,15,20) ) # draw y axis with required label)


barplot((mydata2$V2), main="-log10 Bonferroni P-Value\n",horiz=TRUE, xlab="", log="x",
names.arg=mydata2$V1,col="grey", space=0.1, cex.axis=1, las=1,cex=1, xaxt="n", lend=2, xlim=c(1,3.22E-19),cex.main=1,font=2)
axis(3,at=c(1.0,1e-05,1E-10) ) # draw y axis with required labels 
dev.off()

mydata = read.table("ids.txt", sep="\t", header=T)
mydata2 = read.table("new.txt", sep="\t", header=T)

pdf("Improvements.pdf")
par(mfrow = c(2,2))
par(oma=c(6,1,1,1))
barplot(mydata2$Number, names.arg=mydata2$Main, las=2, col="lightgreen", main="General Improvments for Level E", ylab="Number of PhD / PostDocs")
barplot(mydata$Number, names.arg=mydata$Main, las=2, col=c("lightgreen","lightgreen","lightgreen","lightgreen","lightgreen","purple","purple","purple","purple","purple","skyblue4","skyblue4","skyblue4","brown3","brown3","midnightblue","midnightblue","midnightblue","midnightblue"), main="Detailed Improvments for Level E", ylab="Number of PhD / PostDocs")
dev.off()



pdf("Pre_Post_Top20KME_Targets.pdf")
par(mfrow = c(2,2))
par(oma=c(1,1,1,1))
ranged <- c(floor(min(mydata$MM)), ceiling(max(mydata$MM)))
ranged2 <- c(floor(min(mydata2$MM)), ceiling(max(mydata2$MM)))

barplot((mydata$Main), main="logFC\n", xlab="",
names.arg=mydata$Number,col="grey", cex.axis=0.6, las=1,cex=0.6, xaxt="n", lend=1, xlim=ranged,cex.main=1, font=2)
axis(3,at=c(2.5,2.0,1.5, 1.0, 0.5, 0) ) # draw y axis with required labels
abline(v=0.5, lty=2)


barplot((mydata$Pvalue), main="-log10 P-Value\n",horiz=TRUE, xlab="", log="x",
names.arg=mydata$Gene,col="grey", space=0.1, cex.axis=0.6, las=1,cex=0.6, xaxt="n", lend=2, xlim=c(1,3.22E-30),cex.main=1, font=2)
axis(3,at=c(1.0,1.0E-05,1.0E-10,1.0E-15, 1.0E-20) ) # draw y axis with required labels


barplot((mydata2$MM), main="logFC\n",horiz=TRUE, xlab="",
names.arg=mydata2$Gene,col="grey", cex.axis=0.6, las=1,cex=0.6, xaxt="n", lend=1, xlim=ranged2,cex.main=1, font=2)
axis(3,at=c(2.5,2.0,1.5, 1.0, 0.5, 0) ) # draw y axis with required labels
abline(v=0.5, lty=2)

barplot((mydata2$Pvalue), main="-log10 P-Value\n",horiz=TRUE, xlab="", log="x",
names.arg=mydata2$Gene,col="grey", space=0.1, cex.axis=0.6, las=1,cex=0.6, xaxt="n", lend=2, xlim=c(1,3.22E-30),cex.main=1, font=2)
axis(3,at=c(1.0,1.0E-05,1.0E-10,1.0E-15, 1.0E-20) ) # draw y axis with required labels

dev.off()

#NOW Microarray Data!!!
##############################################################################################################################
######### HERE WE MAP MODULES from RNA-Seq to MicroArray Data, AGAIN INCLUDING CASES AND CONTROLS ##############
###############################################################################################################################
#INPUT MICROARRAY DATA
ControlPre <- read.delim("MRSI_PreDeployment.txt.2.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(ControlPre) =names(ControlPre)
dataExpr1<-as.data.frame(t(ControlPre))
head(dataExpr1)
dim(dataExpr1)
gsg=goodSamplesGenes(dataExpr1,verbose=3)
gsg$allOK

nGenes2 = ncol(dataExpr1)
nSamples2 = nrow(dataExpr1)

traitData2 <- read.delim("MRSI_Clinical_Traits_Pre.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData2)

#MATCH TRAITS TO PD-DEPLOYMENT
rowsExpr1 <- rownames(dataExpr1)
traitRows1 <- match(rowsExpr1,traitData2$studyid.1)
datTraits1 = traitData2[traitRows1, -1];
rownames (datTraits1) = traitData2[traitRows1, 1];
table(rownames(datTraits1)==rownames(dataExpr1)) # EVERYTHING OK?
names(datTraits1)

powers=c(1:40) # in practice this should include powers up to 20.
sft0=pickSoftThreshold(dataExpr1,powerVector=powers, networkType="signed")

pdf("MRSI_PreDeploy_Soft_Threshold.pdf")
par(mfrow=c(1,2))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],labels=powers,col="red")
abline(h=0.90,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

adjacencyMRSIPre = adjacency((dataExpr1),power=30 ,type="signed"); #FOR ONE LARGE NETWORK WE CHOOSE A POWER OF 10
diag(adjacencyMRSIPre)=0
dissTOMMRSIPre   = 1-TOMsimilarity(adjacencyMRSIPre, TOMType="signed")
geneTreeMRSIPre  = flashClust(as.dist(dissTOMMRSIPre), method="average")

pdf("MRSI_PreDeploy_GeneTree.pdf",height=6,width=12)
par(mfrow=c(1,1))
plot(geneTreeMRSIPre,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (Cohort 2 Pre-Deployment)", labels=FALSE,hang=0.04);
dev.off()

#LETS CONSTRUCT MODULES AND USE GENE SIGNFICANCE MEASURE TO ASSOCIATE TO DISEASE AND HEALTHY
mColorh=NULL
for (ds in 0:4){
 tree = cutreeHybrid(dendro = geneTreeMRSIPre, pamStage=FALSE,
   minClusterSize = (60-3*ds), cutHeight = 0.995, 
   deepSplit = ds, distM = dissTOMMRSIPre)
 mColorh=cbind(mColorh,labels2colors(tree$labels));
}
pdf("MRSI_PreDeploy_DeepSplit_Choices_60_2.pdf", height=10,width=25); 
plotDendroAndColors(geneTreeMRSIPre, mColorh, paste("dpSplt =",0:4), main = "Gene Dendrogram and Module Colors - Cohort 2",dendroLabels=FALSE);
dev.off()

#SET DEEP SPLIT CHOICE AND NAME OUR COLORS
modulesMRSIPD =  mColorh[,1]
table(modulesMRSIPD)

#Check to see if disease network modules can be cut and merged...
MEList = moduleEigengenes(dataExpr1, colors=modulesMRSIPD)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = flashClust(as.dist(MEDiss), method ="average")

pdf("MRSI_PreDeploy_Module_Relationships_2.pdf")
plot(METree, main ="Clustering of Cohort 2 Pre-Deployment Module Eigengenes",xlab ="",sub="",mar=c(5.1, 4.1, 4.1, 3))
MEDissThres = 0.1
abline(h=MEDissThres,col="red")
dev.off()

####################################### SHOULD WE MERGE MODULES BASED ON THE ABOVE? #######################################
merge = mergeCloseModules(dataExpr1, modulesMRSIPD, cutHeight = MEDissThres, verbose=3)
mergedColors = merge$colors
table(mergedColors)
table(modulesMRSIPD)
mergedMEs = merge$newME

#COMPARE UNMERGED MODULES TO MERGED MODULES
pdf("PreDeploy_Network_Unmerged_Merged.pdf", w=9)
plotDendroAndColors(geneTreeMRSIPre, cbind(modulesMRSIPD, mergedColors),c("Dynamic Tree Cut","Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

modulesMRSIPD = merge$colors
table(modulesMRSIPD)
MEList = moduleEigengenes(dataExpr1, colors=modulesMRSIPD)
MEs=MEList$eigengenes
##########################################################################################################################
##Our Options
# deep split = 0
# Min Module Size = 60
# Merge past = 0.1
# cutHeight = 0.995
#POWER 20

#CALCULATE PC FOR VISUALIZATION FOR CASE PD-DEPLOYMENT
PCsPD    = moduleEigengenes((dataExpr1),  colors=modulesMRSIPD) 
ME_PD    = PCsPD$eigengenes
distPCPD = 1-abs(cor(ME_PD,use="p"))
distPCPD = ifelse(is.na(distPCPD), 0, distPCPD)
pcTreePD = hclust(as.dist(distPCPD),method=" ") 
MDS_PD   = cmdscale(as.dist(distPCPD),2)
colorsPD = names(table(modulesMRSIPD))
names= row.names((dataExpr1))

pdf("MRSI_PreDeploy_Module_Visualization.pdf",height=8,width=8)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 3) + 0.1, cex=1)
plot(pcTreePD, xlab="",ylab="",main="",sub="")
plot(MDS_PD, col= colorsPD,  main="MDS plot", cex=2, pch=19)

for (which.module in names(table(modulesMRSIPD)))
{
 par(mfrow=c(2,1), mar=c(4, 4.1, 4.1, 2))
 plotMat(t(scale(dataExpr1[,modulesMRSIPD==which.module])),
,cex.axis=2,nrgcols=100,rlabels=F,tck=0, rcols=which.module,main=paste("Heatmap",which.module,"Module"))

  ME = ME_PD[, paste("ME",which.module, sep="")] 
  n<- barplot(ME, col=which.module, cex.main=1, ylab="Eigengene Expression",xlab="")
  axis(1,at=n, labels=row.names(dataExpr1), las=2, cex.axis=0.5, font=2)
};
dev.off();

names(datTraits1)
#DEFINE GENE-SIGNIFICANCE MEASURE FOR DISEASE STATUS
Eventual_PTSD = as.data.frame(datTraits1$Eventual_PTSD)
names(Eventual_PTSD)="Eventual_PTSD"
GS.Eventual_PTSD=as.numeric(cor(dataExpr1,Eventual_PTSD,use="p"))
GS.Eventual_PTSD_abs=abs(GS.Eventual_PTSD)
GS.Eventual_PTSDColor=numbers2colors(GS.Eventual_PTSD,signed=T)

Healthy_Control = as.data.frame(datTraits1$Healthy_Control)
names(Healthy_Control)="Healthy_Control"
GS.Healthy_Control=as.numeric(cor(dataExpr1,Healthy_Control,use="p"))
GS.Healthy_ControlColor=numbers2colors(GS.Healthy_Control,signed=T)

#DETERMINE GENE-SIG. ACROSS MODULES
pdf("MRSI_PreDeploy_Module_GS.pdf")
plotModuleSignificance(GS.Eventual_PTSD, modulesMRSIPD, boxplot = FALSE, main = "Cohort 2 Pre-Deployment\nModule Significance (cor)\n", ylab = "M (cor, Eventual PTSD)", las=2, cex.axis=0.7)
plotModuleSignificance(GS.Eventual_PTSD, modulesMRSIPD, boxplot = TRUE, main = "Cohort 2 Pre-Deployment\nModule Significance (cor)\n", ylab = "MS (cor, Eventual PTSD)", las=2, cex.axis=0.7)
plotModuleSignificance(GS.Eventual_PTSD_abs, modulesMRSIPD, boxplot = FALSE, main = "Cohort 2 Pre-Deployment\nModule Significance (cor)\n", ylab = "MS |(cor, Eventual PTSD)|", las=2, cex.axis=0.7)
dev.off()

#DETERMINE MODULE SIG. USING -LOG PVALUE
datSummary1 <- read.delim("MRSI_PreDeploy_DatSummary.txt.2.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
GS_P=-log10(datSummary1$PValue)
GS_P_Color=numbers2colors(GS_P,signed=T)


pdf("MRSI_PreDeploy_Module_MS_white.pdf")
par(mfrow=c(2,1))
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = FALSE)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H
par(mfrow=c(2,1))
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = FALSE)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H

par(mfrow=c(2,2))
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = FALSE)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = FALSE)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H

par(mfrow=c(3,2))
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = FALSE)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = FALSE)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H

par(mfrow=c(3,3))
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = FALSE)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H
verboseBarplot(GS_P,modulesMRSIPD,color="white",main="Dataset 2 Pre-Deployment Modules\n" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9,KruskalTest = FALSE)
abline(h=1.0,col="black",lty = 2)    #CHOOSE A  R^2 CUT-OFF OF H

dev.off()

##CELL-Type Enrichment Per Module
Gene  = colnames(dataExpr1)

enrichments = userListEnrichment(Gene, modulesMRSIPD,
fnIn = NULL,
catNmIn = NULL,
#fnIn = c("GeneList","ModuleColors"),
#catNmIn =  c("Genes","Modules"),
nameOut = "MRSI_PreDeployment_Cell_Enrichment.csv", 
useBrainLists = TRUE,
useBloodAtlases = TRUE,
omitCategories = "grey", 
outputCorrectedPvalues = TRUE,
useStemCellLists = FALSE, 
outputGenes = TRUE, 
minGenesInCategory = 5, 
useBrainRegionMarkers = FALSE,
useImmunePathwayLists = TRUE,
usePalazzoloWang = TRUE)

# More detailed overlap information is in the pValue output. 
head(enrichments$pValue)

# To see a list of all significant enrichments
enrichments$sigOverlaps
 
# To see all of the overlapping genes between two categories  
enrichments$ovGenes$"black -- IFN alpha/beta__ImmunePathway"
enrichments$ovGenes$"black -- IFN-gamma__ImmunePathway"

IFN_A_B = enrichments$ovGenes$"salmon -- IFN alpha/beta__ImmunePathway"
IFN_G =enrichments$ovGenes$"salmon -- IFN-gamma__ImmunePathway"
colorCT = rep("grey",length(modulesMRSIPD))
colorCT[is.element(Gene,IFN_A_B)] = "red";
colorCT[is.element(Gene,IFN_G)] = "blue";

MEsCT = (moduleEigengenes(dataExpr1, colors=as.character(colorCT), excludeGrey=TRUE))$eigengenes
datCT = t(cor(MEsCT,dataExpr1))
colnames(datCT) = c("INF A_B","INFG")

modsCT  = c("yellow","turquoise")
colorCT = rep("grey",length(colorCT))
for (i in 2:1){  # Use the top 2% of t-scores to determine colors 
 values = as.numeric(datCT[,i])
 colorCT[values>quantile(values,0.95)] = modsCT[i]
}

#INTRODUCE DATAFRAME COLORS USING CORRELATIONS AND PVALUES
datColors0=data.frame(modulesMRSIPD,GS.Eventual_PTSDColor,GS_P_Color,colorCT)
datColors1=data.frame(modulesMRSIPD,GS.Eventual_PTSDColor)

pdf("MRSI_PreDeploy_FinalModule.pdf",height=8,width=14)
plotDendroAndColors(geneTreeMRSIPre, colors=datColors0, main="Cohort 2 Pre-Deployment\nGene Dendrogram and Module Colors (Beta=20)", groupLabels=c("Module colors","GS.Eventual PTSD\n(cor)", "GS.Eventual PTSD\n(-log10 p)", "INF A/B"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
#plotDendroAndColors(geneTreeMRSIPre, colors=datColors1, main="Cohort 2 Pre-Deployment\nGene Dendrogram and Module Colors (Beta=20)", groupLabels=c("Module colors","GS.Eventual PTSD\n(cor)"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
dev.off()


#PLOT RELATIONS AMONG EIGENGENES AND THE TRAITS OF INTEREST
MET=orderMEs(cbind(MEs,Eventual_PTSD,Healthy_Control))
pdf("MRSI_PreDeploy_Modules_Trait.pdf", h=16, w=15)
plotEigengeneNetworks(MET,"",marDendro=c(5.1, 4.1, 4.1, 3), marHeatmap=c(6,6,4,4),cex.lab=0.8,xLabelsAngle=90)
dev.off()

# PLOT MODULE-TRAIT RELATIONSHIP
MEs0  =  moduleEigengenes(dataExpr1,  modulesMRSIPD)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits1, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples2)
data = cbind(moduleTraitCor,moduleTraitPvalue)
write.table(data, "MRSI_Pre_Correlations.txt", sep="\t")
pdf("MRSI_PreDeploy_Correlations_2.pdf")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits1),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.3,
zlim = c(-1,1),
main = paste("Cohort 2 Pre-Deployment\n ME-Trait Relationships"))
dev.off()

#PLOT CONNECTIVITY VS GS
CM_Pre=intramodularConnectivity(adjacencyMRSIPre,colors=modulesMRSIPD)
names(CM_Pre) 
whichmodule="salmon";
restrict1=modulesMRSIPD==whichmodule
pdf("MRSI_salmon_Connectivity.pdf")
verboseScatterplot (CM_Pre$kWithin[restrict1],GS.Eventual_PTSD[restrict1],col=modulesMRSIPD[restrict1],
xlab="Connectivity (k) - C2 ",ylab="G.S. Eventual PTSD - C2") 
dev.off()

# CALCULATE MODULE MEMBERSHIP VALUES (aka. module eigengene based connectivity kME)
datKME=signedKME(dataExpr1, MEs)
colorOfColumn=substring(names(datKME),4)

pdf("MRSI_PreDeploy_Regress_Modules.pdf", h=10, w=16)
par(mfrow = c(3,4))
par(mfrow=c(4,length(selectModules)/4))
 for (module in names(table(modulesMRSIPD))) 
  {
  column = match(module,colorOfColumn) 
  restModule=modulesMRSIPD==module
  verboseScatterplot(datKME[restModule,column],GS.Eventual_PTSD[restModule],
  xlab=paste("Module Membership -C2-",module,"module"),ylab="GS.PTSD Risk",
  main=paste("kME",module,"vs GS"),col="black")
  
  verboseScatterplot(datKME[restModule,column],GS.Healthy_Control[restModule],
  xlab=paste("Module Membership -C2-",module,"module"),ylab="GS.Control",
  main=paste("kME",module,"vs GS"),col="black")
  }
dev.off()

#MODULE MEMBERSHIP (kME) KME is defined as correlation between expression and modules
#USED TO MEASURE CORRELATIONS BETWEEN EACH GENE AND EACH MODULE EIGENGENE

geneModuleMembership1 = signedKME((dataExpr1), MEs)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(dataExpr1)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");

Gene       = rownames(t(dataExpr1))
kMEtable1  = cbind(Gene,Gene,modulesMRSIPD)
for (i in 1:length(colorsPD))
 kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))
write.csv(kMEtable1,"MRSI_PreDeploy_kMEtable.csv",row.names=FALSE)

#WRITE OUTPUT
Pre<-as.data.frame(datTraits1)
names(Pre)<-"Pre"
modNames = substring(names(MEs), 3)
nGenes = ncol(dataExpr1) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
nSamples = nrow(dataExpr1) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
geneModuleMembership<-as.data.frame(cor(dataExpr1,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr1, Pre, use = "p")) #  Generate correlations and p-value for each gene against the trait.  
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Pre), sep="")
names(GSPvalue) = paste("p.GS.", names(Pre), sep="")

Colors=modulesMRSIPD
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="MRSI_PreDeploy_Network_WGCNA_OUTPUT.csv")


#HERE WE EXPORT 20 TOP INTRAMODULAR HUB GENES PER MODULE OF INTEREST
modules = c("blue") #"red", "yellow", "lightgreen") #,"yellow", "lightgreen", "black"
probes=names(dataExpr1)
inModule=is.finite(match(modulesMRSIPD, modules));
modProbes=probes[inModule] ##
modTom=dissTOMMRSIPre[inModule, inModule]
dimnames(modTom)=list(modProbes, modProbes)
nTopHubs = 50
kIN = softConnectivity(dataExpr1[, modProbes])
selectHubs = (rank (-kIN) <= nTopHubs)


#EXPORT ENTIRE NETWORK
cyt=exportNetworkToCytoscape(modTom,
edgeFile=paste("CytoscapeInput_Pre_.05_", paste(modules,collapose="_") , "edges.txt", sep=""),
nodeFile=paste("CytoscapeInput_Pre_.05_", paste(modules,collapose="_") , "nodes.txt", sep=""),
weighted=TRUE,
threshold=0.50,
nodeNames=modProbes,
nodeAttr=modulesMRSIPD[inModule])










library(RColorBrewer)
library(gplots)
library(ClassDiscovery)

exprsFile <- file.path("Dataset2_Post_order.txt")
exprs1 <- as.matrix(read.table(exprsFile, header=TRUE, sep="", row.names=1, as.is=TRUE))
dim(exprs1)

exprsFile <- file.path("Red_Module_Heat2.txt")
exprs2 <- as.matrix(read.table(exprsFile, header=TRUE, sep="", row.names=1, as.is=TRUE))
dim(exprs2)

mycolors <-rep(c( "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black"), times =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

mycolors <-rep(c( "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "grey80", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black"), times=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

dist.pear <- function(x) as.dist(1-cor(t(x)))
dist.euc <- function(x) dist(x,method = "euclidean")
hclust.ave <- function(x) hclust(x, method="average")
hclust.ward <- function(x) hclust(x, method="ward")

pdf("Red_HEATMAP.pdf")
heatmap.2(exprs1, col=bluered(75), scale="row",ylab="70 Genes", key=F, symkey=FALSE, density.info="hist", trace="none", cexRow=0.005, ColSideColors=mycolors,
main="M4A", distfun=dist.euc, hclustfun=hclust.ward, margins=c(8, 10), symbreaks=T, cex=0.008, Colv=F)
myCol <- c("grey80", "black")
legend("top",horiz=TRUE, fill = myCol, cex=0.9,legend = c("Case", "Control"), border=FALSE, bty="n")
dev.off()


mydata = read.table("ids.txt", sep="\t")
mydata2 = read.table("MRSI_Pre.txt", sep="\t")
pdf("Structure_Ontology_Figures.pdf")
par(mfrow = c(2,1))
par(oma=c(1,11,1,1)) 
barplot((mydata$V2), main="-log10 Bonferroni P-Value\n",horiz=TRUE, xlab="", log="x",
names.arg=mydata$V1,col="grey", space=0.1, cex.axis=0.8, las=1,cex=0.8, xaxt="n", lend=2, xlim=c(1,3.22E-19),cex.main=1, font=2)
axis(3,at=c(1.0,1.0E-05,1.0E-10,1.0E-15) ) # draw y axis with required labels

barplot((mydata2$V2), main="-log10 Bonferroni P-Value\n",horiz=TRUE, xlab="", log="x",
names.arg=mydata2$V1,col="grey", space=0.1, cex.axis=0.8, las=1,cex=0.8, xaxt="n", lend=2, xlim=c(1,3.22E-19),cex.main=1,font=2)
axis(3,at=c(1.0,1e-05,1E-10,1E-15) ) # draw y axis with required labels 
dev.off()

########################################################################################################################
#########################################Here we attempt Preservation Statistics########################################
########################################################################################################################


pdf("MRSI_PreDeploy_Colored_Modules.pdf",height=8,width=12)
par(mfrow=c(2,2))
plotDendroAndColors(geneTreePre, colors=modulesPRE, main="Control Pre-Deployment Coloured with Control Pre-Deploy Modules", groupLabels=c("Module Colors"),dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
plotDendroAndColors(geneTreeMRSIPre, colors=modulesPRE, main="Case Pre-Deployment Coloured with Control Pre-Deploy Modules",groupLabels=c("Module Colors"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
dev.off()

#DEFINE GENE-SIGNIFICANCE MEASURE FOR DISEASE STATUS
Eventual_PTSD_A2 = as.data.frame(datTraits1$Eventual_PTSD)
names(Eventual_PTSD_A2)="Eventual_PTSD_A2"
GS.Eventual_PTSD_A2=as.numeric(cor(dataExpr1,Eventual_PTSD_A2,use="p"))
GS.Eventual_PTSD_A2Color=numbers2colors(GS.Eventual_PTSD_A2,signed=T)

Healthy_Control_A2 = as.data.frame(datTraits1$Healthy_Control)
names(Healthy_Control_A2)="Healthy_Control_A2"
GS.Healthy_Control_A2=as.numeric(cor(dataExpr1,Healthy_Control_A2,use="p"))
GS.Healthy_Control_A2Color=numbers2colors(GS.Healthy_Control_A2,signed=T)
datColors0=data.frame(modulesPRE,GS.Eventual_PTSD_A2Color,GS.Healthy_Control_A2Color)

pdf("PreDeploy_FinalModule.pdf",height=8,width=14)
plotDendroAndColors(geneTreePre, colors=datColors0, main="Pre-Deployment Gene Dendrogram and Module Colors", groupLabels=c("Module colors","GS.Eventual PTSD","GS.Healthy Cntl"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
dev.off()

#Assess Module Preservation to quanitfy our visualization results
multiExpr  = list(CASEPD=list(data=(dataExpr0)),CONTROLPD=list(data=(dataExpr1)))
multiColor = list(CASEPD = modulesMRSIPD)
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed", nPermutations=30,maxGoldModuleSize=100,maxModuleSize=1000)
stats = mp$preservation$Z$ref.CASEPD$inColumnsAlsoPresentIn.CONTROLPD
stats[order(-stats[,2]),c(1:2)]

modColors=rownames(stats)
moduleSize = stats$moduleSize
selectModules = !(modColors %in% c("grey", "gold"))
point.label = modColors[selectModules]

#Composite preservation statistics
Zsummary=stats$Zsummary.pres

# plot Zsummary versus module size
pdf("Zsummary_Post2Pre.pdf")
plot(moduleSize[selectModules],Zsummary[selectModules], col = 1,
bg=modColors[selectModules],pch = 21,main="Zsummary Preservation\nPost-Deployment Modules in Pre-Deployment Samples",
cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1,offs=0.03)
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)
dev.off()


#Find kME values for Cases Pre-Deployment
geneModuleMembership1 = signedKME((dataExpr1), ME_PD)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(dataExpr1)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");
Gene= rownames(t(dataExpr1))
kMEtable1 = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
             kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))
write.csv(kMEtable1,"Control_Pre_kMEtable_for_ControlModules.csv",row.names=FALSE)

#Find kME Values for Controls Pre-Deployment
# First calculate MEs for A2, since we haven't done that yet
PCsAD = moduleEigengenes((dataExpr0),  colors=modulesPRE) 
MEsAD = PCsAD$eigengenes
geneModuleMembership2 = signedKME((dataExpr0), MEsAD)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 
MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(dataExpr0)[[2]]); 
colnames(MMPvalue2)=paste("PC",colorsPD,".pval",sep="");
kMEtable2  = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
 kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i])
colnames(kMEtable2)=colnames(kMEtable1)
write.csv(kMEtable2,"CasePre_kMEtable_for_ControlModules.csv",row.names=FALSE)


pdf("ControlPre_vs_CasePre_kMEtable_for_ControlModules.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsPD[c],
                    xlab="kME in Control Pre-Deploy",ylab="kME in Case Pre-Deploy")
}; dev.off()

#ASSESS HUB GENE CONSERVATION
pdf("inModule_ControlPre_vs_CasePre_kMEtable_for_ControlModules.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 inMod = modulesPRE== colorsPD[c]
 verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colorsPD[c],
                    xlab="kME in Control Pre-Deploy",ylab="Case Pre-Deploy")
}; dev.off()


#DETERMINE WHICH GENES ARE HUB GENES IN BOTH NETWORKS (GENES WITH HIGH KME VALUES IN BOTH NETWORKS), TOP 10 GENES PER MODULE BASED ON KME IN BOTH NETWORKS
topGenesKME = NULL
for (c in 1:length(colorsPD)){
 kMErank1    = rank(-geneModuleMembership1[,c])
 kMErank2    = rank(-geneModuleMembership2[,c])
 maxKMErank  = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
 topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; colnames(topGenesKME) = colorsPD
#PRINT TOP 10 GENES PER MODULE BASED ON KME IN BOTH NETWORKS
topGenesKME


# PLOT MODULE-TRAIT RELATIONSHIP
MEs0  =  moduleEigengenes(dataExpr0,  modulesPRE)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf("CasePre_Correlations.pdf")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.3,
zlim = c(-1,1),
main = paste("Case Pre-Deployment \n Module Eigengene Expression-Trait Relationships"))
dev.off()

#GENERATE One Network OUTPUT
Pre<-as.data.frame(datTraits1)
names(Pre)<-"Pre"
modNames = substring(names(MEs), 3)
nGenes = ncol(dataExpr1) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
nSamples = nrow(dataExpr1) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
geneModuleMembership<-as.data.frame(cor(dataExpr1,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr1, Pre, use = "p")) #  Generate correlations and p-value for each gene against the trait.  
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Pre), sep="")
names(GSPvalue) = paste("p.GS.", names(Pre), sep="")

#OUTPUT SPDADSHEET INFORMATION
Colors=modulesPRE
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="ControlPre_Network_WGCNA_OUTPUT.csv")

#Export to Cytoscape
modules = c("grey60","purple")
# Select module probes
inModule=is.finite(match(modulesPRE,modules))
Genes = row.names(t(dataExpr0))
modGenes=Genes[inModule]
match1=match(modGenes,GeneAnnotation$substanceBXH)
modGenes=GeneAnnotation$gene_symbol[match1]

modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files for Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
edgeFile=paste("CytoEdge",paste(modules,collapse="-"),".txt",sep=""),
nodeFile=paste("CytoNode",paste(modules,collapse="-"),".txt",sep=""),
weighted = TRUE, threshold = 0.02,nodeNames=modProbes,
altNodeNames = modGenes, nodeAttr = moduleColorsFemale[inModule])















#NOW MICROARRAY DATA!!!!
##############################################################################################################################
######### HERE WE CREATE ONE LARGE NETWORK INCLUDING CASES AND CONTROLS AND DETERMINE PDSERVATION ACROSS COHORTS  ##########
###############################################################################################################################
library(WGCNA)
library(cluster)
options(stringsAsFactors  =  FALSE)
allowWGCNAThreads()

#INPUT CASE PD
PD <- read.delim("MRSI_PreDeploy.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(PD) =names(PD)
dataExpr0<-as.data.frame(t(PD))
head(dataExpr0)
dim(dataExpr0)
gsg=goodSamplesGenes(dataExpr0,verbose=3)
gsg$allOK

nGenes = ncol(dataExpr0)
nSamples = nrow(dataExpr0)

traitData <- read.delim("MRSI_Clinical_Traits_Pre.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData)

#MATCH TRAITS TO PD-DEPLOYMENT
rowsExpr <- rownames(dataExpr0)
traitRows <- match(rowsExpr,traitData$study_id.1)
datTraits = traitData[traitRows, -1];
rownames (datTraits) = traitData[traitRows, 1];
table(rownames(datTraits)==rownames(dataExpr0)) # EVERYTHING OK?
names(datTraits)

powers=c(1:40) # in practice this should include powers up to 20.
sft0=pickSoftThreshold(dataExpr0,powerVector=powers, networkType="signed")

pdf("MRSI_PreDeploy_Soft_Threshold.pdf")
par(mfrow=c(1,2))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],labels=powers,col="red")
abline(h=0.90,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

adjacencyPre = adjacency((dataExpr0),power=20 ,type="signed"); # WE CHOOSE A POWER OF 20
diag(adjacencyPre)=0
dissTOMPre   = 1-TOMsimilarity(adjacencyPre, TOMType="signed")
geneTreePre  = flashClust(as.dist(dissTOMPre), method="average")

pdf("MRSI_PreDeploy_GeneTree.pdf",height=6,width=12)
par(mfrow=c(1,1))
plot(geneTreePre,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (Pre-Deployment) A2", labels=FALSE,hang=0.04);
dev.off()

cmd1=cmdscale(as.dist(dissTOMPre),2) #This may take a while... 
pdf("MRSI_PreDeploy_MDS_Plot.pdf")
par(mfrow=c(1,1))
plot(cmd1,col=moduleColorsFemale,main="MDS plot", xlab="Scaling Dimension 1",ylab="Scaling Dimension 2")
dev.off()

#CONSTRUCT MODULES AND CHOOSE BEST DEEP-SPLIT CUT-OFF OPTION
mColorh=NULL
for (ds in 0:3){
 tree = cutreeHybrid(dendro = geneTreePre, pamStage=FALSE,
   minClusterSize = (20-3*ds), cutHeight = 0.99998,
   deepSplit = ds, distM = dissTOMCasePre)
 mColorh=cbind(mColorh,labels2colors(tree$labels));
}

pdf("MRSI_PreDeploy_DeepSplit_Choices_20_998.pdf", height=10,width=25); 
plotDendroAndColors(geneTreePre, mColorh, paste("dpSplt =",0:3), main = "Pre-Deployment Gene Dendrogram - A2",dendroLabels=FALSE);
dev.off()

#SET DEEP SPLIT CHOICE AND NAME OUR COLORS
modulesPRE =  mColorh[,3]
table(modulesPRE)

#CHECK TO SEE IF MODULES CAN BE MERGED BASED ON SIMILARITY (HEIGHT)
MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = flashClust(as.dist(MEDiss), method ="average")

pdf("MRSI_PreDeploy_Module_Relationships.pdf")
plot(METree, main ="Clustering of Pre-Deployment Module Eigengenes - A2",xlab ="",sub="")
MEDissThres = 0.1
abline(h=MEDissThres,col="red")
dev.off()
####################################### SHOULD WE MERGE MODULES BASED ON THE ABOVE? #######################################
merge = mergeCloseModules(dataExpr0, modulesPRE, cutHeight = MEDissThres, verbose=3)
mergedColors = merge$colors
table(mergedColors)
table(modulesPRE)
mergedMEs = merge$newME

#COMPARE UNMERGED MODULES TO MERGED MODULES
pdf("MRSI_PreDeploy_Network_Unmerged_Merged.pdf", w=9)
plotDendroAndColors(geneTreeControlPre, cbind(modulesPRE, mergedColors),c("Dynamic Tree Cut","Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

modulesPRE = merge$colors
table(modulesPRE)
MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes
##########################################################################################################################
##Our Options
# deep split = 2
# Min Module Size = 30
# Merge past = 0.1


#CALCULATE PC FOR VISUALIZATION FOR CASE PD-DEPLOYMENT
PCsPD    = moduleEigengenes((dataExpr0),  colors=modulesPRE) 
ME_PD    = PCsPD$eigengenes
distPCPD = 1-abs(cor(ME_PD,use="p"))
distPCPD = ifelse(is.na(distPCPD), 0, distPCPD)
pcTreePD = hclust(as.dist(distPCPD),method="average") 
MDS_PD   = cmdscale(as.dist(distPCPD),2)
colorsPD = names(table(modulesPRE))
names= row.names((dataExpr0))

pdf("MRSI_PreDeploy_Module_Visualization.pdf",height=8,width=8)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 3) + 0.1, cex=1)
plot(pcTreePD, xlab="",ylab="",main="",sub="")
plot(MDS_PD, col= colorsPD,  main="MDS plot", cex=2, pch=19)

ordergenes = geneTreePD$order
plot.mat(scale(log(dataExpr0[ordergenes,])) ,
rlabels= modulesPRE[ordergenes], clabels= colnames(dataExpr0), rcols=modulesPRE[ordergenes])

for (which.module in names(table(modulesPRE)))
{ 
  ME = ME_PD[, paste("ME",which.module, sep="")] 
  n<- barplot(ME, col=which.module, main=paste(which.module, "module"), cex.main=1, ylab="eigengene expression",xlab="Pre-Deployment Samples")
  axis(1,at=n, labels=row.names(dataExpr0), las=2, cex.axis=0.5, font=2)
    axis(1,at=n[c(0, 49)], las=1, cex.axis=0.005, font=2, tck=1, col.ticks = 'darkred')
    #axis(1,at=n[c(1, 48)], las=1, font=2, tck=1)
};
dev.off();


Eventual_PTSD = as.data.frame(datTraits$Eventual_PTSD)
names(Eventual_PTSD)="Eventual_PTSD"
GS.Eventual_PTSD=as.numeric(cor(dataExpr0,Eventual_PTSD,use="p"))
GS.Eventual_PTSDColor=numbers2colors(GS.Eventual_PTSD,signed=T)

Healthy_Control = as.data.frame(datTraits$Healthy_Control)
names(Healthy_Control)="Healthy_Control"
GS.Healthy_Control=as.numeric(cor(dataExpr0,Healthy_Control,use="p"))
GS.Healthy_ControlColor=numbers2colors(GS.Healthy_Control,signed=T)

datColors0=data.frame(modulesPRE,GS.Eventual_PTSDColor,GS.Healthy_ControlColor)

pdf("PreDeploy_FinalModule.pdf",height=8,width=14)
plotDendroAndColors(geneTreePre, colors=datColors0, main="Pre-Deployment Gene Dendrogram and Module Colors", groupLabels=c("Module colors","GS.Eventual PTSD","GS.Healthy Cntl"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
dev.off()


#PLOT RELATIONS AMONG EIGENGENES AND THE TRAITS OF INTEREST
MET=orderMEs(cbind(MEs,Eventual_PTSD,Healthy_Control))
pdf("PreDeploy_Modules_Trait.pdf", h=16, w=15)
plotEigengeneNetworks(MET,"",marDendro=c(1,4,1,2), marHeatmap=c(6,6,4,4),cex.lab=0.8,xLabelsAngle=90)
dev.off()

# PLOT MODULE-TRAIT RELATIONSHIP
MEs0  =  moduleEigengenes(dataExpr0,  modulesPRE)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf("PreDeploy_Correlations.pdf")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.3,
zlim = c(-1,1),
main = paste("Pre-Deployment \n Module Eigengene Expression-Trait Relationships"))
dev.off()


# CALCULATE MODULE MEMBERSHIP VALUES (aka. module eigengene based connectivity kME)
datKME=signedKME(dataExpr0, MEs)
colorOfColumn=substring(names(datKME),4)

pdf("PreDeploy_Regress_Modules.pdf", h=10, w=14)
par(mfrow = c(3,2))
par(mfrow=c(4,length(selectModules)/4))
 for (module in names(table(modulesPRE))) 
  {
  column = match(module,colorOfColumn) 
  restModule=modulesPRE==module
  verboseScatterplot(datKME[restModule,column],GS.Eventual_PTSD[restModule],
  xlab=paste("Module Membership ",module,"module"),ylab="GS.Eventual_PTSD",
  main=paste("kME",module,"vs GS"),col=module)
  
  verboseScatterplot(datKME[restModule,column],GS.Healthy_Control[restModule],
  xlab=paste("Module Membership ",module,"module"),ylab="GS.Healthy_Control",
  main=paste("kME",module,"vs GS"),col=module)
  }
dev.off()

#MODULE MEMBERSHIP (kME) TO COMPARE NETWORKS
#USED TO MEASURE CORRELATIONS BETWEEN EACH GENE AND EACH MODULE EIGENGENE
#PD
geneModuleMembership1 = signedKME((dataExpr0), MEs)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 

MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(dataExpr0)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");

Gene       = rownames(t(dataExpr0))
kMEtable1  = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
 kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))
write.csv(kMEtable1,"PreDeploy_kMEtable.csv",row.names=FALSE)

Pre<-as.data.frame(datTraits)
names(Pre)<-"Pre"
modNames = substring(names(MEs), 3)
nGenes = ncol(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
nSamples = nrow(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
geneModuleMembership<-as.data.frame(cor(dataExpr0,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr0, Pre, use = "p")) #  Generate correlations and p-value for each gene against the trait.  
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Pre), sep="")
names(GSPvalue) = paste("p.GS.", names(Pre), sep="")

#OUTPUT SPDADSHEET INFORMATION
Colors=modulesPRE
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="PreDeploy_Network_WGCNA_OUTPUT.csv")

























#Next we want to Visualize how preserved Disease modules are in Healthy accordingly to both time-points.
pdf("PreDeploy_Colored_Modules.pdf",height=8,width=12)
par(mfrow=c(2,2))
plotDendroAndColors(geneTreeControlPre, colors=modulesPRE, main="Control Pre-Deployment Coloured with Control Pre-Deploy Modules", groupLabels=c("Module Colors"),dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
dev.off()

#Assess Module Preservation to quanitfy our visualization results
multiExpr  = list(CASEPD=list(data=(dataExpr1)),CONTROLPD=list(data=(dataExpr0)))
multiColor = list(CASEPD = modulesPRE)
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed", nPermutations=30,maxGoldModuleSize=100,maxModuleSize=1000)
stats = mp$preservation$Z$ref.CASEPD$inColumnsAlsoPresentIn.CONTROLPD
stats[order(-stats[,2]),c(1:2)]

modColors=rownames(stats)
moduleSize = stats$moduleSize
selectModules = !(modColors %in% c("grey", "gold"))
point.label = modColors[selectModules]

#Composite preservation statistics
Zsummary=stats$Zsummary.pres

# plot Zsummary versus module size
pdf("Zsummary_ControlPre.pdf")
plot(moduleSize[selectModules],Zsummary[selectModules], col = 1,
bg=modColors[selectModules],pch = 21,main="Zsummary Preservation",
cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
labelPoints(moduleSize[selectModules],Zsummary[selectModules],point.label,cex=1,offs=0.03)
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 2)
dev.off()


#Find kME values for Cases Pre-Deployment
geneModuleMembership1 = signedKME((dataExpr1), ME_PD)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(dataExpr1)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");
Gene= rownames(t(dataExpr1))
kMEtable1 = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
             kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))
write.csv(kMEtable1,"Control_Pre_kMEtable_for_ControlModules.csv",row.names=FALSE)

#Find kME Values for Controls Pre-Deployment
# First calculate MEs for A2, since we haven't done that yet
PCsAD = moduleEigengenes((dataExpr0),  colors=modulesPRE) 
MEsAD = PCsAD$eigengenes
geneModuleMembership2 = signedKME((dataExpr0), MEsAD)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 
MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(dataExpr0)[[2]]); 
colnames(MMPvalue2)=paste("PC",colorsPD,".pval",sep="");
kMEtable2  = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
 kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i])
colnames(kMEtable2)=colnames(kMEtable1)
write.csv(kMEtable2,"CasePre_kMEtable_for_ControlModules.csv",row.names=FALSE)


pdf("ControlPre_vs_CasePre_kMEtable_for_ControlModules.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsPD[c],
                    xlab="kME in Control Pre-Deploy",ylab="kME in Case Pre-Deploy")
}; dev.off()

#ASSESS HUB GENE CONSERVATION
pdf("inModule_ControlPre_vs_CasePre_kMEtable_for_ControlModules.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 inMod = modulesPRE== colorsPD[c]
 verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colorsPD[c],
                    xlab="kME in Control Pre-Deploy",ylab="Case Pre-Deploy")
}; dev.off()


#DETERMINE WHICH GENES ARE HUB GENES IN BOTH NETWORKS (GENES WITH HIGH KME VALUES IN BOTH NETWORKS), TOP 10 GENES PER MODULE BASED ON KME IN BOTH NETWORKS
topGenesKME = NULL
for (c in 1:length(colorsPD)){
 kMErank1    = rank(-geneModuleMembership1[,c])
 kMErank2    = rank(-geneModuleMembership2[,c])
 maxKMErank  = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
 topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; colnames(topGenesKME) = colorsPD
#PRINT TOP 10 GENES PER MODULE BASED ON KME IN BOTH NETWORKS
topGenesKME


# PLOT MODULE-TRAIT RELATIONSHIP
MEs0  =  moduleEigengenes(dataExpr0,  modulesPRE)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf("CasePre_Correlations.pdf")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.3,
zlim = c(-1,1),
main = paste("Case Pre-Deployment \n Module Eigengene Expression-Trait Relationships"))
dev.off()

#GENERATE One Network OUTPUT
Pre<-as.data.frame(datTraits1)
names(Pre)<-"Pre"
modNames = substring(names(MEs), 3)
nGenes = ncol(dataExpr1) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
nSamples = nrow(dataExpr1) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
geneModuleMembership<-as.data.frame(cor(dataExpr1,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr1, Pre, use = "p")) #  Generate correlations and p-value for each gene against the trait.  
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Pre), sep="")
names(GSPvalue) = paste("p.GS.", names(Pre), sep="")

#OUTPUT SPDADSHEET INFORMATION
Colors=modulesPRE
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="ControlPre_Network_WGCNA_OUTPUT.csv")

#Export to Cytoscape
modules = c("grey60","purple")
# Select module probes
inModule=is.finite(match(modulesPRE,modules))
Genes = row.names(t(dataExpr0))
modGenes=Genes[inModule]
match1=match(modGenes,GeneAnnotation$substanceBXH)
modGenes=GeneAnnotation$gene_symbol[match1]

modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files for Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
edgeFile=paste("CytoEdge",paste(modules,collapse="-"),".txt",sep=""),
nodeFile=paste("CytoNode",paste(modules,collapse="-"),".txt",sep=""),
weighted = TRUE, threshold = 0.02,nodeNames=modProbes,
altNodeNames = modGenes, nodeAttr = moduleColorsFemale[inModule])





























##################################################################
##################################################################
##################################################################
##############################  4   ##############################
##################################################################
##################################################################
##################################################################


###########################################################################################################
########## HERE WE CREATE FOUR LARGE NETWORK THEN REGRESS MODULES ONTO CONDITION AND TIME STATUS ##########
###########################################################################################################
library(WGCNA)
library(cluster)
options(stringsAsFactors  =  FALSE)
allowWGCNAThreads()

#INPUT CASE PD
PD <- read.delim("Case_Pre.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(PD) =names(PD)
dataExpr0<-as.data.frame(t(PD))
head(dataExpr0)
dim(dataExpr0)
gsg=goodSamplesGenes(dataExpr0,verbose=3)
gsg$allOK

nGenes = ncol(dataExpr0)
nSamples = nrow(dataExpr0)

traitData <- read.delim("Clinical_Traits_CasePre.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData)

#MATCH TRAITS TO PD-DEPLOYMENT
rowsExpr <- rownames(dataExpr0)
traitRows <- match(rowsExpr,traitData$studyid.1)
datTraits = traitData[traitRows, -1];
rownames (datTraits) = traitData[traitRows, 1];
table(rownames(datTraits)==rownames(dataExpr0)) # EVERYTHING OK?
names(datTraits)

powers=c(1:15) # in practice this should include powers up to 20.
sft0=pickSoftThreshold(dataExpr0,powerVector=powers, networkType="signed")

pdf("CasePre_Soft_Threshold.pdf")
par(mfrow=c(1,2))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],labels=powers,col="red")
abline(h=0.90,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

#INPUT CASE POST
CasePost <- read.delim("Case_Post.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(CasePost) =names(CasePost)
dataExpr1<-as.data.frame(t(CasePost))
head(dataExpr1)
dim(dataExpr1)
gsg=goodSamplesGenes(dataExpr1,verbose=3)
gsg$allOK

nGenes1 = ncol(dataExpr1)
nSamples1 = nrow(dataExpr1)

traitData1 <- read.delim("Clinical_Traits_CasePost.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData1)

#MATCH TRAITS TO PD-DEPLOYMENT
rowsExpr1 <- rownames(dataExpr1)
traitRows1 <- match(rowsExpr1,traitData1$studyid.1)
datTraits1 = traitData1[traitRows1, -1];
rownames (datTraits1) = traitData1[traitRows1, 1];
table(rownames(datTraits1)==rownames(dataExpr1)) # EVERYTHING OK?
names(datTraits1)

powers=c(1:15) # in practice this should include powers up to 20.
sft0=pickSoftThreshold(dataExpr1,powerVector=powers, networkType="signed")

pdf("CasePost_Soft_Threshold.pdf")
par(mfrow=c(1,2))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],labels=powers,col="red")
abline(h=0.90,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

#INPUT CONTROL PD
ControlPre <- read.delim("Control_Pre.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(ControlPre) =names(ControlPre)
dataExpr2<-as.data.frame(t(ControlPre))
head(dataExpr2)
dim(dataExpr2)
gsg=goodSamplesGenes(dataExpr2,verbose=3)
gsg$allOK

nGenes2 = ncol(dataExpr2)
nSamples2 = nrow(dataExpr2)

traitData2 <- read.delim("Clinical_Traits_ControlPre.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData2)

#MATCH TRAITS TO PD-DEPLOYMENT
rowsExpr2 <- rownames(dataExpr2)
traitRows2 <- match(rowsExpr2,traitData2$studyid.1)
datTraits2 = traitData2[traitRows2, -1];
rownames (datTraits2) = traitData2[traitRows2, 1];
table(rownames(datTraits2)==rownames(dataExpr2)) # EVERYTHING OK?
names(datTraits2)

powers=c(1:15) # in practice this should include powers up to 20.
sft0=pickSoftThreshold(dataExpr2,powerVector=powers, networkType="signed")

pdf("ControlPre_Soft_Threshold.pdf")
par(mfrow=c(1,2))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],labels=powers,col="red")
abline(h=0.90,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

#INPUT CONTROL POST
ControlPost <- read.delim("Control_Post.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(ControlPost) =names(ControlPost)
dataExpr3<-as.data.frame(t(ControlPost))
head(dataExpr3)
dim(dataExpr3)
gsg=goodSamplesGenes(dataExpr3,verbose=3)
gsg$allOK

nGenes3 = ncol(dataExpr3)
nSamples3 = nrow(dataExpr3)

traitData3 <- read.delim("Clinical_Traits_ControlPost.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData3)

#Match Traits to Control POST
rowsExpr3<- rownames(dataExpr3)
traitRows3 <- match(rowsExpr3,traitData3$studyid.1)
datTraits3 = traitData3[traitRows3, -1];
rownames (datTraits3) = traitData3[traitRows3, 1];
table(rownames(datTraits3)==rownames(dataExpr3)) # EVERYTHING OK?
names(datTraits3)

powers=c(1:15) # in practice this should include powers up to 20.
sft0=pickSoftThreshold(dataExpr3,powerVector=powers, networkType="signed")

pdf("ControlPost_Soft_Threshold.pdf")
par(mfrow=c(1,2))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],labels=powers,col="red")
abline(h=0.90,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

adjacencyCasePre = adjacency((dataExpr0),power=10 ,type="signed"); #FOR ONE LARGE NETWORK WE CHOOSE A POWER OF 10
diag(adjacencyCasePre)=0
dissTOMCasePre   = 1-TOMsimilarity(adjacencyCasePre, TOMType="signed")
geneTreeCasePre  = flashClust(as.dist(dissTOMCasePre), method="average")

adjacencyCasePost = adjacency((dataExpr1),power=10 ,type="signed"); #FOR ONE LARGE NETWORK WE CHOOSE A POWER OF 10
diag(adjacencyCasePost)=0
dissTOMCasePost   = 1-TOMsimilarity(adjacencyCasePost, TOMType="signed")
geneTreeCasePost  = flashClust(as.dist(dissTOMCasePost), method="average")

adjacencyControlPre = adjacency((dataExpr2),power=10 ,type="signed"); #FOR ONE LARGE NETWORK WE CHOOSE A POWER OF 10
diag(adjacencyControlPre)=0
dissTOMControlPre   = 1-TOMsimilarity(adjacencyControlPre, TOMType="signed")
geneTreeControlPre  = flashClust(as.dist(dissTOMControlPre), method="average")

adjacencyControlPost = adjacency((dataExpr3),power=10 ,type="signed"); #FOR ONE LARGE NETWORK WE CHOOSE A POWER OF 10
diag(adjacencyControlPost)=0
dissTOMControlPost   = 1-TOMsimilarity(adjacencyControlPost, TOMType="signed")
geneTreeControlPost  = flashClust(as.dist(dissTOMControlPost), method="average")

pdf("Four_Networks_GeneTree.pdf",height=6,width=12)
par(mfrow=c(2,2))
plot(geneTreeCasePre,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (Case Pre-Deployment)", labels=FALSE,hang=0.04);
plot(geneTreeCasePost,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (Case Post-Deployment)", labels=FALSE,hang=0.04);
plot(geneTreeControlPre,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (Control Pre-Deployment)", labels=FALSE,hang=0.04);
plot(geneTreeControlPost,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (Control Post-Deployment)", labels=FALSE,hang=0.04);
dev.off()

#LETS FIND DISEASE RELATED MODULES AND DETERMINE IF THEY ARE PDSERVED IN HEALTHY
#WE USE BOTH CASE AT PD AND POST DEPLOYMENT TO COMPARE TO CONTROLS

mColorh=NULL
for (ds in 0:3){
 tree = cutreeHybrid(dendro = geneTreeCasePre, pamStage=FALSE,
   minClusterSize = (15-3*ds), cutHeight = 0.99, 
   deepSplit = ds, distM = dissTOMCasePre)
 mColorh=cbind(mColorh,labels2colors(tree$labels));
}

mColorh2=NULL
for (ds in 0:3){
 tree2 = cutreeHybrid(dendro = geneTreeCasePost, pamStage=FALSE,
   minClusterSize = (15-3*ds), cutHeight = 0.99, 
   deepSplit = ds, distM = dissTOMCasePost)
 mColorh2=cbind(mColorh2,labels2colors(tree2$labels));
}

pdf("Two_Case_Networks_DeepSplit_Choices_30.pdf", height=10,width=25); 
plotDendroAndColors(geneTreeCasePre, mColorh, paste("dpSplt =",0:3), main = "Case Pre-Deployment Co-Expression Network",dendroLabels=FALSE);
plotDendroAndColors(geneTreeCasePost, mColorh2, paste("dpSplt =",0:3), main = "Case Post-Deployment Co-Expression Network",dendroLabels=FALSE);
dev.off()

#SET DEEP SPLIT CHOICE AND NAME OUR COLORS
modulesPRE =  mColorh[,2]
table(modulesPRE)

modulesPOST =  mColorh2[,2]
table(modulesPOST)

#Check to see if disease network modules can be cut and merged...

MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = flashClust(as.dist(MEDiss), method ="average")

MEList1 = moduleEigengenes(dataExpr1, colors=modulesPOST)
MEs1=MEList1$eigengenes
MEDiss1 = 1-cor(MEs1)
METree1 = flashClust(as.dist(MEDiss1), method ="average")

pdf("Two_Case_Networks_Module_Relationships.pdf")
plot(METree, main ="Clustering of Case Pre-Deployment Module Eigengenes",xlab ="",sub="")
MEDissThres = 0.1
abline(h=MEDissThres,col="red")
plot(METree1, main ="Clustering of Case Post-Deployment Module Eigengenes",xlab ="",sub="")
MEDissThres1 = 0.1
abline(h=MEDissThres1,col="red")
dev.off()

####################################### SHOULD WE MERGE MODULES BASED ON THE ABOVE? #######################################
merge = mergeCloseModules(dataExpr0, modulesPRE, cutHeight = MEDissThres, verbose=3)
mergedColors = merge$colors
table(mergedColors)
table(modulesPRE)
mergedMEs = merge$newME

merge1 = mergeCloseModules(dataExpr1, modulesPOST, cutHeight = MEDissThres1, verbose=3)
mergedColors1 = merge1$colors
table(mergedColors1)
table(modulesPOST)
mergedMEs1 = merge1$newME

#COMPARE UNMERGED MODULES TO MERGED MODULES
pdf("Two_Case_Network_Unmerged_Merged.pdf", w=9)
plotDendroAndColors(geneTreeCasePre, cbind(modulesPRE, mergedColors),c("Dynamic Tree Cut","Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(geneTreeCasePost, cbind(modulesPOST, mergedColors1),c("Dynamic Tree Cut","Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

modulesPRE = merge$colors
table(modulesPRE)
MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes

modulesPOST = merge1$colors
table(modulesPOST)
MEList1 = moduleEigengenes(dataExpr1, colors=modulesPOST)
MEs1=MEList1$eigengenes
##########################################################################################################################
#CALCULATE PC FOR VISUALIZATION FOR CASE PD-DEPLOYMENT
PCsPD    = moduleEigengenes((dataExpr0),  colors=modulesPRE) 
ME_PD    = PCsPD$eigengenes
distPCPD = 1-abs(cor(ME_PD,use="p"))
distPCPD = ifelse(is.na(distPCPD), 0, distPCPD)
pcTreePD = hclust(as.dist(distPCPD),method="average") 
MDS_PD   = cmdscale(as.dist(distPCPD),2)
colorsPD = names(table(modulesPRE))
names= row.names((dataExpr0))

pdf("CasePre_Module_Visualization.pdf",height=8,width=8)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 3) + 0.1, cex=1)
plot(pcTreePD, xlab="",ylab="",main="",sub="")
plot(MDS_PD, col= colorsPD,  main="MDS plot", cex=2, pch=19)

ordergenes = geneTreePD$order
plot.mat(scale(log(dataExpr0[ordergenes,])) ,
rlabels= modulesPRE[ordergenes], clabels= colnames(dataExpr0), rcols=modulesPRE[ordergenes])

for (which.module in names(table(modulesPRE)))
{ 
  ME = ME_PD[, paste("ME",which.module, sep="")] 
  n<- barplot(ME, col=which.module, main=paste(which.module, "module"), cex.main=1.5, ylab="eigengene expression",xlab="Case Pre-Deployment Samples")
  axis(1,at=n, labels=row.names(dataExpr0), las=2, cex.axis=0.5, font=2)
    #axis(1,at=n[c(0, 48,96,144)], las=1, cex.axis=0.005, font=2, tck=1, col.ticks = 'darkred')
    #axis(1,at=n[c(1, 48)], las=1, font=2, tck=1)
};
dev.off();

PCsPOST    = moduleEigengenes((dataExpr1),  colors=modulesPOST) 
ME_POST    = PCsPOST$eigengenes
distPCPOST = 1-abs(cor(ME_POST,use="p"))
distPCPOST = ifelse(is.na(distPCPOST), 0, distPCPOST)
pcTreePOST = hclust(as.dist(distPCPOST),method="average") 
MDS_POST   = cmdscale(as.dist(distPCPOST),2)
colorsPOST = names(table(modulesPOST))
names= row.names((dataExpr1))

pdf("CasePost_Module_Visualization.pdf",height=8,width=8)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 3) + 0.1, cex=1)
plot(pcTreePOST, xlab="",ylab="",main="",sub="")
plot(MDS_POST, col= colorsPOST,  main="MDS plot", cex=2, pch=19)
ordergenes = geneTreePOST$order
plot.mat(scale(log(dataExpr1[ordergenes,])) ,
rlabels= modulesPOST[ordergenes], clabels= colnames(dataExpr1), rcols=modulesPOST[ordergenes])
for (which.module in names(table(modulesPOST)))
{ 
  ME = ME_POST[, paste("ME",which.module, sep="")] 
  n<- barplot(ME, col=which.module, main=paste(which.module, "module"), cex.main=1.5, ylab="eigengene expression",xlab="Case Pre-Deployment Samples")
  axis(1,at=n, labels=row.names(dataExpr1), las=2, cex.axis=0.5, font=2)
    #axis(1,at=n[c(0, 48,96,144)], las=1, cex.axis=0.005, font=2, tck=1, col.ticks = 'darkred')
    #axis(1,at=n[c(1, 48)], las=1, font=2, tck=1)
};
dev.off();

#Next we want to Visualize how preserved Disease modules are in Healthy accordingly to both time-points.
pdf("Four_Networks_Colored_Modules.pdf",height=8,width=12)
par(mfrow=c(2,2))
plotDendroAndColors(geneTreeCasePre, colors=modulesPRE, main="Case Pre-Deployment Coloured with Case Pre-Deploy Modules",groupLabels=c("Module Colors"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
plotDendroAndColors(geneTreeCasePost, colors=modulesPOST, main="Case Post-Deployment Coloured with Case Post-Deploy Modules",groupLabels=c("Module Colors"),dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
plotDendroAndColors(geneTreeControlPre, colors=modulesPRE, main="Control Pre-Deployment Coloured with Case Pre-Deploy Modules", groupLabels=c("Module Colors"),dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
plotDendroAndColors(geneTreeControlPost, colors=modulesPOST, main="Control Post-Deployment Coloured with Case Post-Deploy Modules", groupLabels=c("Module Colors"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
dev.off()

#Assess Module Preservation to quanitfy our visualization results
multiExpr  = list(CASEPD=list(data=(dataExpr0)),CONTROLPD=list(data=(dataExpr2)))
multiColor = list(CASEPD = modulesPRE)
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed", nPermutations=30,maxGoldModuleSize=100,maxModuleSize=1000)
stats = mp$preservation$Z$ref.CASEPD$inColumnsAlsoPresentIn.CONTROLPD
stats[order(-stats[,2]),c(1:2)]

multiExpr2  = list(Case_Post=list(data=(dataExpr1)),Control_Post=list(data=(dataExpr3)))
multiColor2 = list(Case_Post = modulesPOST)
mp2=modulePreservation(multiExpr2,multiColor2,referenceNetworks=1,verbose=3,networkType="signed", nPermutations=30,maxGoldModuleSize=100,maxModuleSize=1000)
stats2 = mp2$preservation$Z$ref.Case_Post$inColumnsAlsoPresentIn.Control_Post
stats2[order(-stats2[,2]),c(1:2)]

#Find kME values for Cases Pre-Deployment
geneModuleMembership1 = signedKME((dataExpr0), ME_PD)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(dataExpr0)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");
Gene= rownames(t(dataExpr0))
kMEtable1 = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
             kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))
write.csv(kMEtable1,"Case_Pre_kMEtable.csv",row.names=FALSE)

#Find kME Values for Controls Pre-Deployment
# First calculate MEs for A2, since we haven't done that yet
PCsAD = moduleEigengenes((dataExpr2),  colors=modulesPRE) 
MEsAD = PCsAD$eigengenes
geneModuleMembership2 = signedKME((dataExpr2), MEsAD)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 
MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(dataExpr2)[[2]]); 
colnames(MMPvalue2)=paste("PC",colorsPD,".pval",sep="");
kMEtable2  = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
 kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i])
colnames(kMEtable2)=colnames(kMEtable1)
write.csv(kMEtable2,"ControlPre_kMEtable.csv",row.names=FALSE)


pdf("CasePre_vs_ControlPre_kMEtable.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsPD[c],
                    xlab="kME in Control Pre-Deploy",ylab="kME in Case Pre-Deploy")
}; dev.off()

#ASSESS HUB GENE CONSERVATION
pdf("inModule_CasePre_Vs_ControlPre_kMEtable.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 inMod = modulesPRE== colorsPD[c]
 verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colorsPD[c],
                    xlab="kME in Control Pre-Deploy",ylab="Case Pre-Deploy")
}; dev.off()


#DETERMINE WHICH GENES ARE HUB GENES IN BOTH NETWORKS (GENES WITH HIGH KME VALUES IN BOTH NETWORKS), TOP 10 GENES PER MODULE BASED ON KME IN BOTH NETWORKS
topGenesKME = NULL
for (c in 1:length(colorsPD)){
 kMErank1    = rank(-geneModuleMembership1[,c])
 kMErank2    = rank(-geneModuleMembership2[,c])
 maxKMErank  = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
 topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; colnames(topGenesKME) = colorsPD
#PRINT TOP 10 GENES PER MODULE BASED ON KME IN BOTH NETWORKS
topGenesKME

#Find kME Values for Cases Post-Deployment
geneModuleMembership3 = signedKME((dataExpr1), MEs1)
colnames(geneModuleMembership3)=paste("PC",colorsPOST,".cor",sep=""); 
MMPvalue3=corPvalueStudent(as.matrix(geneModuleMembership3),dim(dataExpr1)[[2]]); 
colnames(MMPvalue3)=paste("PC",colorsPOST,".pval",sep="");
Gene= rownames(t(dataExpr1))

kMEtable3 = cbind(Gene,Gene,modulesPOST)
for (i in 1:length(colorsPOST))
             kMEtable3 = cbind(kMEtable3, geneModuleMembership3[,i], MMPvalue3[,i])
colnames(kMEtable3)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership3), colnames(MMPvalue3))))
write.csv(kMEtable3,"CasePost_kMEtable.csv",row.names=FALSE)


#GRAB kME Values for Controls Pre-Deployment
# First calculate MEs for A2, since we haven't done that yet
PCs3A = moduleEigengenes((dataExpr3),  colors=modulesPRE) 
MEsCPost = PCs3A$eigengenes
geneModuleMembership4 = signedKME((dataExpr3), MEsCPost)
colnames(geneModuleMembership3)=paste("PC",colorsPOST,".cor",sep=""); 
MMPvalue4=corPvalueStudent(as.matrix(geneModuleMembership4),dim(dataExpr3)[[2]]); 
colnames(MMPvalue4)=paste("PC",colorsPOST,".pval",sep="");
kMEtable4  = cbind(Gene,Gene,modulesPOST)
for (i in 1:length(colorsPOST))
 kMEtable4 = cbind(kMEtable4, geneModuleMembership4[,i], MMPvalue4[,i])
colnames(kMEtable4)=colnames(kMEtable1)
write.csv(kMEtable4,"ControlPost_kMEtable.csv",row.names=FALSE)


#PLOT KME VALUES OF EACH GENE IN PD AGAINST THE CORRESPONDING KME VALUES OF EACH GENE IN POST.
pdf("CasePre_vs_ControlPre_kMEtable.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsPD[c],
                    xlab="kME in POST",ylab="kME in PD")
}; dev.off()

pdf("CasePost_vs_ControlPost_kMEtable.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsPD[c],
                    xlab="kME in POST",ylab="kME in PD")
}; dev.off()

#ASSESS HUB GENE CONSERVATION
pdf("inModule_POST_kMEtable2_vs_PD_kMEtable.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 inMod = modulesPRE== colorsPD[c]
 verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colorsPD[c],
                    xlab="kME in POST",ylab="kME in PD")
}; dev.off()


# PLOT MODULE-TRAIT RELATIONSHIP
MEs0  =  moduleEigengenes(dataExpr0,  modulesPRE)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf("CasePre_Correlations.pdf")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.3,
zlim = c(-1,1),
main = paste("Pre- and Post-Deployment \n Module Eigengene Expression-Trait Relationships"))
dev.off()

#GENERATE One Network OUTPUT
Pre<-as.data.frame(datTraits)
names(Pre)<-"Pre"
modNames = substring(names(MEs), 3)
nGenes = ncol(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
nSamples = nrow(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
geneModuleMembership<-as.data.frame(cor(dataExpr0,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr0, Pre, use = "p")) #  Generate correlations and p-value for each gene against the trait.  
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Pre), sep="")
names(GSPvalue) = paste("p.GS.", names(Pre), sep="")

#OUTPUT SPDADSHEET INFORMATION
Colors=modulesPRE
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="CasePre_Network_WGCNA_OUTPUT.csv")







datKME=signedKME(dataExpr0, MEs)
colorOfColumn=substring(names(datKME),4)

pdf("One_Network_Regress_Modules.pdf", h=10, w=14)
par(mfrow = c(3,4))
par(mfrow=c(4,length(selectModules)/4))
 for (module in names(table(modulesPRE))) 
  {
  column = match(module,colorOfColumn) 
  restModule=modulesPRE==module
  verboseScatterplot(datKME[restModule,column],GS.PTSD_V0[restModule],
  xlab=paste("Module Membership ",module,"module"),ylab="GS.PTSD_V0",
  main=paste("kME",module,"vs GS"),col=module)
  
  verboseScatterplot(datKME[restModule,column],GS.HC_V0[restModule],
  xlab=paste("Module Membership ",module,"module"),ylab="GS.HC_V0",
  main=paste("kME",module,"vs GS"),col=module)

  verboseScatterplot(datKME[restModule,column],GS.PTSD_V2[restModule],
  xlab=paste("Module Membership ",module,"module"),ylab="GS.PTSD_V2",
  main=paste("kME",module,"vs GS"),col=module)

  verboseScatterplot(datKME[restModule,column],GS.HC_V2[restModule],
  xlab=paste("Module Membership ",module,"module"),ylab="GS.HC_V2",
  main=paste("kME",module,"vs GS"),col=module)
  }
  dev.off()







#MODULE MEMBERSHIP (kME) TO COMPARE NETWORKS
#USED TO MEASURE CORRELATIONS BETWEEN EACH GENE AND EACH MODULE EIGENGENE
#PD
geneModuleMembership1 = signedKME((dataExpr0), ME_PD)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 

MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(dataExpr0)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");

Gene       = rownames(t(dataExpr0))
kMEtable1  = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
 kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))

write.csv(kMEtable1,"PD_kMEtable.csv",row.names=FALSE)



#PLOT KME VALUES OF EACH GENE IN PD AGAINST THE CORRESPONDING KME VALUES OF EACH GENE IN POST.
#MODULES WITH POINTS SHOWING A HIGH CORRELATION ARE HIGHLY PDSERVED.

pdf("PD_kMEtable_vs_POST_kMEtable.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsPD[c],
                    xlab="kME in POST",ylab="kME in PD")
}; dev.off()

#ASSESS HUB GENE CONSERVATION
pdf("inModule_POST_kMEtable2_vs_PD_kMEtable.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 inMod = modulesPRE== colorsPD[c]
 verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colorsPD[c],
                    xlab="kME in POST",ylab="kME in PD")
}; dev.off()








#Here we have settled with a Deep Split of 1, min modules 15, and merging modules at 0.1 height

#PLOT RELATIONS AMONG EIGENGENES AND THE TRAITS OF INTEREST
MET=orderMEs(cbind(MEs,PTSD_V0,HC_V0,PTSD_V2,HC_V2))
pdf("One_Network_Modules_Trait.pdf", h=16, w=15)
plotEigengeneNetworks(MET,"",marDendro=c(1,4,1,2), marHeatmap=c(6,6,4,4),cex.lab=0.8,xLabelsAngle=90)
dev.off()

# PLOT MODULE-TRAIT RELATIONSHIP
MEs0  =  moduleEigengenes(dataExpr0,  modulesPRE)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf("One_Network_Correlations.pdf")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.3,
zlim = c(-1,1),
main = paste("Pre- and Post-Deployment \n Module Eigengene Expression-Trait Relationships"))
dev.off()

#GENERATE One Network OUTPUT
Pre<-as.data.frame(datTraits)
names(Pre)<-"Pre"
modNames = substring(names(MEs), 3)
nGenes = ncol(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
nSamples = nrow(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
geneModuleMembership<-as.data.frame(cor(dataExpr0,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr0, Pre, use = "p")) #  Generate correlations and p-value for each gene against the trait.  
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Pre), sep="")
names(GSPvalue) = paste("p.GS.", names(Pre), sep="")

#OUTPUT SPDADSHEET INFORMATION
Colors=modulesPRE
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="One_Network_WGCNA_OUTPUT.csv")

#CALCULATE PC FOR VISUALIZATION FOR PD-DEPLOYMENT
PCsPD    = moduleEigengenes((dataExpr0),  colors=modulesPRE) 
ME_PD    = PCsPD$eigengenes
distPCPD = 1-abs(cor(ME_PD,use="p"))
distPCPD = ifelse(is.na(distPCPD), 0, distPCPD)
pcTreePD = hclust(as.dist(distPCPD),method="average") 
MDS_PD   = cmdscale(as.dist(distPCPD),2)
colorsPD = names(table(modulesPRE))
names= row.names((dataExpr0))

pdf("One_Network_Module_Visualization.pdf",height=8,width=8)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 3) + 0.1, cex=1)
plot(pcTreePD, xlab="",ylab="",main="",sub="")
plot(MDS_PD, col= colorsPD,  main="MDS plot", cex=2, pch=19)

ordergenes = geneTreePD$order
plot.mat(scale(log(dataExpr0[ordergenes,])) ,
rlabels= modulesPRE[ordergenes], clabels= colnames(dataExpr0), rcols=modulesPRE[ordergenes])

for (which.module in names(table(modulesPRE)))
{ 
  ME = ME_PD[, paste("ME",which.module, sep="")] 
  n<- barplot(ME, col=which.module, main=paste(which.module, "module"), cex.main=1.5, ylab="eigengene expression",xlab="      PTSD Pre-Deploy  PTSD Post-Deploy  H-Cntrl Pre-Deploy  H-Cntrl Post-Deploy")
  axis(1,at=n, labels=row.names(dataExpr0), las=2, cex.axis=0.3, font=2)
    axis(1,at=n[c(0, 48,96,144)], las=1, cex.axis=0.005, font=2, tck=1, col.ticks = 'darkred')
    #axis(1,at=n[c(1, 48)], las=1, font=2, tck=1)

};
dev.off();

# CALCULATE MODULE MEMBERSHIP VALUES (aka. module eigengene based connectivity kME)
datKME=signedKME(dataExpr0, MEs)
colorOfColumn=substring(names(datKME),4)

pdf("One_Network_Regress_Modules.pdf", h=10, w=14)
par(mfrow = c(3,4))
par(mfrow=c(4,length(selectModules)/4))
 for (module in names(table(modulesPRE))) 
  {
  column = match(module,colorOfColumn) 
  restModule=modulesPRE==module
  verboseScatterplot(datKME[restModule,column],GS.PTSD_V0[restModule],
  xlab=paste("Module Membership ",module,"module"),ylab="GS.PTSD_V0",
  main=paste("kME",module,"vs GS"),col=module)
  
  verboseScatterplot(datKME[restModule,column],GS.HC_V0[restModule],
  xlab=paste("Module Membership ",module,"module"),ylab="GS.HC_V0",
  main=paste("kME",module,"vs GS"),col=module)

  verboseScatterplot(datKME[restModule,column],GS.PTSD_V2[restModule],
  xlab=paste("Module Membership ",module,"module"),ylab="GS.PTSD_V2",
  main=paste("kME",module,"vs GS"),col=module)

  verboseScatterplot(datKME[restModule,column],GS.HC_V2[restModule],
  xlab=paste("Module Membership ",module,"module"),ylab="GS.HC_V2",
  main=paste("kME",module,"vs GS"),col=module)
  }
  dev.off()







#MODULE MEMBERSHIP (kME) TO COMPARE NETWORKS
#USED TO MEASURE CORRELATIONS BETWEEN EACH GENE AND EACH MODULE EIGENGENE
#PD
geneModuleMembership1 = signedKME((dataExpr0), ME_PD)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 

MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(dataExpr0)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");

Gene       = rownames(t(dataExpr0))
kMEtable1  = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
 kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))

write.csv(kMEtable1,"PD_kMEtable.csv",row.names=FALSE)



#PLOT KME VALUES OF EACH GENE IN PD AGAINST THE CORRESPONDING KME VALUES OF EACH GENE IN POST.
#MODULES WITH POINTS SHOWING A HIGH CORRELATION ARE HIGHLY PDSERVED.

pdf("PD_kMEtable_vs_POST_kMEtable.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsPD[c],
                    xlab="kME in POST",ylab="kME in PD")
}; dev.off()

#ASSESS HUB GENE CONSERVATION
pdf("inModule_POST_kMEtable2_vs_PD_kMEtable.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 inMod = modulesPRE== colorsPD[c]
 verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colorsPD[c],
                    xlab="kME in POST",ylab="kME in PD")
}; dev.off()





































#INPUT PD
PD <- read.delim("MRSII_PreDeploy.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(PD) =names(PD)
dataExpr0<-as.data.frame(t(PD))
head(dataExpr0)
dim(dataExpr0)
gsg=goodSamplesGenes(dataExpr0,verbose=3)
gsg$allOK

nGenes = ncol(dataExpr0)
nSamples = nrow(dataExpr0)

#INPUT POST
POST <- read.delim("MRSII_PostDeploy.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(POST) =names(POST)
dataExpr1<-as.data.frame(t(POST))
head(dataExpr1)
dim(dataExpr1)
gsg=goodSamplesGenes(dataExpr1,verbose=3)
gsg$allOK

nGenes1 = ncol(dataExpr1)
nSamples1 = nrow(dataExpr1)

#INPUT CLINICAL TRIATS
traitData <- read.delim("Clinical_Traits_MRSII.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData)

#MATCH TRAITS TO PD-DEPLOYMENT
rowsExpr <- rownames(dataExpr0)
traitRows <- match(rowsExpr,traitData$studyid.1)
datTraits = traitData[traitRows, -1];
rownames (datTraits) = traitData[traitRows, 1];
table(rownames(datTraits)==rownames(dataExpr0)) # EVERYTHING OK?
names(datTraits)

#MATCH TRAITS TO POST-DEPLOYMENT
rowsExpr1 <- rownames(dataExpr1)
traitRows1 <- match(rowsExpr1,traitData$studyid.1)
datTraits1 = traitData[traitRows1, -1];
rownames (datTraits1) = traitData[traitRows1, 1];
table(rownames(datTraits1)==rownames(dataExpr1)) # EVERYTHING OK?
names(datTraits1)


#CONFIGURE SOFT-THRESHOLDS
powers=c(1:15) # in practice this should include powers up to 20.
sft0=pickSoftThreshold(dataExpr0,powerVector=powers, networkType="signed")

powers=c(1:15) # in practice this should include powers up to 20.
sft1=pickSoftThreshold(dataExpr1,powerVector=powers, networkType="signed")


#VIEW SOFT-THRESHOLD RESULTS - NEGATIVE SLOPE
pdf("SOFT_THRESHOLD_PreDeploy.pdf")
par(mfrow=c(1,2))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],labels=powers,col="red")
abline(h=0.90,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

pdf("SOFT_THRESHOLD_PostDeploy.pdf")
par(mfrow=c(1,2))
plot(sft1$fitIndices[,1],-sign(sft1$fitIndices[,3])*sft1$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft1$fitIndices[,1],-sign(sft1$fitIndices[,3])*sft1$fitIndices[,2],labels=powers,col="red")
abline(h=0.90,col="red")#CHOOSE A  R^2 CUT-OFF OF H
plot(sft1$fitIndices[,1],sft1$fitIndices[,5],type="n",
xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft1$fitIndices[,1],sft1$fitIndices[,5],labels=powers,col="red")
dev.off()

#CREATE ADJACENCY MATRICES AND GENE-TREES
adjacencyPreDeploy = adjacency((dataExpr0),power=8 ,type="signed");
diag(adjacencyPreDeploy)=0
dissTOMPreDeploy   = 1-TOMsimilarity(adjacencyPreDeploy, TOMType="signed")
geneTreePreDeploy  = flashClust(as.dist(dissTOMPreDeploy), method="average")

adjacencyPostDeploy = adjacency((dataExpr1),power=8 ,type="signed");
diag(adjacencyPostDeploy)=0
dissTOMPostDeploy   = 1-TOMsimilarity(adjacencyPostDeploy, TOMType="signed")
geneTreePostDeploy  = flashClust(as.dist(dissTOMPostDeploy), method="average")

pdf("View_Resulting_Gene_Trees.pdf",height=6,width=12)
par(mfrow=c(2,2))
plot(geneTreePreDeploy,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (Pre-Deployment)", labels=FALSE,hang=0.04);
plot(geneTreePostDeploy,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (Post-Deployment)", labels=FALSE,hang=0.04); 
dev.off()

mColorh=NULL
for (ds in 0:3){
 tree = cutreeHybrid(dendro = geneTreePreDeploy, pamStage=FALSE,
   minClusterSize = (30-3*ds), cutHeight = 0.99, 
   deepSplit = ds, distM = dissTOMPreDeploy)
 mColorh=cbind(mColorh,labels2colors(tree$labels));
}

pdf("PreDeploy_DeepSplit_Choices.pdf", height=10,width=25); 
plotDendroAndColors(geneTreePreDeploy, mColorh, paste("dpSplt =",0:3), main = "Pre-Deploy",dendroLabels=FALSE);
dev.off()

#SET DEEP SPLIT CHOICE AND NAME OUR COLORS
modulesPRE =  mColorh[,3]
table(modulesPRE)

MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = flashClust(as.dist(MEDiss), method ="average")

pdf("PreDeploy_Module_Relationships.pdf")
plot(METree, main ="Clustering of module eigengenes",xlab ="",sub="")
MEDissThres = 0.13
abline(h=MEDissThres,col="red")
dev.off()

merge = mergeCloseModules(dataExpr0, modulesPRE, cutHeight = MEDissThres, verbose=3)
modulesPRE = merge$colors
mergedMEs = merge$newME

#COMPARE UNMERGED MODULES TO MERGED MODULES
pdf("PreDeploy_Compare_Unmerged_2_Merged_Modules.pdf", w=9)
plotDendroAndColors(geneTreePreDeploy, cbind(modulesPRE, mergedColors),c("Dynamic Tree Cut","Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes


######## NOW POST-DEPLOYMENT ############

mColorh2=NULL
for (ds in 0:3){
 tree2 = cutreeHybrid(dendro = geneTreePostDeploy, pamStage=FALSE,
   minClusterSize = (30-3*ds), cutHeight = 0.99, 
   deepSplit = ds, distM = dissTOMPostDeploy)
 mColorh2=cbind(mColorh2,labels2colors(tree2$labels));
}

pdf("PostDeploy_DeepSplit_Choices.pdf", height=10,width=25); 
plotDendroAndColors(geneTreePostDeploy, mColorh2, paste("dpSplt =",0:3), main = "Post-Deploy",dendroLabels=FALSE);
dev.off()

#SET DEEP SPLIT CHOICE
modulesPOST =  mColorh2[,3]
table(modulesPOST)

MEList1 = moduleEigengenes(dataExpr1, colors=modulesPOST)
MEs1=MEList1$eigengenes
MEDiss1 = 1-cor(MEs1)
METree1 = flashClust(as.dist(MEDiss1), method ="average")

pdf("PostDeploy_Module_Relationships.pdf")
plot(METree1, main ="Clustering of module eigengenes",xlab ="",sub="")
MEDissThres1 = 0.1
abline(h=MEDissThres1,col="red")
dev.off()

merge1 = mergeCloseModules(dataExpr1, modulesPOST, cutHeight = MEDissThres1, verbose=3)
mergedColors1 = merge1$colors
mergedMEs1 = merge1$newME

#COMPARE UNMERGED MODULES TO MERGED MODULES
pdf("PostDeploy_Compare_Unmerged_2_Merged_Modules.pdf", w=9)
plotDendroAndColors(geneTreePostDeploy, cbind(modulesPOST, mergedColors1),c("Dynamic Tree Cut","Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

MEList1 = moduleEigengenes(dataExpr0, colors=mnodulesPOST)
MEs1=MEList$eigengenes

###############################################################
#VIEW EACH NETWORK AND GS MEASURE FOR EACH VARIABLE OF INTEREST.
PTSD_V2 = as.data.frame(datTraits$V2_PTSD)
names(PTSD_V2)="PTSD_V2"
GS.PTSD_V2=as.numeric(cor(dataExpr0,PTSD_V2,use="p"))
GS.PTSD_V2Color=numbers2colors(GS.PTSD_V2,signed=T)

HC_V2 = as.data.frame(datTraits1$V2_HealthyCntl)
names(HC_V2)="HC_V2"
GS.HC_V2=as.numeric(cor(dataExpr1,HC_V2,use="p"))
GS.HC_V2Color=numbers2colors(GS.HC_V2,signed=T)

rowsExpr <- rownames(dataExpr0)
traitRows <- match(rowsExpr,traitData$studyid.1)
datTraits = traitData[traitRows, -1];


datColors0=data.frame(modulesPRE,GS.PTSD_V0Color,GS.HC_V0Color)
datColors1=data.frame(modulesPRE,GS.PTSD_V2Color,GS.HC_V2Color)
datColors2=data.frame(modulesPOST,GS.PTSD_V2Color,GS.HC_V2Color)
datColors3=data.frame(modulesPOST,GS.PTSD_V0Color,GS.HC_V0Color)

pdf("FINAL_MODULES.pdf",height=8,width=12)
plotDendroAndColors(geneTreePreDeploy, colors=datColors0, main="Pre-Deployment Coloured with Pre-Deploy Modules", groupLabels=c("Module colors","GS.PTSD","GS.Healthy Cntl"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
plotDendroAndColors(geneTreePostDeploy, colors=datColors1, main="Post-Deployment Coloured with Pre-Deploy Modules", groupLabels=c("Module colors","GS.PTSD","GS.Healthy Cntl"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
plotDendroAndColors(geneTreePostDeploy, colors=datColors2, main="Post-Deployment Coloured with Post-Deploy Modules", groupLabels=c("Module colors","GS.PTSD","GS.Healthy Cntl"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
plotDendroAndColors(geneTreePreDeploy, colors=datColors3, main="Pre-Deployment Coloured with Post-Deploy Modules", groupLabels=c("Module colors","GS.PTSD","GS.Healthy Cntl"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05) 
dev.off()

#Calculate Module Preservation Statistics
multiExpr  = list(PD=list(data=(dataExpr0)),POST=list(data=(dataExpr1)))
multiColor = list(PD = modulesPRE)
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
nPermutations=30,maxGoldModuleSize=100,maxModuleSize=1000)
stats = mp$preservation$Z$ref.PD$inColumnsAlsoPresentIn.POST
#THE HIGHER Zsummary.pres = The More Preserved the module.
# If Z > 5 and < 10 = Moderate Preservation.
# If Z > 10  = High Preservation.
stats[order(-stats[,2]),c(1:2)]

multiExpr2  = list(POST=list(data=(dataExpr1)),PD=list(data=(dataExpr0)))
multiColor2 = list(POST = modulesPOST)
mp2=modulePreservation(multiExpr2,multiColor2,referenceNetworks=1,verbose=3,networkType="signed",
nPermutations=30,maxGoldModuleSize=100,maxModuleSize=1000)
stats2 = mp2$preservation$Z$ref.POST$inColumnsAlsoPresentIn.PD
#THE HIGHER Zsummary.pres = The More Preserved the module.
# If Z > 5 and < 10 = Moderate Preservation.
# If Z > 10  = High Preservation.
stats2[order(-stats2[,2]),c(1:2)]



#PLOT RELATIONS AMONG EIGENGENES AND THE TRAITS OF INTEREST
MET=orderMEs(cbind(MEs,PTSD_V0,HC_V0))
pdf("PreDeploy_MODULE_RELATIONS_PLUS_TRAIT.pdf", h=16, w=15)
plotEigengeneNetworks(MET,"",marDendro=c(1,4,1,2), marHeatmap=c(6,6,4,4),cex.lab=0.8,xLabelsAngle=90)
dev.off()

MET1=orderMEs(cbind(MEs1,PTSD_V2,HC_V2))
pdf("PostDeploy_MODULE_RELATIONS_PLUS_TRAIT.pdf", h=16, w=15)
plotEigengeneNetworks(MET1,"",marDendro=c(1,4,1,2), marHeatmap=c(6,6,4,4),cex.lab=0.8,xLabelsAngle=90)
dev.off()


# PLOT MODULE-TRAIT RELATIONSHIP
MEs0  =  moduleEigengenes(dataExpr0,  modulesPRE)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf("PreDeploy_MODULETRAIT_RELATIONS.pdf")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.3,
zlim = c(-1,1),
main = paste("Pre-Deployment \n Module-Trait Relationships"))
dev.off()

MEs2  =  moduleEigengenes(dataExpr1,  modulesPOST)$eigengenes
MEs3 = orderMEs(MEs2)
moduleTraitCor = cor(MEs3, datTraits1, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf("PostDeploy_MODULETRAIT_RELATIONS.pdf")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs3),
ySymbols = names(MEs3),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.3,
zlim = c(-1,1),
main = paste("Post-Deployment \n Module-Trait Relationships"))
dev.off()

#GENERATE Pre-Deployment OUTPUT
Pre<-as.data.frame(datTraits)
names(Pre)<-"Pre"
modNames = substring(names(MEs), 3)
nGenes = ncol(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
nSamples = nrow(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
geneModuleMembership<-as.data.frame(cor(dataExpr0,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr0, Pre, use = "p")) #  Generate correlations and p-value for each gene against the trait.  
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Pre), sep="")
names(GSPvalue) = paste("p.GS.", names(Pre), sep="")

#OUTPUT SPDADSHEET INFORMATION
Colors=modulesPRE
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="PreDeploy_WGCNA_OUTPUT_Common.csv")


#GENERATE Post-Deployment OUTPUT
Post<-as.data.frame(datTraits1)
names(Post)<-"Post"
modNames = substring(names(MEs3), 3)
nGenes = ncol(dataExpr1) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
nSamples = nrow(dataExpr1) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
geneModuleMembership<-as.data.frame(cor(dataExpr1,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr1, Pre, use = "p")) #  Generate correlations and p-value for each gene against the trait.  
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Pre), sep="")
names(GSPvalue) = paste("p.GS.", names(Pre), sep="")

#OUTPUT SPDADSHEET INFORMATION
Colors=modulesPOST
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="PostDeploy_WGCNA_OUTPUT_Common.csv")






#CALCULATE PC FOR VISUALIZATION FOR PD-DEPLOYMENT
PCsPD    = moduleEigengenes((dataExpr0),  colors=modulesPRE) 
ME_PD    = PCsPD$eigengenes
distPCPD = 1-abs(cor(ME_PD,use="p"))
distPCPD = ifelse(is.na(distPCPD), 0, distPCPD)
pcTreePD = hclust(as.dist(distPCPD),method="average") 
MDS_PD   = cmdscale(as.dist(distPCPD),2)
colorsPD = names(table(modulesPRE))
names= row.names((dataExpr0))

pdf("MODULE_EIGENGENE_VISUALIZATION_PDDEPLOY.pdf",height=8,width=8)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 3) + 0.1, cex=1)
plot(pcTreePD, xlab="",ylab="",main="",sub="")
plot(MDS_PD, col= colorsPD,  main="MDS plot", cex=2, pch=19)

ordergenes = geneTreePD$order
plot.mat(scale(log(dataExpr0[ordergenes,])) ,
rlabels= modulesPRE[ordergenes], clabels= colnames(dataExpr0), rcols=modulesPRE[ordergenes])

for (which.module in names(table(modulesPRE))){
  ME = ME_PD[, paste("ME",which.module, sep="")] 
  barplot(ME, col=which.module, main="", cex.main=2, 
      ylab="eigengene expression",xlab="array sample") 
}; 

dev.off();



#MODULE MEMBERSHIP (kME) TO COMPARE NETWORKS
#USED TO MEASURE CORRELATIONS BETWEEN EACH GENE AND EACH MODULE EIGENGENE
#PD
geneModuleMembership1 = signedKME((dataExpr0), ME_PD)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 

MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(dataExpr0)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");

Gene       = rownames(t(dataExpr0))
kMEtable1  = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
 kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))

write.csv(kMEtable1,"PD_kMEtable.csv",row.names=FALSE)

#POST
# First calculate MEs for A2, since we haven't done that yet
PCsPOST = moduleEigengenes((dataExpr1),  colors=modulesPRE) 
ME_POST = PCsPOST$eigengenes

geneModuleMembership2 = signedKME((dataExpr1), ME_POST)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 

MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(dataExpr1)[[2]]); 
colnames(MMPvalue2)=paste("PC",colorsPD,".pval",sep="");

kMEtable2  = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
 kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i])
colnames(kMEtable2)=colnames(kMEtable1)

write.csv(kMEtable2,"POST_kMEtable.csv",row.names=FALSE)


#PLOT KME VALUES OF EACH GENE IN PD AGAINST THE CORRESPONDING KME VALUES OF EACH GENE IN POST.
#MODULES WITH POINTS SHOWING A HIGH CORRELATION ARE HIGHLY PDSERVED.

pdf("PD_kMEtable_vs_POST_kMEtable.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsPD[c],
                    xlab="kME in POST",ylab="kME in PD")
}; dev.off()

#ASSESS HUB GENE CONSERVATION
pdf("inModule_POST_kMEtable2_vs_PD_kMEtable.pdf",height=8,width=8)
for (c in 1:length(colorsPD)){
 inMod = modulesPRE== colorsPD[c]
 verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colorsPD[c],
                    xlab="kME in POST",ylab="kME in PD")
}; dev.off()


#DETERMINE WHICH GENES ARE HUB GENES IN BOTH NETWORKS (GENES WITH HIGH KME VALUES IN BOTH NETWORKS)
#TOP 10 GENES PER MODULE BASED ON KME IN BOTH NETWORKS
topGenesKME = NULL
for (c in 1:length(colorsPD)){
 kMErank1    = rank(-geneModuleMembership1[,c])
 kMErank2    = rank(-geneModuleMembership2[,c])
 maxKMErank  = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
 topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; colnames(topGenesKME) = colorsPD
#PRINT TOP 10 GENES PER MODULE BASED ON KME IN BOTH NETWORKS
topGenesKME


#EXPORT FOR VISUALIZATION
for (co in colorsPD[colorsPD!="grey"])
  visantPrepOverall(modulesPRE, co, t(dataExpr0), rownames(dataExpr0), 500, softPower, TRUE)

for (co in colorsPD[colorsPD!="grey"])
  visantPrepOverall(modulesA1, co, t(dataExpr1), rownames(dataExpr1), 500, softPower, TRUE)







# CALCULATE MODULE MEMBERSHIP VALUES (aka. module eigengene based connectivity kME)
datKME=signedKME(dataExpr0, MEs)
colorOfColumn=substring(names(datKME),4)
pdf("MODULE_CORRELATIONS_Common.pdf", h=10, w=12)
par(mfrow = c(4,4))
selectModules=c("lightgreen","tan","red","purple", "turquoise", "greenyellow", "brown", "midnightblue" , "pink", "cyan", "magenta", "yellow")
par(mfrow=c(4,length(selectModules)/4))
for (module in selectModules) {
column = match(module,colorOfColumn)
restModule=moduleColors==module
verboseScatterplot(datKME[restModule,column],GS.V0_EP[restModule],
xlab=paste("Module Membership ",module,"module"),ylab="GS.V0_EP",
main=paste("kME.",module,"vs. GS"),col=module)}
dev.off()


#HERE WE CREATE INFORMATION FOR CREATING CO-EXPDSSION NETWORKS FOR CYTOSCAPE UTILITY
#RECALCULATE OVERLAP
TOM = TOMsimilarityFromExpr(dataExpr0,power= 3, networkType="signed", corType="bicor") # This takes a while....
#SELECT YOUR MODULE OF CHOICE
modules=("yellow")
probes=names(dataExpr0)
inModule=is.finite(match(modulesPRE, modules));
modProbes=probes[inModule]
modTom=TOM[inModule, inModule]
dimnames(modTom)=list(modProbes, modProbes)
#EXPORT
cyt=exportNetworkToCytoscape(modTom,
edgeFile=paste("CytoscapeInput-edges-", paste(modules, collapose="-") , ".txt", sep=""),
nodeFile=paste("CytoscapeInput-nodes-", paste(modules, collapose="-") , ".txt", sep=""),
weighted=TRUE,
threshold=0.02,
nodeNames=modProbes,
#altNodeNAmes=modGenes,
nodeAttr=moduleColors[inModule])
write.table (cyt$edgeData, "CytoScape_edge_MRSII.txt", sep="\t")
write.table (cyt$nodeData, "CytoScape_node_MRSII.txt", sep="\t")





























Color <- read.delim("list2.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
head(Color)
modulesNoLit = Color$Colors
overlap1 = overlapTable(modulesNoLit, Color$Depression,ignore = "silver")
overlap2 = overlapTable(modulesNoLit, Color$PTSD,ignore = "silver")
overlap3 = overlapTable(modulesNoLit, Color$TE,ignore = "silver")
overlap4 = overlapTable(modulesNoLit, Color$PTSD.Dep,ignore = "silver")

overlap5 = overlapTable(modulesNoLit, Color$DeReubis,ignore = "silver")
overlap6 = overlapTable(modulesNoLit, Color$Sanders,ignore = "silver")
overlap7 = overlapTable(modulesNoLit, Color$ParkID,ignore = "silver")
overlap8 = overlapTable(modulesNoLit, Color$Fromer,ignore = "silver")

overlap9 = overlapTable(modulesNoLit, Color$Astrocytes,ignore = "silver")
overlap10 = overlapTable(modulesNoLit, Color$Neuron,ignore = "silver")
overlap11 = overlapTable(modulesNoLit, Color$Oligodendrocytes,ignore = "silver")


overlap1 = overlapTable(modulesNoLit, Color$Astrocytes,ignore = "silver")
overlap2 = overlapTable(modulesNoLit, Color$Neuron,ignore = "silver")
overlap3 = overlapTable(modulesNoLit, Color$Oligodendrocytes,ignore = "silver")
overlap4 = overlapTable(modulesNoLit, Color$Zeisel_CA1_pyrNeurons,ignore = "silver")
overlap5 = overlapTable(modulesNoLit, Color$Zeisel_S1_pyrNeurons,ignore = "silver")
overlap6 = overlapTable(modulesNoLit, Color$Zeisel_Interneurons,ignore = "silver")
overlap7 = overlapTable(modulesNoLit, Color$Zeisel_Oligodendrocytes,ignore = "silver")
overlap8 = overlapTable(modulesNoLit, Color$Zeisel_Astrocytes,ignore = "silver")
overlap9 = overlapTable(modulesNoLit, Color$Zeisel_Microglial,ignore = "silver")
overlap10 = overlapTable(modulesNoLit, Color$Zeisel_Endothelial,ignore = "silver")


overlap1 = overlapTable(modulesNoLit, Color$DeReubis,ignore = "silver")
overlap1 = overlapTable(modulesNoLit, Color$Sanders,ignore = "silver")
overlap2 = overlapTable(modulesNoLit, Color$ParkID,ignore = "silver")
overlap3 = overlapTable(modulesNoLit, Color$Fromer,ignore = "silver")
overlap4 = overlapTable(modulesNoLit, Color$DD,ignore = "silver")
overlap5 = overlapTable(modulesNoLit, Color$Epilepsy,ignore = "silver")


numMat1=-log10(overlap1$pTable)
numMat1[numMat1 >100] =100;
numMat2=-log10(overlap2$pTable)
numMat2[numMat2 >500] =500;
numMat3=-log10(overlap3$pTable)
numMat3[numMat3 >500] =500;
numMat4=-log10(overlap4$pTable)
numMat4[numMat4 >500] =500;
numMat5=-log10(overlap5$pTable)
numMat5[numMat5 >500] =500;
numMat6=-log10(overlap6$pTable)
numMat6[numMat6 >500] =500;
numMat7=-log10(overlap7$pTable)
numMat7[numMat7 >500] =500;
numMat8=-log10(overlap8$pTable)
numMat8[numMat8 >500] =500;
numMat9=-log10(overlap9$pTable)
numMat9[numMat9 >500] =500;
numMat10=-log10(overlap10$pTable)
numMat10[numMat10 >500] =500;
numMat11=-log10(overlap11$pTable)
numMat11[numMat11 >500] =500;

textMat1 =paste(overlap1$countTable, "\n", signif(overlap1$pTable, 2));
textMat2 =paste(overlap2$countTable, "\n", signif(overlap2$pTable, 2));
textMat3 =paste(overlap3$countTable, "\n", signif(overlap3$pTable, 2));
textMat4 =paste(overlap4$countTable, "\n", signif(overlap4$pTable, 2));
textMat5 =paste(overlap5$countTable, "\n", signif(overlap5$pTable, 2));
textMat6 =paste(overlap6$countTable, "\n", signif(overlap6$pTable, 2));
textMat7 =paste(overlap7$countTable, "\n", signif(overlap7$pTable, 2));
textMat8 =paste(overlap8$countTable, "\n", signif(overlap8$pTable, 2));
textMat9 =paste(overlap9$countTable, "\n", signif(overlap9$pTable, 2));
textMat10 =paste(overlap10$countTable, "\n", signif(overlap10$pTable, 2));
textMat11 =paste(overlap11$countTable, "\n", signif(overlap11$pTable, 2));


dim(textMat1) =dim(numMat1)
dim(textMat2) =dim(numMat2)
dim(textMat3) =dim(numMat3)
dim(textMat4) =dim(numMat4)
dim(textMat5) =dim(numMat5)
dim(textMat6) =dim(numMat6)
dim(textMat7) =dim(numMat7)
dim(textMat8) =dim(numMat8)
dim(textMat9) =dim(numMat9)
dim(textMat10) =dim(numMat10)
dim(textMat11) =dim(numMat11)

largetextMat = cbind(textMat1,textMat2,textMat3,textMat4,textMat5,textMat6,textMat7,textMat8,textMat9,textMat10,textMat11)
numMatlarge = cbind(numMat1,numMat2,numMat3,numMat4,numMat5,numMat6,numMat7,numMat8,numMat9,numMat10,numMat11)

yLabels = paste ("M", sort(unique(modulesNoLit)));
ySymbols = paste ("M", sort(unique(modulesNoLit)), ": ", table (modulesNoLit), sep="")

fcex=.2;
pcex=.5;
fcex1=.5;
pcexl=.5;
fp=FALSE

#xLabels =c("DE", "Treated", "Contrast", "Zhang et al 2005", "Seelan et al 2008", "McEachin et al 2010", "Lowthert et al 2012", "Beech et al 2014", "Wantanabe et al 2014", "Hunsberger et al 2015", "Gupta et al 2012", "DisGenNet 2015")

xLabels =c("Depression","PTSD",  "TE", "PTSD/Dep","DeRubeis","Sanders","ParkID","Fromer","Astocytes","Neuron","Oligodendrocytes")

xLabels =c("Astrocytes", "Neuron", "Oligodendrocytes", "Zeisel_CA1_pyrNeurons", "Zeisel_S1_pyrNeurons", "Zeisel_Interneurons", "Zeisel_Oligodendrocytes", "Zeisel_Astrocytes", "Zeisel_Microglial", "Zeisel_Endothelial")

ySymbols=c("M1: 333","M2: 1341","M3: 945","M4: 74",
 "M5: 596","M6: 3896","M7: 212","M8: 287",
"M9: 171","M10: 369","M11: 92","M12: 95",
"M13: 1359", "M14: 935")
 
par(mar=c(4,7,7,4));
colfunc <- colorRampPalette(c("white","white", "tomato"))
labeledHeatmap(Matrix=t(numMatlarge),
               xLabels=yLabels, xSymbols=ySymbols,
               yLabels=xLabels, #ySymbols=ySymbols,
               colorLabels=F, 
               colors=colfunc(50),
               textMatrix=t(largetextMat),cex.text=0.001, setStdMargins=FALSE,
               cex.lab=0.6,
               xColorWidth=1,
               main="", cex.main=1)

labeledHeatmap(Matrix=(numMatlarge),
               yLabels=yLabels, ySymbols=ySymbols,
               xLabels=xLabels, #ySymbols=ySymbols,
               colorLabels=F, 
               colors=colfunc(50),
               textMatrix=(largetextMat),cex.text=0.5, setStdMargins=FALSE,
               cex.lab=0.6,
               xColorWidth=1,
               main="", cex.main=1)


numMatlarge <- read.delim("list2.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
head(numMatlarge)
yLabels=c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14")
ySymbols=c("CD14+ Monocytes", "CD56+ NK Cells", "CD4+ T Cells", "CD8+ T Cells", "CD19+ B cells", "CD33+ Myeloid", "CD34+", "CD71+ Early Erythroid", "Astrocytes", "Neurons", "Oligodendrocytes", "Microglia", "Endothelial")
labeledHeatmap(Matrix=(numMatlarge),
               yLabels=ySymbols, ySymbols=ySymbols,
               xLabels=xLabels, #ySymbols=ySymbols,
               colorLabels=T, 
               colors=blueWhiteRed(100)[50:100],
               textMatrix=(numMatlarge),cex.text=0.001, setStdMargins=FALSE,
               cex.lab=0.5,
               xColorWidth=1,
               main="", cex.main=1)



matrix <- read.delim("list2.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
head(matrix)
par(mfrow=c(3,3))
barplot(-log10(matrix$Depression), col="white", cex.axis=0.7, cex=0.7,names.arg=rownames(matrix), las=2,main="Dep")
abline(h=-log10(0.05), col="red")
barplot(-log10(matrix$PTSD), col="white", cex.axis=0.7, cex=0.7,names.arg=rownames(matrix), las=2,main="PTSD")
abline(h=-log10(0.05), col="red")
barplot(-log10(matrix$Neuron), col="white", cex.axis=0.7, cex=0.7,names.arg=rownames(matrix), las=2,main="Neuron")
abline(h=-log10(0.05), col="red")
barplot(-log10(matrix$Sanders), col="white", cex.axis=0.7, cex=0.7,names.arg=rownames(matrix), las=2,main="Sanders")
abline(h=-log10(0.05), col="red")


par(mfrow=c(3,3))
barplot(-log10(matrix$Depression), col=matrix$Module, border=matrix$Module,cex.axis=0.8, cex=0.8,names.arg=rownames(matrix), las=2,main="Dep")
abline(h=-log10(0.05), col="black")
barplot(-log10(matrix$PTSD), col=matrix$Module, border=matrix$Module,cex.axis=0.8, cex=0.8,names.arg=rownames(matrix), las=2,main="PTSD",xaxt='n',yaxt='n')
par(new=TRUE)
barplot(-log10(matrix$PTSD_c), col="white",border=matrix$Module,cex.axis=0.8, cex=0.8,names.arg=rownames(matrix), las=2,main="PTSD", ylim=c(0,43))
abline(h=-log10(0.05), col="black")
par(new=TRUE)
barplot(-log10(matrix$PTSD), col="white", border=matrix$Module,cex.axis=0.8, cex=0.7,names.arg=rownames(matrix), las=2,main="PTSD")
abline(h=-log10(0.05), col="black")

colx=col=c("blue", "blue","blue", "blue","blue","turquoise", "turquoise","turquoise", "turquoise","turquoise", "brown", "brown","brown", "brown","brown", "yellow", "yellow","yellow", "yellow","yellow")

barplot(t, beside=T, col=colx,border="grey30", las=2, cex.axis=0.7, ylim=c(-.6,.6))
abline(v=1.5, lty=3, col="black")
abline(v=2.5, lty=3, col="grey30")
abline(v=3.5, lty=3, col="grey50")
abline(v=4.5, lty=3, col="grey70")
abline(v=5.5, lty=3, col="grey80")

abline(v=7.5, lty=3, col="black")
abline(v=8.5, lty=3, col="grey30")
abline(v=9.5, lty=3, col="grey50")
abline(v=10.5, lty=3, col="grey70")
abline(v=11.5, lty=3, col="grey80")

abline(v=13.5, lty=3, col="black")
abline(v=14.5, lty=3, col="grey30")
abline(v=15.5, lty=3, col="grey50")
abline(v=16.5, lty=3, col="grey70")
abline(v=17.5, lty=3, col="grey80")

abline(v=19.5, lty=3, col="black")
abline(v=20.5, lty=3, col="grey30")
abline(v=21.5, lty=3, col="grey50")
abline(v=22.5, lty=3, col="grey70")
abline(v=23.5, lty=3, col="grey80")

