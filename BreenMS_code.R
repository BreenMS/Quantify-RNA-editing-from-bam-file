#################################################################
################# Load the Libraries ############################
#################################################################
library(DESeq)
library(RColorBrewer)
library(gplots)
library(ggplot2)

library(RUVSeq)
library(EDASeq)
library(limma)

####READ IN RAW EXPRESSION DATA
x <-read.ilmn(files="Drakenstein_raw.txt")
head(x$E)
dim(x$E)

####READ IN META-DATA
targets <- read.delim("Targets.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(targets)
Batch <- factor((targets$Batch))
Time <- factor((targets$Time))
Strain <- factor((targets$Strain))
Labels <- factor((targets$Combined))
colors <- c("firebrick1","brown","steelblue4", "orange")

###PERFORM BACKGROUND CORRECTION,QUANTILE NORMALIZATION AND LOG2 TRANSFORMATION WITH NEQC
y <- neqc(x)
dim(y)

###REMOVE LOWLY EXPRESSED GENES (NON-SPECIFIC FILTERING)
###We keep probes that are expressed in at least 50 arrays with detection p-value < 0.05:
expressed <- rowSums(y$other$Detection < 0.05) >= 50
y <- y[expressed,]
dim(y)

## MEDIAN CONDENSE PROBES MATACHING TO THE SAME GENE SYMBOL
##Begin median summary

write.csv(y, "Norm_expression.csv",row.names=FALSE)
read.csv("Norm_expression.csv",header=TRUE)->mat_med
dim(mat_med)

#####MEDIAN SUMMARY
require(plyr)
mat_medians <- ddply(.data = mat_med,
                     .variables = "SYMBOL",
                     function(df){
                       # Remember, exclude the "SYMBOL" column in the median calculation!
                       apply(X = df[, !colnames(df) %in% "SYMBOL"], # apply to each df subset
                             MARGIN = 2, 
                             FUN = median)
                     })

nrow(mat_medians)
length(unique(mat_med$SYMBOL))
write.table(mat_medians, "Norm_median.txt", sep="\t")

##OUTLIER ANALYSIS by PCA
library(rgl)
matrix <-  read.table("Norm_median_clean.txt",header=TRUE, row.names=1, sep="\t")
head(matrix)
dim(matrix)

#CREATE 3D INTERACTIVE PCA PLOT
tmatrix<-t(matrix)
pcs<-prcomp(tmatrix)
summary(pcs)
PC1<-pcs$x[,1]
PC2<-pcs$x[,2]
PC3<-pcs$x[,3]

PCA_details <- cbind(PC1, PC2, PC3)

plot3d(PC1, PC2, PC3, main="Condition", xlab="PC1", ylab="PC2", zlab="PC3",
       ,box=FALSE, size=1, type="s")

labels<-rownames(tmatrix)
text3d(PC1+6, PC2+6, PC3+6,text=labels, font=1, cex=1)
 
#ADD ELLIPSOIDS AS STANDARD DEVIATIONS
mean.vec <- c(mean(PC1), mean(PC2), mean(PC3))
allcomp <- cbind(PC1, PC2, PC3)
sigma <- cov(allcomp) #sigma <-  c(sd(PC1), sd(PC2), sd(PC3))
plot3d( ellipse3d(x = sigma, centre=mean.vec,scale = c(1.5, 1.5, 1.5),), col="PeachPuff", alpha=0.50, add = TRUE, level=0.95, smooth=TRUE)
plot3d( ellipse3d(x = sigma, centre=mean.vec,scale = c(2, 2, 2),), col="GhostWhite", alpha=0.50, add = TRUE, level=0.95,smooth=TRUE)



#UNSUPERVISED CLUSTERING WITH PEARSONS CORRELATION
library(WGCNA)
library(cluster)
options(stringsAsFactors  =  FALSE)

#DATA INPUT AND FILTERING
M1 <- read.delim("Norm_median.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(M1) =names(M1)
dataExpr0<-as.data.frame(t(M1))
head(dataExpr0)
dim(dataExpr0)
gsg=goodSamplesGenes(dataExpr0,verbose=3)
gsg$allOK

nGenes = ncol(dataExpr0)
nSamples = nrow(dataExpr0)

#APPLY Meta-Data
traitData <- read.delim("Targets.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData)
head(traitData)
rowsExpr <- rownames(dataExpr0)
traitRows <- match(rowsExpr,traitData$SID.1)
datTraits = traitData[traitRows, -1];
rownames (datTraits) = traitData[traitRows, 1];
table(rownames(datTraits)==rownames(dataExpr0)) # EVERYTHING OK?

#SAMPLE NETWORK BASED ON SQUARED EUCLIDEAN DISTANCE
A=adjacency(t(dataExpr0))
k=as.numeric(apply(A,2,sum))-1
Z.k=scale(k)

#DESIGNATE SAMPLES AS OUTLIERS WITH RED COLORS
thresholdZ.k= -2 # 
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")

#CLUSTER ALL SAMPLES
#pearsons distance = as.dist(1-cor(t(x)))
sampleTree = hclust(as.dist(1-A), method = "average")
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits),sep="")
datColors=data.frame(Outliers=outlierColor,traitColors)

plotDendroAndColors(sampleTree,groupLabels=names(datColors), cex.dendroLabels = 0.6, las=1,
colors=datColors,main="Clustergram and trait heatmap")


M1 <- read.delim("Norm_median_clean.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)

###QUALITY CONTROL PLOTS
par(mfrow=c(2,2))
boxplot(log2(x$E),range=0,ylab="log2 intensity - Raw", las=2, cex.axis=0.4, main="Raw")  #boxplots log raw expression
hist(log2(x$E),main="Raw", xlab="log2 intensity - Raw") #histogram log raw expression
plotRLE(log2(x$E),ylim=c(-0.4,0.4), cex.axis=0.5, las=2, outline=FALSE, ylab="RLE", main="Raw") # RLE plot log raw expression
#legend("top", fill =colors, cex=0.6, legend = c("BN-SS", "F344", "WKY", "LEW"), horiz=T)
rawq = as.matrix(sapply(log2(x$E), as.numeric)) #set up QQ plot
qqnorm(rawq, pch=20, las=1, cex.axis=0.8, col="grey30", main="Raw");qqline(rawq, lwd=1.5, lty=3) #QQ plot raw expression
boxplot((M1),range=0,ylab="log2 intensity - Norm", las=2, cex.axis=0.4,main="Norm") #boxplots log normalized expression
hist(M1),main="Norm", xlab="log2 intensity - Norm") #histogram log normalized expression
plotRLE((y$E),ylim=c(-0.2,0.2), cex.axis=0.5, las=2, outline=FALSE, ylab="RLE", main="Norm")
normq = as.matrix(sapply(log2(y$E), as.numeric))  
qqnorm(normq, pch=20, las=1, cex.axis=0.8, col="grey30", main="Norm");qqline(normq, lwd=1.5, lty=3) #QQplot normalized expression

par(mfrow=c(2,2))
plotMDS(y,labels=targets$Batch, main="Batch",col=colors[Batch], cex=0.5)
plotMDS(y,labels=targets$Time, main="Time",col=colors[Time],cex=0.5)
plotMDS(y,labels=targets$Strain, main="Strain",col=colors[Strain],cex=0.5)
plot(hclust(dist(t(y$E), method="euclidean"),method="average"), cex=0.5)


###QUALITY CONTROL PLOTS - variancePartion
library(variancePartition)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

exprsFile <- file.path("Norm_median_clean.txt") #load normalized matrix
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
dim(exprs)

info <- read.delim("Targets.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(info)

form <- ~ RIN + (1|TE) + (1|Alcohol) + (1|Nicotine) + (1|Ethnicity) + (1|Delivery) + (1|Gender) + (1|PrevMed) + (1|PregMed) + (1|Batch)  + (1|Depression) + (1|PTSD) 

varPart <- fitExtractVarPartModel(exprs, form, info )
vp <- sortCols( varPart )
plotVarPart(vp,label.angle=50)


##OUTLIER ANALYSIS by PCA
library(rgl)
matrix <-  read.table("Norm_median_clean.txt",header=TRUE, row.names=1, sep="\t")
head(matrix)
dim(matrix)

targets <- read.delim("Targets.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
Group <- factor((targets$Dx))

#CREATE 3D INTERACTIVE PCA PLOT
tmatrix<-t(matrix)
pcs<-prcomp(tmatrix)
summary(pcs)
PC1<-pcs$x[,1]
PC2<-pcs$x[,2]
PC3<-pcs$x[,3]

colors=c("steelblue4", "steelblue3", "aquamarine4","indianred4","salmon")
colfunc <- colorRampPalette(c("aquamarine4", "white","steelblue4"))

plot(PC1, PC2, col=colors[Group], pch=16, cex=0.8, las=1, cex.axis=0.7, xlab="PC1: 16% variance", ylab="PC2: 12% variance")
abline(h=0, lty=2, lwd=0.5)
abline(v=0, lty=2, lwd=0.5)
legend("bottomleft", fill =colors, cex=0.6, legend = c("Control", "TE", "Dep", "PTSD","PTSD/Dep"), bty="n")

###QUALITY CONTROL PLOTS - PC1 association with meta-data
matrix <-  read.table("Norm_median_clean.txt",header=TRUE, row.names=1, sep="\t")
nSamples = ncol(matrix)
tmatrix<-t(matrix)
pcs<-prcomp(tmatrix)
summary(pcs)
PC1<-pcs$x[,1]

#Load meta-data
targets <- read.delim("Targets.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")

moduleTraitCor = cor(PC1, targets, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
targets <- read.delim("testing.txt", stringsAsFactors=FALSE, row.names=1, header=F, sep="\t")
barplot((-log10(targets$V2)), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(targets))
abline(v=-log10(0.05), col="red")


#####DIFFERENTIAL EXPRESSION ANALYSIS WITH LIMMA ######
library(limma)
exprs <- read.delim("PTSD_Dep_Matrix_sum_c.txt", check.names=FALSE, stringsAsFactors=FALSE, header=T, row.names=1)
head(exprs)
dim(exprs)

targets <- read.delim("PTSD_Dep_Matrix_sum_c_targets.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(targets)
Group <- factor(targets$PTSD, levels=c("A","B"))
Group
design <- model.matrix(~0+Group)
colnames(design) <- c("A","B")
colnames(design)

Nicotine <- factor(targets$Nicotine)
Alcohol <- factor(targets$Alcohol)
Gender <- factor(targets$Gender)
Delivery <- factor(targets$Delivery)
Ethnicity <- factor(targets$Ethnicity)
Chip <- (targets$Chip)
RIN <- (targets$RIN)
Batch <- factor(targets$Batch)

design <- model.matrix(~Group+Nicotine+Alcohol+RIN+Batch+Gender+Delivery+Ethnicity)
colnames(design)

fit <- lmFit(exprs, design)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")

Diff<- topTable(fit2, coef="GroupB", n=25000)
write.table(Diff, file="PTSD_Dep_c_DEG_corrected.txt", quote=FALSE, row.names=TRUE)



####SVA CORRECTION####
library(sva)
library(limma)

exprsFile <- file.path("PTSD_Matrix_sum.txt")
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="", row.names=1, as.is=TRUE))
pheno <- read.delim("PTSD_Matrix_sum_targets.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(pheno)

mod = model.matrix(~PTSD, data=pheno) #Cancer as factor in Pheno File, treated as a factor variable
colnames(mod) <- c("A","B")
mod0 = model.matrix(~1,data=pheno)

#Find latent factors that need to be estimated
n.sv = num.sv(exprs,mod,method="be",B = 50)
n.sv

#Estimate Surrogate Variables
svobj = sva(exprs,mod,mod0,n.sv=n.sv)
modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(exprs,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

#SVA through limma
fit = lmFit(exprs,modSv)
contrast.matrix <- cbind("PTSD"=c(1,-1,rep(0,svobj$n.sv)))

fit2 = contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)

#output the differentially expressed genes from Comparison 1
topTable(fit2, coef="PTSD", adjust="BH")
DEG <- topTable(fit2, coef="PTSD", adjust="BH",n=30000) #collect 24K genes
write.table(DEG, file="PTSD_DGE_sva.txt", quote=FALSE, row.names=TRUE)


##Cluster matrix
exprs <- read.delim("Depression_Matrix_sum.txt", check.names=FALSE, stringsAsFactors=FALSE, header=T, row.names=1)
dim(exprs)
plot(hclust(dist(t(exprs), method="euclidean"),method="ward.D"), cex=0.5)

dists = dist(t(exprs))
mat = as.matrix( dists )
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, col = hmcol, trace="none", margin=c(7, 7), cexRow=0.5, cexCol=0.7)


##HEATMAP
library(gplots)
library(RColorBrewer);
marker <- as.matrix(read.table("list2.txt", header=TRUE, row.names=1, sep="\t", as.is=TRUE))
dim(marker)

dist.pear <- function(x) as.dist(1-cor(t(x)))
dist.euc <- function(x) dist(x,method = "euclidean")
hclust.ave <- function(x) hclust(x, method="average")
hclust.ward <- function(x) hclust(x, method="ward.D2")
hclust.complete <- function(x) hclust(x, method="complete")

d <- dist.pear(marker)
fit <- hclust.ave(d)
clusters <- cutree(fit, h=0.8)
nofclust.height <-  length(unique(as.vector(clusters)));

selcol2 <- colorRampPalette(brewer.pal(9,"Set1"))
clustcol.height = selcol2(nofclust.height);
clustcol.height[clusters]

par(cex.main=0.8)
heatmap.2(marker,scale="row",hclustfun=hclust.ave, distfun=dist.pear,cexRow=0.7,margins = c(10, 10),cexCol=0.7, #Colv=NA,
key=FALSE, density.info="none",symkey=TRUE, trace="none", RowSideColors=clustcol.height[clusters],
col=colorRampPalette(c("navy", "white", "red"))(1024))


##Create Jaccard clustering
exprsFile <- file.path("list2.txt")
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
dim(exprs)

library(vegan)
dist.mat<-vegdist(t(exprs),method="jaccard")
clust.res<-hclust(dist.mat)
plot(clust.res, hang=0.6, cex=.6, main="", las=1, cex.axis=0.7)
plot(hclust(dist(t(exprs), method="euclidean"),method="ward.D"), las=2, cex=0.7, cex.axis=0.7, main="")

library(dendextend)
dend <- as.dendrogram(clust.res)
dend <- rotate(dend, 1:7)
dend <- color_branches(dend, k=2)
dend <- hang.dendrogram(dend,hang_height=0.1)
dend <- set(dend, "labels_cex", 0.6)
par(bg = "white")
labels_colors(dend) <- "white"
plot(dend, main = "", horiz=T,las=1,nodePar = list(cex = .007),
     cex.axis=0.6,xlim=c(1, 0.8), lwd=5)
legend("topleft", legend =c("Combat trauma", "Interpersonal trauma"))




#correspondence of logFC
matrix1 <- read.delim("list2.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
head(matrix1)

par(mfrow=c(2,4))
verboseScatterplot(matrix1$PTSD_c, matrix1$TE, pch=19, cex=0.8, col=matrix1$Col1,las=1, cex.axis=0.8, abline=T, ylab="TE",xlab="PTSD",abline.lty = 2)
abline(h=0,lty=3,col="grey40")
abline(v=0,lty=3,col="grey40")
verboseScatterplot(matrix1$PTSD_c, matrix1$PTSD.Dep_c, pch=19, cex=0.8, col=matrix1$Col2,las=1, cex.axis=0.7, abline=T, ylab="PTSD.Dep",xlab="PTSD",abline.lty = 2)
abline(h=0,lty=3,col="grey40")
abline(v=0,lty=3,col="grey40")
verboseScatterplot(matrix1$PTSD_c, matrix1$Depression, pch=19, cex=0.8, col=matrix1$Col3,las=1, cex.axis=0.7, abline=T, ylab="Dep",xlab="PTSD",abline.lty = 2)
abline(h=0,lty=3,col="grey40")
abline(v=0,lty=3,col="grey40")
verboseScatterplot(matrix1$PTSD.Dep_c, matrix1$Depression, pch=19, cex=0.8, col=matrix1$Col4,las=1, cex.axis=0.7, abline=T, ylab="Dep",xlab="PTSD.Dep",abline.lty = 2)
abline(h=0,lty=3,col="grey40")
abline(v=0,lty=3,col="grey40")



matrix <- read.delim("list3.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
head(matrix)
verboseScatterplot(matrix1$PTSD_c, matrix1$TE, pch=19, cex=0.8, col="grey",las=1, cex.axis=0.7, abline=T, abline.color="grey40",ylab="TE",xlab="PTSD", ylim=c(-.5,.5),bty='n',fdrak,abline.lty = 2)
abline(h=0,lty=3,col="grey30")
abline(v=0,lty=3,col="grey30")
par(new=TRUE)
verboseScatterplot(matrix$TE_TE, matrix$TE_PTSD, pch=19, cex=0.8, col=matrix$Color,las=1, cex.axis=0.7, abline=T,abline.color="black", ylab="TE",xlab="PTSD", ylim=c(-.5,.5),abline.lty = 2)

t=cbind(matrix$PTSD, matrix$Dep, matrix$TE, matrix$PTSD.Dep)

plot(matrix$PTSD, matrix$Dep, pch=20, cex=0.6, las=1, cex.axis=0.7)#,col=matrix$Color)

matrix <- read.delim("list3.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
head(matrix)
plot(matrix$logFC, -log10(matrix$P.Value), pch=20, cex=0.5, las=1, cex.axis=0.7, ylim=c(0,4.5))#,col=matrix$Color)
abline(v=0,lty=3)
abline(h=-log10(0.005),lty=2, col="red")
abline(h=-log10(0.001),lty=2, col="orange")


matrix <- read.delim("list3.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=F)
head(matrix)
colfunc <- colorRampPalette(c("salmon", "white"))
barplot(-log10(matrix$V2),horiz=T,col=colfunc(5), cex.axis=0.7)
abline(v=-log10(0.05),col="black")






matrix <- read.delim("list3.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
head(matrix)

library(ggplot2)
ggplot(data=matrix, aes(x=Measure, y=Value, fill=Group)) +
    geom_bar(position=position_dodge(), stat="identity") +
       # geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=.1,  position=position_dodge(.9))+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                                    #axis.line = element_line(colour = "black",size=0.5),
                                    axis.text.x = element_text(colour="black", size=7, angle=45, hjust=1),
                                    axis.text.y = element_text(colour="black", size=7),
                                    axis.ticks.x = element_line(colour="black"),
                                    axis.ticks.y = element_line(colour="black"),
                                    #panel.grid.major.y = element_line(colour = "black", linetype="dashed"),
                                    panel.background = element_rect(fill = "grey90", colour = "white")) +
                                    scale_fill_manual(values=c("blue", "orange","red")) 

matrix <- as.data.frame(read.table("list3.txt", header=T, row.names=NULL, sep="\t", as.is=TRUE))
head(matrix)

tgc <-summarySE(matrix, measurevar="Value", groupvars=c("Group", "Measure"))

tgc2 <- tgc
tgc2$Group <- factor(tgc2$Group)
tgc2$Measure <- factor(tgc2$Measure)

ggplot(tgc2, aes(x=Measure, y=Value, fill=Group)) + 
geom_boxplot(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=Value-se, ymax=Value+se), width=.1,  position=position_dodge(.9))+
    scale_fill_manual(values=c("grey80", "salmon","lightcoral"))#+
    #scale_x_discrete(breaks=c("a", "y", "c"),labels=c("B","Month 6", "Month 24"))


matrix <- as.data.frame(read.table("list3.txt", header=T, row.names=NULL, sep="\t", as.is=TRUE))
head(matrix)

matrix$Measure <- factor(matrix$Measure, levels=unique(matrix$Measure))
matrix$Group <- factor(matrix$Group, levels=unique(matrix$Group))

ggplot(matrix, aes(x=Measure, y=Value, fill=Group)) +
  geom_boxplot(position=position_dodge(0.7)) +   coord_cartesian(ylim=c(70,15)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                                    #axis.line = element_line(colour = "black",size=2),
                                    axis.text.x = element_text(colour="black", size=11),#, angle=45, hjust=1),
                                    axis.text.y = element_text(colour="black", size=11),
                                    axis.ticks.x = element_line(colour="black"),
                                    axis.ticks.y = element_line(colour="black"),
                                    panel.grid.major.y = element_line(colour = "grey90"),
                                    panel.grid.major.x = element_line(colour = "grey90"),
                                    panel.background = element_rect(fill = "white", colour = "grey90"))















exprsFile <- file.path("PTSD_Matrix_sum.txt") 
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
dim(exprs)

targets <- read.delim("PTSD_Matrix_sum_targets.txt", check.names=FALSE, stringsAsFactors=FALSE) #load metadata
head(targets)

#label meta-data
Combined <- factor(targets$Combined)

#create design matrix
design <- model.matrix(~0+Combined)
colnames(design) <- c("Control_hr", "Control_sn", "Control_ss", "Control_su", "CRS_hr", "CRS_sn", "CRS_ss", "CRS_su", "PCRS_hr", "PCRS_sn", "PCRS_ss", "PCRS_su")
colnames(design)
fit <- lmFit(exprs, design)
#plotSA(fit, main="Probe-level")

#specify the comparisons you want to make
control.matrix <- makeContrasts(
Control1=Control_ss-Control_sn,
Control2=Control_ss-Control_su,
Control3=Control_ss-Control_hr,
Control4=Control_sn-Control_su,
Control5=Control_sn-Control_hr,
Control6=Control_su-Control_hr,

CRS1=CRS_ss-CRS_sn,
CRS2=CRS_ss-CRS_su,
CRS3=CRS_ss-CRS_hr,
CRS4=CRS_sn-CRS_su,
CRS5=CRS_sn-CRS_hr,
CRS6=CRS_su-CRS_hr,

PCRS1=PCRS_ss-PCRS_sn,
PCRS2=PCRS_ss-PCRS_su,
PCRS3=PCRS_ss-PCRS_hr,
PCRS4=PCRS_sn-PCRS_su,
PCRS5=PCRS_sn-PCRS_hr,
PCRS6=PCRS_su-PCRS_hr,
levels=design)

#use linear models to make those comparisons
fit2 <- contrasts.fit(fit, control.matrix)
fit2 <- eBayes(fit2)

#Check out results from Comparison1
Control1 <- topTable(fit2, coef="Control1", adjust="BH",n=24000)
Control2 <- topTable(fit2, coef="Control2", adjust="BH",n=24000)
Control3 <- topTable(fit2, coef="Control3", adjust="BH",n=24000)
Control4 <- topTable(fit2, coef="Control4", adjust="BH",n=24000)
Control5 <- topTable(fit2, coef="Control5", adjust="BH",n=24000)
Control6 <- topTable(fit2, coef="Control6", adjust="BH",n=24000)

CRS1 <- topTable(fit2, coef="CRS1", adjust="BH",n=24000)
CRS2 <- topTable(fit2, coef="CRS2", adjust="BH",n=24000)
CRS3 <- topTable(fit2, coef="CRS3", adjust="BH",n=24000)
CRS4 <- topTable(fit2, coef="CRS4", adjust="BH",n=24000)
CRS5 <- topTable(fit2, coef="CRS5", adjust="BH",n=24000)
CRS6 <- topTable(fit2, coef="CRS6", adjust="BH",n=24000)

PCRS1 <- topTable(fit2, coef="PCRS1", adjust="BH",n=24000)
PCRS2 <- topTable(fit2, coef="PCRS2", adjust="BH",n=24000)
PCRS3 <- topTable(fit2, coef="PCRS3", adjust="BH",n=24000)
PCRS4 <- topTable(fit2, coef="PCRS4", adjust="BH",n=24000)
PCRS5 <- topTable(fit2, coef="PCRS5", adjust="BH",n=24000)
PCRS6 <- topTable(fit2, coef="PCRS6", adjust="BH",n=24000)

write.table(Control1, file="DEG_Control_ss_v_sn.txt", quote=FALSE, row.names=TRUE)
write.table(Control2, file="DEG_Control_ss_v_su.txt", quote=FALSE, row.names=TRUE)
write.table(Control3, file="DEG_Control_ss_v_hr.txt", quote=FALSE, row.names=TRUE)
write.table(Control4, file="DEG_Control_sn_v_su.txt", quote=FALSE, row.names=TRUE)
write.table(Control5, file="DEG_Control_sn_v_hr.txt", quote=FALSE, row.names=TRUE)
write.table(Control6, file="DEG_Control_su_v_hr.txt", quote=FALSE, row.names=TRUE)

write.table(CRS1, file="DEG_CRS_ss_v_sn.txt", quote=FALSE, row.names=TRUE)
write.table(CRS2, file="DEG_CRS_ss_v_su.txt", quote=FALSE, row.names=TRUE)
write.table(CRS3, file="DEG_CRS_ss_v_hr.txt", quote=FALSE, row.names=TRUE)
write.table(CRS4, file="DEG_CRS_sn_v_su.txt", quote=FALSE, row.names=TRUE)
write.table(CRS5, file="DEG_CRS_sn_v_hr.txt", quote=FALSE, row.names=TRUE)
write.table(CRS6, file="DEG_CRS_su_v_hr.txt", quote=FALSE, row.names=TRUE)

write.table(PCRS1, file="DEG_PCRS_ss_v_sn.txt", quote=FALSE, row.names=TRUE)
write.table(PCRS2, file="DEG_PCRS_ss_v_su.txt", quote=FALSE, row.names=TRUE)
write.table(PCRS3, file="DEG_PCRS_ss_v_hr.txt", quote=FALSE, row.names=TRUE)
write.table(PCRS4, file="DEG_PCRS_sn_v_su.txt", quote=FALSE, row.names=TRUE)
write.table(PCRS5, file="DEG_PCRS_sn_v_hr.txt", quote=FALSE, row.names=TRUE)
write.table(PCRS6, file="DEG_PCRS_su_v_hr.txt", quote=FALSE, row.names=TRUE)


control.matrix <- makeContrasts(
SS1=CRS_ss-Control_ss,
SS2=PCRS_ss-Control_ss,
SN1=CRS_sn-Control_sn,
SN2=PCRS_sn-Control_sn,
HR1=CRS_hr-Control_hr,
HR2=PCRS_hr-Control_hr,
SU1=CRS_su-Control_su,
SU2=PCRS_su-Control_su,
levels=design)

control.matrix <- makeContrasts(
"CRS_su-Control_su",
"PCRS_su-Control_su",
"PCRS_su-CRS_su",
levels=design)

control.matrix <- makeContrasts(
one=CRS_su-Control_su,
two=PCRS_su-Control_su,
three=PCRS_su-CRS_su,
levels=design)


fit2 <- contrasts.fit(fit, control.matrix)
fit2 <- eBayes(fit2)

a <- topTable(fit2, coef="one", adjust="BH",n=24000)
b <- topTable(fit2, coef="two", adjust="BH",n=24000)
c <- topTable(fit2, coef="three", adjust="BH",n=24000)
write.table(a, file="DEG_1.txt", quote=FALSE, row.names=TRUE)
write.table(b, file="DEG_2.txt", quote=FALSE, row.names=TRUE)
write.table(c, file="DEG_3.txt", quote=FALSE, row.names=TRUE)

test <- topTable(fit2, adjust="BH",n=24000)
write.table(test, file="TestingDEG.txt", row.names=T)

#Check out results from Comparison1
#topTable(fit2, coef="SS1", adjust="BH")
WKY_CRS <- topTable(fit2, coef="SS1", adjust="BH",n=24000)
WKY_PCRS <- topTable(fit2, coef="SS2", adjust="BH",n=24000)
F344_CRS <- topTable(fit2, coef="SN1", adjust="BH",n=24000)
F344_PCRS <- topTable(fit2, coef="SN2", adjust="BH",n=24000)
BNSS_CRS <- topTable(fit2, coef="HR1", adjust="BH",n=24000)
BNSS_PCRS <- topTable(fit2, coef="HR2", adjust="BH",n=24000)
LEW_CRS <- topTable(fit2, coef="SU1", adjust="BH",n=24000)
LEW_PCRS <- topTable(fit2, coef="SU2", adjust="BH",n=24000)

write.table(WKY_CRS, file="DEG_CRS_WKY.txt", quote=FALSE, row.names=TRUE)
write.table(WKY_PCRS, file="DEG_PCRS_WKY.txt", quote=FALSE, row.names=TRUE)
write.table(F344_CRS, file="DEG_CRS_F344.txt", quote=FALSE, row.names=TRUE)
write.table(F344_PCRS, file="DEG_PCRS_F344.txt", quote=FALSE, row.names=TRUE)
write.table(BNSS_CRS, file="DEG_CRS_BNSS.txt", quote=FALSE, row.names=TRUE)
write.table(BNSS_PCRS, file="DEG_PCRS_BNSS.txt", quote=FALSE, row.names=TRUE)
write.table(LEW_CRS, file="DEG_CRS_LEW.txt", quote=FALSE, row.names=TRUE)
write.table(LEW_PCRS, file="DEG_PCRS_LEW.txt", quote=FALSE, row.names=TRUE)



Control1_2 <- cbind(Genes = rownames(Control1), Control1)
DE_Cntl1 <- factor(Control1_2$Genes[Control1_2$adj.P.Val < 0.05])
DE_Cntl1 <- data.frame(DE_Cntl1)
DE_Cntl1$DE_Cntl1 <- 1
nrow(DE_Cntl1)

Control2_2 <- cbind(Genes = rownames(Control2), Control2)
DE_Cntl2 <- factor(Control2_2$Genes[Control2_2$adj.P.Val < 0.05])
DE_Cntl2 <- data.frame(DE_Cntl2)
DE_Cntl2$DE_Cntl2 <- 1
nrow(DE_Cntl2)

Control3_2 <- cbind(Genes = rownames(Control3), Control3)
DE_Cntl3 <- factor(Control3_2$Genes[Control3_2$adj.P.Val < 0.05])
DE_Cntl3 <- data.frame(DE_Cntl3)
DE_Cntl3$DE_Cntl3 <- 1
nrow(DE_Cntl3)

Control4_2 <- cbind(Genes = rownames(Control4), Control4)
DE_Cntl4 <- factor(Control4_2$Genes[Control4_2$adj.P.Val < 0.05])
DE_Cntl4 <- data.frame(DE_Cntl4)
DE_Cntl4$DE_Cntl4 <- 1
nrow(DE_Cntl4)

Control5_2 <- cbind(Genes = rownames(Control5), Control5)
DE_Cntl5 <- factor(Control5_2$Genes[Control5_2$adj.P.Val < 0.05])
DE_Cntl5 <- data.frame(DE_Cntl5)
DE_Cntl5$DE_Cntl1 <- 1
nrow(DE_Cntl5)

Control6_2 <- cbind(Genes = rownames(Control6), Control6)
DE_Cntl6 <- factor(Control6_2$Genes[Control6_2$adj.P.Val < 0.05])
DE_Cntl6 <- data.frame(DE_Cntl6)
DE_Cntl6$DE_Cntl6 <- 1
nrow(DE_Cntl6)


rawdata <- read.delim("list.txt", check.names=FALSE, stringsAsFactors=FALSE, row.names=1)
head(rawdata)
dim(rawdata)

par(mfrow=c(2,2))
verboseScatterplot(rawdata$WKY1,rawdata$WKY2, xlab="CRS logFC", ylab="PCRS logFC", main="WKY\n",cex = 0.8, cex.axis = 0.8, cex.lab = 0.8, cex.main = 1, abline=TRUE, abline.color=2,abline.lty=2,ylim=c(-3,3), xlim=c(-3,3))
abline(h=0.5)
abline(h=-0.5)
abline(v=-0.5)
abline(v=0.5)
verboseScatterplot(rawdata$F3441,rawdata$F3442, xlab="CRS logFC", ylab="PCRS logFC", main="F344\n",cex = 0.8, cex.axis = 0.8, cex.lab = 0.8, cex.main = 1, abline=TRUE, abline.color=2,abline.lty=2,ylim=c(-3,3), xlim=c(-3,3))
abline(h=0.5)
abline(h=-0.5)
abline(v=-0.5)
abline(v=0.5)
verboseScatterplot(rawdata$BNSS1,rawdata$BNSS2, xlab="CRS logFC", ylab="PCRS logFC", main="BN-SS\n",cex = 0.8, cex.axis = 0.8, cex.lab = 0.8, cex.main = 1, abline=TRUE, abline.color=2,abline.lty=2,ylim=c(-3,3), xlim=c(-3,3))
abline(h=0.5)
abline(h=-0.5)
abline(v=-0.5)
abline(v=0.5)
verboseScatterplot(rawdata$LEW1,rawdata$LEW2, xlab="CRS logFC", ylab="PCRS logFC", main="LEW\n",cex = 0.8, cex.axis = 0.8, cex.lab = 0.8, cex.main = 1, abline=TRUE, abline.color=2,abline.lty=2,ylim=c(-3,3), xlim=c(-3,3))
abline(h=0.5)
abline(h=-0.5)
abline(v=-0.5)
abline(v=0.5)



plot(hclust(dist(t(rawdata), method="euclidean"),method="average"), las=2, cex=0.7, cex.axis=0.7, main="Similarity based on logFC\n ns v PCRS")

plot(rawdata$WKY1, rawdata$WKY2, xlab="CRS", ylab="PCRS", ylim=c(-3,3), xlim=c(-3,3), col="white")
abline(lm(rawdata$WKY2~rawdata$WKY1), col="gray90",lwd=2) # regression line (y~x) 
abline(lm(rawdata$F3442~rawdata$F3441), col="gray40",lwd=2) # regression line (y~x) 
abline(lm(rawdata$BNSS2~rawdata$BNSS1), col="turquoise4",lwd=2) # regression line (y~x) 
abline(lm(rawdata$LEW2~rawdata$LEW1), col="turquoise", lwd=2) # regression line (y~x) 


points(rawdata$F3441,rawdata$F3442)
points(rawdata$BNSS1,rawdata$BNSS2)
points(rawdata$LEW1,rawdata$LEW2)













#################################################################
################# LOAD RAW MATRIX FILE ##########################
#################################################################

rawdata <- read.delim("AMG_Matrix.txt", check.names=FALSE, stringsAsFactors=FALSE, row.names=1)
head(rawdata)
dim(rawdata)

#################################################################
################# Apply non-specific filtering ##########################
#################################################################

filter<-apply(rawdata,1, function(x) length(x[x>1])>=10) #require more than 5 reads in at least 10 samples
filtered <- rawdata[filter,]
dim(filtered)

#################################################################
################# Input meta-data file ##########################
#################################################################

targets <- read.delim("Targets.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(targets)
Group <- factor((targets$Group))

#################################################################
################# Create an Eset (expressionSet) ##########################
#################################################################

set<-newSeqExpressionSet(as.matrix(filtered),
phenoData=data.frame(Group,row.names=colnames(filtered)))

#################################################################
################# Perform Upper Quantile Normalization ##########################
#################################################################

set_uq <- betweenLaneNormalization(set,which="upper")

# extract and transpose expression matrix
upper_quartile_matrix <- normCounts(set_uq)
upper_quartile_matrix_t =t(upper_quartile_matrix)
write.table(upper_quartile_matrix, "Blood_Matrix_UQ.txt")

#################################################################
################# Input some color codes for plots ##########################
#################################################################

colors <- c("firebrick1","brown","steelblue4", "steelblue2")
colnames(pData(set_uq))
table(pData(set)[,"Group"]);
table(pData(set_uq)[,"Group"]);

#################################################################
################# create diagnostic plots ##########################
#################################################################

par(mfrow=c(2,3))
boxplot(counts(set),col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main=" Raw Counts", outline=FALSE,ylab="Raw Expression")
boxplot(set_uq,col=colors[Group], las=2, cex=0.5, cex.axis=0.5, main="Norm Counts (UQ)", outline=FALSE, ylab="log UQ Expression",ylim=c(0,9))
hist(log(counts(set_uq)), main="Histogram of log Norm UQ Counts", xlab="log UQ Expression")
plotRLE(set_uq,ylim=c(-2,2), cex.axis=0.5, las=2,col=colors[Group], outline=FALSE, ylim="RLE", main="RLE on Norm Counts UQ")
plotPCA(set_uq,cex=0.6,col=colors[Group], main="PCA UQ", ylim=c(-0.3, 0.5), xlim=c(-0.25, 0.25), las=2)
plot(hclust(dist(upper_quartile_matrix_t, method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7)
meanVarPlot(set_uq,log=TRUE,ylim=c(0,20), main="UQ Over-Dispersion")


#################################################################
################# Differential gene expression analysis ##########################
#################################################################
library(limma)
exprsFile <- file.path("AMG_Matrix_voom.txt") #load normalized matrix
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
dim(exprs)

targets <- read.delim("Targets.txt", check.names=FALSE, stringsAsFactors=FALSE) #load metadata
head(targets)

#label meta-data
Group <- factor(targets$Combined, levels=c("Control_BNSS", "Control_F344", "Control_LEW", "Control_WKY", "CRS_BNSS", "CRS_F344", "CRS_LEW", "CRS_WKY", "PCRS_BNSS", "PCRS_F344", "PCRS_LEW", "PCRS_WKY"))

#create design matrix
design <- model.matrix(~0+Group)
colnames(design) <- c("Control_BNSS", "Control_F344", "Control_LEW", "Control_WKY", "CRS_BNSS", "CRS_F344", "CRS_LEW", "CRS_WKY", "PCRS_BNSS", "PCRS_F344", "PCRS_LEW", "PCRS_WKY")
colnames(design)
fit <- lmFit(exprs, design)

#specify the comparisons you want to make
contrast.matrix <- makeContrasts(
Comparison1=Control_WKY-PCRS_WKY,
Comparison2=Control_F344-PCRS_F344,
Comparison3=Control_LEW-PCRS_LEW,
Comparison4=Control_BNSS-PCRS_BNSS,
levels=design)

#use linear models to make those comparisons
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Check out results from Comparison1
topTable(fit2, coef="WKY", adjust="BH")
WKY_DEG <- topTable(fit2, coef="WKY", adjust="BH",n=24000)
write.table(WKY_DEG, file="AMG_DEG_WKY.txt", quote=FALSE, row.names=TRUE)

#Check out results from Comparison2
topTable(fit2, coef="F344", adjust="BH")
F344_DEG <- topTable(fit2, coef="F344", adjust="BH",n=24000)
write.table(F344_DEG, file="AMG_DEG_F344.txt", quote=FALSE, row.names=TRUE)

#Check out results from Comparison3
topTable(fit2, coef="LEW", adjust="BH")
LEW_DEG <- topTable(fit2, coef="LEW", adjust="BH",n=24000)
write.table(LEW_DEG, file="AMG_DEG_LEW.txt", quote=FALSE, row.names=TRUE)

#Check out results from Comparison4
topTable(fit2, coef="BNSS", adjust="BH")
BNSS_DEG <- topTable(fit2, coef="BNSS", adjust="BH",n=24000)
write.table(BNSS_DEG, file="AMG_DEG_BNSS.txt", quote=FALSE, row.names=TRUE)



































library(limma)
library(EDASeq)
####READ IN DATA
x <-read.ilmn(files="smplProbes_amg.txt")
head(x$E)

####BACKGROUND CORRECT WITH (http://svitsrv25.epfl.ch/R-doc/library/limma/html/backgroundcorrect.html)
y <- backgroundCorrect(x, method="normexp",offset=5)
head(y$E)

###QUANTILE NORMALIZATION WITH (http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/normalizebetweenarrays.html)
y <- normalizeBetweenArrays(y, method="quantile")
head(y$E)

###ALTERNATIVELY, PERFORM BACKGROUND CORRECTION,QUANTILE NORMALIZATION AND LOG2 TRANSFORMATION WITH NEQC
z <- neqc(x, offset=5)
dim(z)

###Input Target file with labels, then assign colours
colors <- c("firebrick1","brown","steelblue4")
targets <- read.delim("Targets.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(targets)
Group <- factor((targets$Group))


###PLOT BOXPLOTS, HISTOGRAMS, RLE PLOTS
pdf("QC_plot.pdf")
par(mfrow=c(3,3))
boxplot(log2(x$E),range=0,ylab="log2 intensity", las=2, cex.axis=0.6, main="Raw Expression", col=colors[Group])
boxplot(log2(y$E),range=0,ylab="log2 intensity", las=2, cex.axis=0.6, main="Quantile Normalized",col=colors[Group])
boxplot(log2(z$E),range=0,ylab="log2 intensity", las=2, cex.axis=0.6, main="Neqc Normalized",col=colors[Group])
hist(log2(x$E))
hist(log2(y$E))
hist(log2(z$E))
plotRLE(log2(x$E),ylim=c(-0.15,0.15), cex.axis=0.5, las=2, outline=FALSE, ylab="RLE", main="Raw Expression",col=colors[Group])
legend("top", fill =colors, cex=1, legend = c("Control", "CRS", "PCRS"), horiz=T)
plotRLE(log2(y$E),ylim=c(-0.15,0.15), cex.axis=0.5, las=2, outline=FALSE, ylab="RLE", main="Quantile Normalized",col=colors[Group])
legend("top", fill =colors, cex=1, legend = c("Control", "CRS", "PCRS"), horiz=T)
plotRLE(log2(z$E),ylim=c(-0.15,0.15), cex.axis=0.5, las=2, outline=FALSE, ylab="RLE", main="Neqc Normalized",col=colors[Group])
legend("top", fill =colors, cex=1, legend = c("Control", "CRS", "PCRS"), horiz=T)
#dev.off()

###PLOT PCA, CLUSTERING AND QQPLOTS
#pdf("QC_plot2.pdf")
par(mfrow=c(3,3))
plotPCA(x$E,cex=0.5,main="Raw Expression",col=colors[Group])
plotPCA(y$E,cex=0.5,main="Quantile Expression",col=colors[Group])
plotPCA(z$E,cex=0.5,main="Neqc Expression",col=colors[Group])
plot(hclust(dist(t(x$E), method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7)
plot(hclust(dist(t(y$E), method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7)
plot(hclust(dist(t(z$E), method="euclidean"),method="average"), las=2, cex=0.5, cex.axis=0.7)
qqnorm(x$E, pch=20, las=1, cex.axis=0.8, col="grey30");qqline(x$E, lwd=1.5, lty=3)
qqnorm(y$E, pch=20, las=1, cex.axis=0.8, col="grey30");qqline(y$E, lwd=1.5, lty=3)
qqnorm(z$E, pch=20, las=1, cex.axis=0.8, col="grey30");qqline(z$E, lwd=1.5, lty=3)
dev.off()

###REMOVE LOWLY EXPRESSED GENES (NON-SPECIFIC FILTERING)
###We keep probes that are expressed in at least three arrays according to a detection p-values of 5%:
expressed <- rowSums(z$other$Detection < 0.1) >= 1
z <- z[expressed,]
dim(z)

write.table(z$E, "Quantile_Normalized_AMG.txt", sep="\t")


library(limma)
exprsFile <- file.path("AMG_Matrix_voom.txt") #load normalized matrix
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
dim(exprs)

Group <- factor(targets$Combined, levels=c("Control_BNSS", "Control_F344", "Control_LEW", "Control_WKY", "CRS_BNSS", "CRS_F344", "CRS_LEW", "CRS_WKY", "PCRS_BNSS", "PCRS_F344", "PCRS_LEW", "PCRS_WKY"))

#create design matrix
design <- model.matrix(~0+Group)
colnames(design) <- c("Control_BNSS", "Control_F344", "Control_LEW", "Control_WKY", "CRS_BNSS", "CRS_F344", "CRS_LEW", "CRS_WKY", "PCRS_BNSS", "PCRS_F344", "PCRS_LEW", "PCRS_WKY")
colnames(design)
fit <- lmFit(exprs, design)

#specify the comparisons you want to make
contrast.matrix <- makeContrasts(
Comparison1=Control_WKY-PCRS_WKY,
Comparison2=Control_F344-PCRS_F344,
Comparison3=Control_LEW-PCRS_LEW,
Comparison4=Control_BNSS-PCRS_BNSS,
levels=design)

#use linear models to make those comparisons
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Check out results from Comparison1
topTable(fit2, coef="WKY", adjust="BH")
WKY_DEG <- topTable(fit2, coef="WKY", adjust="BH",n=24000)
write.table(WKY_DEG, file="AMG_DEG_WKY.txt", quote=FALSE, row.names=TRUE)

#Check out results from Comparison2
topTable(fit2, coef="F344", adjust="BH")
F344_DEG <- topTable(fit2, coef="F344", adjust="BH",n=24000)
write.table(F344_DEG, file="AMG_DEG_F344.txt", quote=FALSE, row.names=TRUE)

#Check out results from Comparison3
topTable(fit2, coef="LEW", adjust="BH")
LEW_DEG <- topTable(fit2, coef="LEW", adjust="BH",n=24000)
write.table(LEW_DEG, file="AMG_DEG_LEW.txt", quote=FALSE, row.names=TRUE)

#Check out results from Comparison4
topTable(fit2, coef="BNSS", adjust="BH")
BNSS_DEG <- topTable(fit2, coef="BNSS", adjust="BH",n=24000)
write.table(BNSS_DEG, file="AMG_DEG_BNSS.txt", quote=FALSE, row.names=TRUE)














mydata = read.table("Test.txt", sep="\t")
#par(mfrow=c(2,2))
colors=c("indianred3","steelblue")
#colors=c("red","blue")

plot <-cbind(mydata$V2,mydata$V3,mydata$V4,mydata$V5)
names=c("WKY", "F344", "BN-SS", "LEW")

x=barplot(plot,col=colors,las=2,cex.axis=0.6, cex=0.6,names.arg=names, angle=45,ylim=c(0,400), main="CRS", beside=T)
text(x=x, y=mydata$V2, labels=(mydata$V2), pos=3, cex=0.9, font=1, las=2)
legend("topleft", fill =c("indianred3","steelblue"), cex=0.6, legend =mydata$V1, bty='n')
abline(h=905, lty=2, col="grey")
abline(h=442, lty=2, col="steelblue4")


matrix <- read.delim("Test.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
par(mfrow=c(4,2))
color =matrix[,3]
color1 =matrix1[,7]
color2 =matrix3[,7]

 plot(matrix[2:1], pch=20, main="", cex.axis=0.5, cex=0.8,col=color,las=1, xant=F)
 abline(v=0.5)
 abline(v=-0.5)
 abline(h=1.30)


mydata = read.table("Test.txt", sep="\t")
par(mar=c(8,18,6,3))
x <-barplot(mydata$V2, las=2,horiz=T,xaxt="n", xlab="-log P",
border="black",col="grey", space=0.2,cex.axis=0.7, cex.names=0.7,names.arg=mydata$V1, las=1)
axis(1,at=c(0,1.0,2.0),cex.axis=0.7,lwd = 1.3,col="gray25", las=1, font=1) # draw y axis with required labels
abline(v=1.3, lty=3, col="indianred3")

text(x=x, y=mydata$V2, labels=(mydata$V2), pos=3, cex=0.7, font=1, las=2)


x <-barplot(mydata$V2, las=2,horiz=T,xaxt="n",
border="black",col="grey", space=0.3,cex.axis=0.6, cex.names=0.6,names.arg=mydata$V1)
axis(1,at=c(0,1.0,2.05),cex.axis=0.7,lwd = 1.3,col="gray25", las=2, font=1) # draw y axis with required labels
text(x=x, y=mydata$V2, labels=(mydata$V2), pos=3, cex=0.7, font=1, las=2)
text(x=x+0.5, y=-3,labels=(mydata$V1), xpd=TRUE, srt=60, pos=2, family="sans", cex=0.8)

text(x=x, y=mydata$V3, labels=(mydata$V3), pos=3, cex=0.7, font=2)


axis(1,at=x,labels=mydata$V1,las=2,cex.axis=0.7)

axis(2,seq(0,15,3),c(0,3,6,9,12,15))

colors1=c("gray94","azure3","azure4")



faov





mydata = read.table("Test.txt", sep="\t")
par(mar=c(8,4,4,4))
colors <- c("grey80","grey80","grey80", "grey40", "grey40", "grey40")
#plot <-cbind(mydata$V2,mydata$V3,mydata$V4,mydata$V5,mydata$V6,mydata$V7,mydata$V8,mydata$V9,mydata$V10,mydata$V11,mydata$V12,mydata$V13,mydata$V14,mydata$V15,mydata$V16,mydata$V17)
plot <-cbind(mydata$V2,mydata$V3,mydata$V4,mydata$V5,mydata$V6,mydata$V7)
#names=c("Pair1_CTRL_c1", "Pair1_CTRL_c2", "Pair1_CTRL_c3", "Pair1_PMS_c1", "Pair1_PMS_c2", "Pair2_CTRL_c1", "Pair2_CTRL_c2", "Pair2_PMS_c1", "Pair2_PMS_c2", "Pair3_CTRL_c1", "Pair3_CTRL_c2", "Pair3_CTRL_c3", "Pair3_PMS_c1", "Pair3_PMS_c2", "Pair3_PMS_c3", "Pair3_PMS_c4")
names=c("Cntl 1", "Cntl 2", "Cntl 3", "PMS 1 (DS:1)", "PMS 2 (DS: 85298)", "PMS 3 (DS: 6902568)")
x=barplot(plot,col=colors,las=2,cex.axis=0.6, cex=0.6,srt=45,names.arg=names, main=mydata$V1, beside=T,space = c(0,.2), ylim=c(0,10))
par(new=TRUE)
plot(mydata$V1, mydata$V3,type="o",col="indianred3",xaxt="n",yaxt="n",xlab="",ylab="", frame.plot=F,ylim=c(3,7))
axis(4, ylim=range(mydata$V3), col="black",col.axis="indianred3",las=1,cex.axis=0.6)

par(mfrow=c(4,3))
x=barplot(plot,col=colors[Group],las=2,cex.axis=0.6, cex=0.6,names.arg=names, main=mydata$V1, beside=T,space = c(0,.2))


text(x=x, y=mydata$V2, labels=(mydata$V2), pos=3, cex=0.9, font=1, las=2)
legend("topright", fill =c("grey80","grey40"), cex=0.6, legend =c("Cntl","PMS"), bty='n', horiz=F)
abline(h=905, lty=2, col="grey")
abline(h=442, lty=2, col="steelblue4")




rawdata <- read.delim("Test.txt", check.names=FALSE, stringsAsFactors=FALSE, row.names=1)
head(rawdata)
dim(rawdata)

par(mfrow=c(2,2))
verboseScatterplot(rawdata$one,rawdata$two, xlab="deletion size", ylab="SHANK3 expression", main="",cex = 0.5, cex.axis = 0.5, cex.lab = 0.5, cex.main = 1, abline=TRUE, abline.color=2,abline.lty=2)
