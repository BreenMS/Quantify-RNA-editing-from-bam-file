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
targets <- read.delim("MetaData.txt", check.names=FALSE, stringsAsFactors=FALSE)
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

## MEDIAN SUMMARIZE PROBES MATACHING TO THE SAME GENE SYMBOL
##Begin median summary

write.csv(y, "Norm_expression.csv",row.names=FALSE)
read.csv("Norm_expression.csv",header=TRUE)->mat_med
dim(mat_med)

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

##OUTLIER ANALYSIS by 3D PCA
library(rgl)
matrix <-  read.table("Norm_median_clean.txt",header=TRUE, row.names=1, sep="\t")
head(matrix)
dim(matrix)

#CREATE 3D INTERACTIVE PCA PLOT FOR OUTLIER IDENTIFICATION
tmatrix<-t(matrix)
pcs<-prcomp(tmatrix)
summary(pcs)
PC1<-pcs$x[,1]
PC2<-pcs$x[,2]
PC3<-pcs$x[,3]

PCA_details <- cbind(PC1, PC2, PC3)

#Create 3D PCA
plot3d(PC1, PC2, PC3, main="Condition", xlab="PC1", ylab="PC2", zlab="PC3",
       ,box=FALSE, size=1, type="s")

labels<-rownames(tmatrix)
text3d(PC1+6, PC2+6, PC3+6,text=labels, font=1, cex=1)
 
#ADD ELLIPSOIDS AS STANDARD DEVIATIONS
mean.vec <- c(mean(PC1), mean(PC2), mean(PC3))
allcomp <- cbind(PC1, PC2, PC3)
sigma <- cov(allcomp) #sigma <-  c(sd(PC1), sd(PC2), sd(PC3))
plot3d( ellipse3d(x = sigma, centre=mean.vec,scale = c(1.5, 1.5, 1.5),), col="PeachPuff", alpha=0.50, add = TRUE, level=0.95, smooth=TRUE) #1.5 SD
plot3d( ellipse3d(x = sigma, centre=mean.vec,scale = c(2, 2, 2),), col="GhostWhite", alpha=0.50, add = TRUE, level=0.95,smooth=TRUE) #2 SD

#Create 2D PCA
colors=c("steelblue4", "steelblue3", "aquamarine4","indianred4","salmon")
plot(PC1, PC2, col=colors[Group], pch=16, cex=0.8, las=1, cex.axis=0.7, xlab="PC1: 16% variance", ylab="PC2: 12% variance")
abline(h=0, lty=2, lwd=0.5)
abline(v=0, lty=2, lwd=0.5)
legend("bottomleft", fill =colors, cex=0.6, legend = c("Control", "TE", "Dep", "PTSD","PTSD/Dep"), bty="n")

#UNSUPERVISED CLUSTERING WITH PEARSONS CORRELATION, FOR OUTLIER IDENTIFICATION
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
traitData <- read.delim("MetaData.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")
dim(traitData)
head(traitData)
rowsExpr <- rownames(dataExpr0)
traitRows <- match(rowsExpr,traitData$SID.1)
datTraits = traitData[traitRows, -1];
rownames (datTraits) = traitData[traitRows, 1];
table(rownames(datTraits)==rownames(dataExpr0)) # EVERYTHING OK?

A=adjacency(t(dataExpr0))
k=as.numeric(apply(A,2,sum))-1
Z.k=scale(k)

#DESIGNATE SAMPLES AS OUTLIERS WITH RED COLORS
thresholdZ.k= -2 #2 SD from average
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


###QUALITY CONTROL PLOTS - variancePartion
library(variancePartition)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

exprsFile <- file.path("Norm_median_clean.txt") #load normalized matrix
exprs <- as.matrix(read.table(exprsFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE))
dim(exprs)

info <- read.delim("MetaData", check.names=FALSE, stringsAsFactors=FALSE)
head(info)

form <- ~ RIN + (1|TE) + (1|Alcohol) + (1|Nicotine) + (1|Ethnicity) + (1|Delivery) + (1|Gender) + (1|PrevMed) + (1|PregMed) + (1|Batch)  + (1|Depression) + (1|PTSD) 

varPart <- fitExtractVarPartModel(exprs, form, info )
vp <- sortCols( varPart )
plotVarPart(vp,label.angle=50)

#####DIFFERENTIAL EXPRESSION ANALYSIS WITH LIMMA ######
library(limma)
exprs <- read.delim("Matrix.txt", check.names=FALSE, stringsAsFactors=FALSE, header=T, row.names=1)
head(exprs)
dim(exprs)

targets <- read.delim("MetaData.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(targets)
Group <- factor(targets$PTSD, levels=c("Case","Control"))
Group
design <- model.matrix(~0+Group)
colnames(design) <- c("Case","Control")
colnames(design)

#Add covariates
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
write.table(Diff, file="DGE_result.txt", quote=FALSE, row.names=TRUE)

#create volcano plots of DGE signaturs
plot(Diff$logFC, -log10(Diff$P.value))

###QUALITY CONTROL PLOTS of DGE PC1 association with meta-data
matrix <-  read.table("Norm_median_clean.txt",header=TRUE, row.names=1, sep="\t") # a matrix of DGE signatures only
nSamples = ncol(matrix)
tmatrix<-t(matrix)
pcs<-prcomp(tmatrix)
summary(pcs)
PC1<-pcs$x[,1]

#Load meta-data
targets <- read.delim("MetaData.txt", stringsAsFactors=FALSE, row.names=1, header=TRUE, sep="\t")

moduleTraitCor = cor(PC1, targets, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples) #requries library(WGCNA)
targets <- read.delim("testing.txt", stringsAsFactors=FALSE, row.names=1, header=F, sep="\t")
barplot((-log10(targets$V2)), horiz=T, las=1, col="white", cex.axis=0.7, cex=0.7, names.arg=rownames(targets))
abline(v=-log10(0.05), col="red")
