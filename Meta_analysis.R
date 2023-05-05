##########meta analysis gastric##############################

library(pca3d)
library(RColorBrewer)
library(affy)
library(AnnotationHub)
library(base)
library(bindrcpp)
library(Biobase)
library(BiocUpgrade)
library(Biostrings)
library(bit64)
library(bitops)
library(data.table)
library(DelayedArray)
library(DESeq2)
library(devtools)
library(digest)
library(dplyr)
library(EDASeq)
library(edgeR)
library(EnsDb.Hsapiens.v75)
library(ensembldb)
library(exprAnalysis)
library(extrafont)
library(factoextra)
library(genefilter)
library(GenomeInfoDb)
library(GenomeInfoDbData)
library(GenomicFeatures)
library(GenomicRanges)
library(GEOquery)
library(ggbiplot)
library(ggfortify)
library(ggplot2)
library(ggrepel) 
library(ggthemes)
library(gplots)
library(gProfileR)
library(grid)
library(gridExtra)
library(grl)
library(igraph)
library(imputeTS)
library(IRnges)
library(kableExtra)
library(knitr)
library(lattice)
library(limma)
library(LSD)
library(matrixStats)
library(metaRNASeq)
library(pca3d)
library(pcaGoPromoter)
library(pheatmap)
library(plyr)
library(preprocessCore)
library(RColorBrewer)
library(readr)
library(readxl)
library(reshape)
library(reshape2)
library(rlang)
library(Rsamtools)
library(RSQLite)
library(rtracklayer)
library(samr)
library(scater)
library(scatterplot3d)
library(shiny)
library(shinyWidgets)
library(statmod)
library(stats)
library(STRINGdb)
library(stringr)
library(SummarizedExperiment)
library(sva)
library(tidyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tximport)
library(tximportData)
library(utils)
library(VennDiagram)
library(vsn)
library(wasabi)
library(XML)
require(cowplot)
require(dplyr)
require(factoextra)
require(FactoMineR)
require(GenomicAlignments)
require(ggplot2)
require(ggplot2,ggthemes)
require(MASS)
require(reshape2)
require(Rsamtools)
require(tidyr)
require(VennDiagram)



########difrrential expression analysis and batch effect removal######## 


###difrrential expression analysis and batch effect removal 

dds< DESeqDataSetFromMatrix(countData, colData, formula(~ batch +  ~ condition))
dds<-DESeqDataSetFromMatrix(countData=countTable3,colData=coldata,design = ~batch+cond1*cond2)

dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ batch + treatment)
dds$treatment <- factor(dds$treatment, levels=c("control", "BD"))
dds$batch <- factor(dds$batch, levels=c("A", "B"))
dds <- DESeq(dds, full=design(dds), reduced = ~ batch)

dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ batch + condition) #asli 

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ batch + condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_Tumor_vs_normal")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="condition_Tumor_vs_normal", type="apeglm")

class(res)
### creat group of samples

setwd("C:/Users/susi/Desktop/GAST_NEW/newmerge")
adogr = as.matrix(read.table("phenoTable.txt"))
RAWmerge = as.matrix(read.table("mergebig1.txt"))
dim(RAWmerge)
RAWmerge
colData = as.matrix(read.table("colData.txt"))
colData = as.data.frame(colData)
write.csv(colData, file="colData.csv", row.names=T) 

batch1 <- rep(1, ncol(SRP170025))
batch2 <- rep(2, ncol(SRP133891))
batch3 <- rep(3, ncol(SRP135952))
batch4 <- rep(4, ncol(SRP073446))
Batch.type <- as.factor(c(batch1,batch2,batch3,batch4)) 
Batch.type
class(Batch.type)
class(Batch.type)

condition = factor(c( rep("Tumor", 3), rep("normal", 3), rep("Tumor", 9), rep("normal", 12), rep("Tumor", 6), rep("normal", 6), 
                      rep("Tumor", 6), rep("normal", 3)))

head(condition)

### Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
colData <- data.frame(condition= condition, Batch=Batch.type)
colData
class(colData)
head(colData)
dds <- DESeqDataSetFromMatrix(countData = RAWmerge, colData = colData, design = ~ Batch +  condition)

class(RAWmerge)

#storage.mode(RAWmerge) = "integer"

###normalize with deseq2
library(DESeq2)
cds <- DESeq(cds)
resultsNames(cds)
library(DESeq2)
cds <- DESeq(cds)
res <- results(cds)
rld <- rlog(cds, blind=F)
vsd <- varianceStabilizingTransformation(cds, blind=TRUE)
cntado <- counts(cds, normalized=T)
cntadolog <- log2(1+counts(cds, normalized=T))
dir()
dif <- results(dds, pAdjustMethod = "BH", alpha = 0.001)
dif$padj <- p.adjust(dif$pvalue, method="BH")
dif <- dif[order(dif$padj),]
down = subset(dif, log2FoldChange < -2 & padj< 0.001)
up = subset(dif, log2FoldChange > 2 & padj< 0.001)
nrow(up)
nrow(down)
dif
dif = as.data.frame(dif)
up = as.data.frame(read.delim("up.txt"))
down = as.matrix(read.delim("down.txt"))
resSig05 = as.matrix(read.delim("resSig05.txt"))

write.table(up , "up.txt" , sep = "\t" , quote = F , row.names = F , col.names = T)
write.table(down1 , "down1.txt" , sep = "\t" , quote = F , row.names = T , col.names = T)
write.table(resSig05 , "resSig05.txt" , sep = "\t" , quote = F , row.names = T , col.names = T)
write.table(up1 , "up1.txt" , sep = "\t" , quote = F , row.names = F , col.names = F)


#####################start to analysis by taximport


stringsAsFactors = FALSE
setwd("C:/Users/susi/Desktop/GAST_NEW/tximerge")
dir = "C:/Users/susi/Desktop/GAST_NEW/tximerge"
samples <- read.table(file.path(dir, "salmon1.txt"), header = TRUE)
samples
class(samples)
samples = as.data.frame(samples)



files <- c("quant766.sf", "quant764.sf","quant762.sf","quant761.sf",
           "quant763.sf", "quant765.sf", "quant512.sf",  "quant514.sf", "quant529.sf",
           "quant534.sf", "quant536.sf", "quant544.sf", "quant550.sf" , "quant551.sf",
           "quant554.sf", "quant513.sf", "quant515.sf", "quant517.sf", "quant523.sf",
           "quant530.sf", "quant531.sf","quant533.sf", "quant535.sf", "quant543.sf" ,
           "quant549.sf","quant553.sf", "quant557.sf","quant674.sf",  "quant676.sf",
           "quant678.sf","quant680.sf","quant682.sf", "quant684.sf", "quant683.sf", 
           "quant681.sf", "quant679.sf", "quant677.sf","quant673.sf","quant675.sf",
           "quant113.sf", "quant112.sf","quant111.sf","quant110.sf", "quant109.sf",
           "quant108.sf", "quant107.sf", "quant106.sf", "quant105.sf")  






names(files) <- paste0("samp", 1:48)
all(file.exists(files))




#2) Now we can import the salmon quantification.
#create a tx2gene.txt table

#create a tx2gene.txt table
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
Tx.ensemble <- transcripts(edb,columns = c("tx_id", "gene_id", "gene_name"), return.type = "DataFrame")
nrow(Tx.ensemble)
tx2gene<- Tx.ensemble[,c(1,2)]
head(tx2gene)

tx2gene = as.matrix(read.delim("tx2gene.txt"))
class(tx2gene)
tx2gene = as.data.frame(tx2gene)
tx2gene
write.csv(tx2gene, file="tx2gene.csv", row.names=T)

#3) creat txi 

txi <- tximport(files, type="salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
names(txi)
head(txi$counts)
dim(txi$counts)
txi$counts
class(txi)
tximerge = txi$counts
dim(tximerge)
tximerge <- tximerge[!rowSums(tximerge==0)>=1, ]
tximerge<-as.data.frame(tximerge)
dim(tximerge)
#SRP073446  <- SRP073446[rowMeans(SRP073446) > 50, ]

#remove NA from matrix 
any(is.na(tximerge))
tximerge <- na.omit(tximerge)
tximerge[is.na(tximerge)] <- 0
class(tximerge)
dim(tximerge)
tximerge = as.data.frame(tximerge)
head(tximerge)




dim(tximergelog)

tximergelog = log2(tximerge + 1)

tximergelog <- tximergelog[!rowSums(tximergelog==0)>=1, ]



setwd("C:/Users/susi/Desktop/GAST_NEW/tximerge")
write.csv(tximerge, file="tximerge.csv", row.names=T)
write.csv(tximergelog, file="tximergelog.csv", row.names=T)

SRP073446log = as.matrix(read.delim("SRP073446log.txt"))
SRP073446 = as.matrix(read.delim("SRP073446.txt"))
#////////////////////////////////////////////////////////
dds <- makeExampleDESeqDataSet(betaSD=1,interceptMean=10)
dds$batch <- factor(rep(c("A","B"),each=6))
#VST, remove batch effect, then plotPCA:

vsd <- vst(dds)
plotPCA(vsd, "batch")
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
plotPCA(vsd, "batch")


((((((((#before meta analysis))))))))))))))))))
  
  
  #Assign experimental variables:
  #Deseq2
  #https://fishycat.netlify.com/en/2017/09/rnaseqmm/
  
  
  #creat design matrix
  batch1 <- rep(1, ncol(SRP170025))
  batch2 <- rep(2, ncol(SRP133891))
  batch3 <- rep(3, ncol(SRP135952))
  batch4 <- rep(4, ncol(SRP073446))
  Batch.type <- as.factor(c(batch1,batch2,batch3,batch4)) 
  Batch.type
  class(Batch.type)
  class(Batch.type)
  
  condition = factor(c( rep("Tumor", 3), rep("normal", 3), rep("Tumor", 9), rep("normal", 12), rep("Tumor", 6), rep("normal", 6), 
                        rep("Tumor", 6), rep("normal", 3)  ))
  
  colData <- data.frame(condition= condition, Batch=Batch.type)
  colData
  class(colData)
  head(colData)
  
  
  #run dseq2
  
  ddsTxi <- DESeqDataSetFromTximport(txi,  colData = colData, design = ~ Batch +  condition)
  dds <- DESeq(ddsTxi, betaPrior=FALSE)
  
  #Pre-filtering
  dds <- dds[ rowSums(counts(dds)) > 10, ]
  dim(dds)
  
  class(dds)
  #Get counts:
  
  counttxi = counts( dds, normalized=TRUE )
  dim(counttxi)
  
  #Filter out low expression transcripts:
  
  counttxi <- counttxi[!rowSums(counttxi==0)>=1, ]
  counttxi<-as.data.frame(counttxi)
  dim(counttxi)
  
  
  counttxilog <- log2(1+counttxi)
  dim(counttxilog)
  counttxilog<-counttxilog[!rowSums(counttxilog==0)>=1, ]
  counttxilog<-as.data.frame(counttxilog)
  dim(counttxilog)
  
  ########### Rename a column in R
  
  names(count17344)[1]<-"SRR3400107"
  names(count17344)[2]<-"SRR3400106"
  names(count17344)[3]<-"SRR3400105"
  names(count17344)[4]<-"SRR3400113"
  names(count17344)[5]<-"SRR3400112"
  names(count17344)[6]<-"SRR3400111"
  names(count17344)[7]<-"SRR3400110"
  names(count17344)[8]<-"SRR3400109"
  names(count17344)[9]<-"SRR3400108"
  
  
  ########### Rename a column in R
  names(count17344log)[1]<-"SRR3400107"
  names(count17344log)[2]<-"SRR3400106"
  names(count17344log)[3]<-"SRR3400105"
  names(count17344log)[4]<-"SRR3400113"
  names(count17344log)[5]<-"SRR3400112"
  names(count17344log)[6]<-"SRR3400111"
  names(count17344log)[7]<-"SRR3400110"
  names(count17344log)[8]<-"SRR3400109"
  names(count17344log)[9]<-"SRR3400108"
  #reasult
  setwd("C:/Users/susi/Desktop/gastricdif/quant173446")
  write.csv(count17344, file="count17344.csv", row.names=T)
  write.csv(counttxilog, file="counttxilog.csv", row.names=T)
  
  read.csv(quant173446, header = TRUE, sep = "\t", dec = "." )
  read.csv(quant173446log, header = TRUE, sep = "\t", dec = ".")
  
  quant173446 = as.matrix(read.delim("quant173446.txt"))
  quant173446log = as.matrix(read.delim("count17344log.txt"))
  
  
  
  
  # get Degs with dif 
  setwd("C:/Users/susi/Desktop/GAST_NEW/tximerge")
  
  allc1 = as.matrix(read.table("allc1.txt"))
  dim(allc1)
  
  
  
  dif <- results(res, c("condition", "normal", "tumor"))
  dif <- dif[complete.cases(dif),]  #remove any rows with NA
  dif = res[dif[, "baseMean"] != 0, ]  # remove zeros
  mcols(dif, use.names=TRUE)
  summary(dif)
  nrow(dif)
  head(dif)
  class(dif)
  dim(dif)
  
  #View results
  #dif <- results(dds, pAdjustMethod = "BH", alpha = 0.01)
  dif$padj <- p.adjust(dif$pvalue, method="BH")
  dif <- dif[order(dif$padj),]
  head(dif)
  class(dif)
  dim(dif)
  dif = as.data.frame(dif)
  deg <- subset(as.data.frame(dif), padj < 0.01) #alldegs
  down = subset(deg1, log2FoldChange < -1 & padj< 0.01)
  up = subset(deg1, log2FoldChange > 1 & padj< 0.01)
  nrow(up)
  nrow(down)
  nrow(deg1)
  head(up)
  head(down)
  head(deg)
  nrow(deg1)
  
  #resultsHow many adjusted p-values were less than 0.001?
  sum(dif$padj < 0.01, na.rm=TRUE)
  
  
  
  
  #convert ensemble id to gene name 
  
  #Add Gene Symbols to Results output and export
  library(dplyr)
  library( org.Hs.eg.db ) 
  library(AnnotationDbi) 
  deg$symbol <- mapIds(org.Hs.eg.db, keys=row.names(deg), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  deg <- deg[, c(7, 1, 2, 3, 4, 5, 6)]
  deg <- as.data.frame(deg)
  head(deg)
  rownames(deg)= NULL
  deg1 <- as.matrix(deg[,-1])
  head(deg1)
  rownames(deg1) <- deg[,1]
  head(deg1)
  dim(deg1)
  deg1
  #remove NA from matrix 
  any(is.na(deg1))
  deg1 <- na.omit(deg1)
  deg1[is.na(deg1)] <- 0
  class(deg1)
  dim(deg1)
  deg1 = as.data.frame(deg1)
  deg = row.names(deg1)
  
  write.csv(deg1, file="deg1.csv", row.names=T)
  write.csv(deg, file="deg.csv", row.names=T)
  
  
  up1 = row.names(up)
  down1 = row.names(down)
  write.csv(result, file="result.csv", row.names=T)
  write.csv(deg, file="deg.csv", row.names=T)
  write.csv(up1, file="up1.csv", row.names=T)
  write.csv(up, file="up.csv", row.names=T)
  write.csv(down, file="down.csv", row.names=T)
  write.csv(down1, file="down1.csv", row.names=T)
  write.csv(dif, file="dif.csv", row.names=T)
  library(dplyr)
  library( org.Hs.eg.db ) 
  library(AnnotationDbi) 
  dif$symbol <- mapIds(org.Hs.eg.db, keys=row.names(dif), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  head(dif)
  
  write.csv(allc, file="allc.csv", row.names=T)
  
  
  #########################metaseq
  library("metaSeq")
  library("snow")
  condition = factor(c( rep("Tumor", 3), rep("normal", 3), rep("Tumor", 9), rep("normal", 12), rep("Tumor", 6), rep("normal", 6), 
                        rep("Tumor", 6), rep("normal", 3)  ))
  
  flag1 <- c(1,1,1,0,0,0, 1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0, 0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0)
  flag2 <- c("A","A","A","A","A","A","B","B", "B","B","B","B", "B","B","B","B", "B","B","B","B", "B","B","B","B","B","B","B",
             "C","C","C","C","C","C","C","C","C","C","C","C", "D","D","D","D","D","D","D","D","D")
  
  dim(tximerge)
  
  cds <- meta.readData(data = tximerge, factor = flag1, studies = flag2)
  cl <- makeCluster(4, "SOCK")
  
  result <- meta.oneside.noiseq(cds, k = 0.5, norm = "tmm", replicates = "biological",
                                factor = flag1, conditions = c(1, 0), studies = flag2, cl = cl)
  
  data(Result.Meta)
  
  result <- Result.Meta
  class(result)
  
  head(F$Upper)
  head(F$Lower)
  
  F = (F$Upper)
  
  
  
  ############metaseq
  pc = prcomp(allc)
  pcr = data.frame(pc$r[,1: 3], Batch.type)
  
  Dataset = c(rep("SRP170025", 6), rep("SRP133891", 21), rep("SRP135952", 12), rep("SRP073446", 9))
  
  condition = factor(c( rep("Tumor", 3), rep("normal", 3), rep("Tumor", 9), rep("normal", 12), rep("Tumor", 6), rep("normal", 6), 
                        rep("Tumor", 6), rep("normal", 3)  ))
  
  ggplot(pcr, aes(PC1, PC2, color=Dataset, shape=Batch.type))+geom_point(size=3)+theme_bw()
  
  