rm(list = ls(all=T))
graphics.off()

#General Bioconductor packages
library(Biobase)
library(oligoClasses)

#Annotation and data import packages
library(ArrayExpress)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)

#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)

#Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)

#Formatting/documentation packages
#library(rmarkdown)
#library(BiocStyle)
library(dplyr)
library(tidyr)

#Helpers:
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
#library(devtools)
setwd("D:/GitHub/D_FBA/Data/Gerosa/RawData/ArrayExpress/")
raw_data_dir = "D:/GitHub/D_FBA/Data/Gerosa/RawData/ArrayExpress/"
AEset = ArrayExpress("E-MTAB-3392")
phenodata = AEset@phenoData@data
featuredata = AEset@featureData@data

for (i in 1:8) {
  for (j in 1:8) {
    if (i!=j) {
      setwd("D:/GitHub/D_FBA/Data/Gerosa/RawData/ArrayExpress/")
      #anno_AE <- getAE("E-MTAB-3392", path = raw_data_dir, type = "raw")
      SDRF <- read.delim("E-MTAB-3392.sdrf.txt",check.names=FALSE,stringsAsFactors=FALSE)
      source = SDRF[,"Factor Value[growth condition]"]
      refcond = unique(source)[[i]]
      expcond = unique(source)[[j]]
      ref = which(source==refcond)
      exp = which(source==expcond)
      SDRF = SDRF[union(ref, exp),]
      source = SDRF[,"Factor Value[growth condition]"]
      Treatment <- factor(source, levels=c(refcond, expcond))
      
      x <- read.maimages(SDRF[,"Array Data File"],source="agilent", green.only=TRUE, other.columns="gIsWellAboveBG")
      y <- backgroundCorrect(x, method="normexp")
      y <- normalizeBetweenArrays(y, method="quantile")
      Control <- y$genes$ControlType==1L
      NegativeControl = y$genes$GeneName=="NegativeControl"
      Noname <- is.na(y$genes$GeneName)
      IsExpr <- rowSums(y$other$gIsWellAboveBG > 0) >= 4
      
      yfilt <- y[!Control & !NegativeControl & !Noname & IsExpr, ]
      dim(yfilt)
      
      design <- model.matrix(~ Treatment + 0)
      colnames(design) <- levels(Treatment)
      fit2 <- lmFit(yfilt,design)
      cont.matrix <- makeContrasts(eval(parse(text=paste(expcond,"-", refcond, sep = "", collapse = ""))), levels=design)
      fit2 <- contrasts.fit(fit2, cont.matrix)
      fit2 <- eBayes(fit2)
      tT2 = topTable(fit2,adjust="fdr", sort.by="B", number=10000)
      tT2 <- subset(tT2, select=c("GeneName","SystematicName","adj.P.Val","P.Value","t","B","logFC"))
      tT2 = tT2[tT2$adj.P.Val<0.05,]
      
      setwd("D:/GitHub/D_FBA/Data/Gerosa/RawData/")
      write.xlsx2(tT2, paste(expcond,"-", refcond,".xlsx", sep = "", collapse = ""), row.names=F)
    }
    
  }
  
}
