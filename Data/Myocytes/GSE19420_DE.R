rm(list = ls(all=T))
graphics.off()
library(here)
library(Biobase)
library(GEOquery)
library(limma)
setwd(paste0(here(),"Data/Myocytes/"))

gset <- getGEO("GSE19420", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL571", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

phenodata = gset@phenoData@data
featuredata = gset@featureData@data
unique(phenodata$characteristics_ch1)
ngt = which(phenodata$characteristics_ch1==unique(phenodata$characteristics_ch1)[2])
t2d = which(phenodata$characteristics_ch1==unique(phenodata$characteristics_ch1)[3])

gsms = rep("X", times = length(phenodata$title))
gsms[ngt] = "0"
gsms[t2d] = "1"
gsms = paste(gsms, collapse = "")

sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=10000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
tT = tT[tT$P.Value<0.05,]

library(biomaRt)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)


ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
searchFilters(mart = ensembl, pattern = "ensembl")
df = getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), 
           filters = 'hgnc_symbol', 
           values = as.character(tT$Gene.symbol), 
           mart = ensembl) 

tT$Ensembl = NA

for (i in 1:length(tT$ID)) {
  ind = which(df$hgnc_symbol==tT$Gene.symbol[i])
  if (length(ind)>0) {
    ind = ind[1]
    tT$Ensembl[i] = df$ensembl_gene_id[ind]
  }
}
tT = tT[complete.cases(tT),]
write.table(tT,"GSE19420_DE.txt" , row.names=F, sep="\t")

