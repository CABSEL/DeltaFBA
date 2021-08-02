rm(list = ls(all=T))
graphics.off()
library(here)
library(tidyverse)
library(viridis)
library(stringr)
library(fgsea)
library(plotly)
library(heatmaply)
library(R.matlab)
library(xlsx)
library(tdr) #tdStats
library(Metrics)
library(ggthemes)
library(extrafont)
library(likert)
library(stringr)
library(reshape2)

setwd(paste0(here(),"/Data/Ishii/"))
wd = read.xlsx("Fluxomics.xlsx", 1)
conditions = colnames(wd)[-c(1,2)]
conditions[1] = "D = 0.1/h"
conditions[2] = "D = 0.4/h"
conditions[3] = "D = 0.5/h"
conditions[4] = "D = 0.7/h"
fluxes = as.character(wd[,1])
wtflux = wd[,2]
M = wd[,-c(1,2)]

setwd(paste0(here(),"/Simulations/Ishii/RT_PCR_Relaxed_DFBA/"))
Md = read.xlsx("Results.xlsx",1, header = F)
Pd = read.xlsx("Results.xlsx",2, header = F)
dat  = matrix(0, nrow = ncol(Md), ncol = nrow(Md))
for (i in 1:length(fluxes)) {
  exp = as.numeric(Md[i,])
  sim = as.numeric(Pd[i,])
  rmse = c()
  for (j in 1:length(sim)) {
    e = tdStats(sim[j], exp[j], functions=c("rmse"))
    rmse[j] = e/diff(range(as.numeric(Md[,j])))
  }
  dat[,i] = rmse
}
colnames(dat) = fluxes

dat = melt(dat)

rm(list= ls()[!(ls() %in% c('dat'))])


setwd(paste0(here(),"/Data/Ishii/"))
wd = read.xlsx("Fluxomics.xlsx", 1)
conditions = colnames(wd)[-c(1,2)]
conditions[1] = "D = 0.1/h"
conditions[2] = "D = 0.4/h"
conditions[3] = "D = 0.5/h"
conditions[4] = "D = 0.7/h"
fluxes = as.character(wd[,1])
wtflux = wd[,2]
M = wd[,-c(1,2)]

setwd(paste0(here(),"/Simulations/Ishii/RT_PCR_Relaxed_DFBA/"))
Md = read.xlsx("Results.xlsx",1, header = F)
Pd = read.xlsx("Results.xlsx",2, header = F)

m1 = c()

for (i in 1:length(conditions)) {
  sim = Pd[,i]
  exp = Md[,i]
  
  m1[i] = tryCatch(tdStats(sim, exp, functions= c("nrmse")) , error  = function(e) NA)
  
}
adddat = data.frame(Var1 = 1:length(m1), Var2 = "dFBA", value = m1)
dat = rbind(adddat,dat)

cols =  c("#4379a7",rep("#f8766d",length(fluxes)))
dat$Var2 = as.character(dat$Var2)
dat$Var2[which(dat$Var2=="dFBA")] = paste0("\u0394", "FBA")
dat$Var2 = factor(dat$Var2, levels = unique(dat$Var2))
  

setwd(paste0(here(),"/Plotting/Ishii/"))

tiff('IndividualMetrics.tiff', units="px", width=(1400*3), height=(700*3), res=300)
p <- ggplot(dat, aes(x = Var2, y = value, fill= Var2)) + 
  geom_boxplot(outlier.size = 2) + 
  stat_boxplot(geom = "errorbar",width = 0.4,size = 0.75) + 
  geom_boxplot(lwd=0.55) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = cols) +
  theme(legend.position = "none") +
  #ylim(-2,2) +
  
  xlab("Fluxes") + ylab("") +
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold",angle = 90, hjust = 1, vjust = 0.4))

p
dev.off()
