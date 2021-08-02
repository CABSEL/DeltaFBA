rm(list = ls(all=T))
graphics.off()
library(here)

options(java.parameters = "-Xmx9000m")
library(data.table)
library(tidyverse)
library(xlsx)
library(tdr) #tdStats
library(Metrics)
library(ggthemes)
library(extrafont)
library(likert)

setwd(paste0(here(),"/Data/Ishii/"))
wd = read.xlsx("Fluxomics.xlsx", 1)
conditions = colnames(wd)[c(9,10,15,21,24,5,6)]
conditions[6] = "D = 0.5/h"
conditions[7] = "D = 0.7/h"
fluxes = as.character(wd[,1])
wtflux = wd[,2]
M = wd[,c(9,10,15,21,24,5,6)]

setwd(paste0(here(),"/Simulations/Ishii/Array_Relaxed_DFBA/"))
Md = read.xlsx("Results.xlsx",1, header = F)
Pd = read.xlsx("Results.xlsx",2, header = F)
Rd = read.xlsx("Results.xlsx",3, header = F)
MRatio = M/wtflux
MRatio = as.matrix(MRatio)
MRatio[is.nan(MRatio)] = 0
MRatio[which(MRatio==Inf)] = 100
MRatio[which(MRatio==-Inf)] = 0.001


m1 = c()
m2 = c()
m3 = c()
m4 = c()
m5 = c()
m6 = c()

for (i in 1:length(conditions)) {
  sim = Pd[,i]
  exp = Md[,i]
  
  m1[i] = tryCatch(tdStats(sim, exp, functions= c("nrmse")) , error  = function(e) NA)
  
  m2[i] =  tryCatch(sim %*%  exp / (norm(as.matrix(sim), type = "2") * norm(as.matrix(exp), type = "2")), error = function(e) NA)
  
  s_Md = as.numeric(exp)
  sel = which(abs(s_Md)>median(abs(s_Md)))
  sel2 = which(abs(s_Md)>quantile(abs(s_Md),0.75))
  s_Md[s_Md>0] = 1;
  s_Md[s_Md<0] = -1;
  
  s_Pd = as.numeric(sim)
  s_Pd[s_Pd>0] = 1;
  s_Pd[s_Pd<0] = -1;
  
  m3[i] = accuracy(s_Md,s_Pd)
  m4[i] = accuracy(s_Md[sel], s_Pd[sel])
  m5[i] = accuracy(s_Md[sel2], s_Pd[sel2])
  
  m6[i] = summary(lm(MRatio[,i]~Rd[,i]))$r.squared
  
}


data <- data.frame(
  condition=rep(conditions,3),
  group=c(rep('NRMSE', length(conditions)), rep("\u03c1", length(conditions)),rep('Sign Acc', length(conditions))) ,
  value= c(m1, m2,m3))

data$test = "Global"

rm(list= ls()[!(ls() %in% c('data'))])
graphics.off()
options(java.parameters = "-Xmx9000m")
library(data.table)
library(tidyverse)
library(xlsx)
library(tdr) #tdStats
library(Metrics)
library(ggthemes)
library(extrafont)
library(likert)

setwd(paste0(here(),"/Data/Ishii/"))
wd = read.xlsx("Fluxomics.xlsx", 1)
conditions = colnames(wd)[c(9,10,15,21,24,5,6)]
conditions[6] = "D = 0.5/h"
conditions[7] = "D = 0.7/h"
fluxes = as.character(wd[,1])
wtflux = wd[,2]
M = wd[,c(9,10,15,21,24,5,6)]

setwd(paste0(here(),"/Simulations/Ishii/RT_PCR_Relaxed_DFBA//"))
Md = read.xlsx("Results.xlsx",1, header = F)
Pd = read.xlsx("Results.xlsx",2, header = F)
Rd = read.xlsx("Results.xlsx",3, header = F)
Md = Md[,c(7,8,13,19,22,3,4)]
Pd = Pd[,c(7,8,13,19,22,3,4)]
Rd = Rd[,c(7,8,13,19,22,3,4)]
MRatio = M/wtflux
MRatio = as.matrix(MRatio)
MRatio[is.nan(MRatio)] = 0
MRatio[which(MRatio==Inf)] = 100
MRatio[which(MRatio==-Inf)] = 0.001

m1 = c()
m2 = c()
m3 = c()
m4 = c()
m5 = c()
m6 = c()

for (i in 1:length(conditions)) {
  sim = Pd[,i]
  exp = Md[,i]
  
  m1[i] = tryCatch(tdStats(sim, exp, functions= c("nrmse")) , error  = function(e) NA)
  
  m2[i] =  tryCatch(sim %*%  exp / (norm(as.matrix(sim), type = "2") * norm(as.matrix(exp), type = "2")), error = function(e) NA)
  
  s_Md = as.numeric(exp)
  sel = which(abs(s_Md)>median(abs(s_Md)))
  sel2 = which(abs(s_Md)>quantile(abs(s_Md),0.75))
  s_Md[s_Md>0] = 1;
  s_Md[s_Md<0] = -1;
  
  s_Pd = as.numeric(sim)
  s_Pd[s_Pd>0] = 1;
  s_Pd[s_Pd<0] = -1;
  
  m3[i] = accuracy(s_Md,s_Pd)
  m4[i] = accuracy(s_Md[sel], s_Pd[sel])
  m5[i] = accuracy(s_Md[sel2], s_Pd[sel2])
  
  m6[i] = summary(lm(MRatio[,i]~Rd[,i]))$r.squared
  
}


data2 <- data.frame(
  condition=rep(conditions,3),
  group=c(rep('NRMSE', length(conditions)), rep("\u03c1", length(conditions)),rep('Sign Acc', length(conditions))) ,
  value= c(m1, m2,m3))
data2$test = "Local"
rm(list= ls()[!(ls() %in% c('data', 'data2', 'p_cond'))])


dat = rbind(data,data2)
dat$test = factor(dat$test, levels = c("Global", "Local"))
dat$group = factor(dat$group, levels = c("NRMSE","\u03c1", "Sign Acc"))

cols = c("#2d506f", "#4983b5")

setwd(paste0(here(),"/Plotting/Ishii//"))

tiff('GlobalvLocal.tiff', units="px", width=(700*3), height=(450*3), res=300)

p <- ggplot(dat, aes(x = test, y = value, fill= test)) + 
  geom_boxplot(outlier.size = 2) + 
  stat_boxplot(geom = "errorbar",width = 0.4,size = 0.75) + 
  geom_boxplot(lwd=0.55) + 
  theme_bw(base_size = 12) +
  facet_grid(.~group)+
  theme(strip.text.x = element_text(angle = 0, family = "Calibri", face = "bold", size = 16))+
  theme(panel.spacing =unit(0.5, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, size = 1.5), 
        strip.background = element_rect(color = "#476b6b", size = 1.5, fill = "#d6dce4"))+
  #theme(strip.text.x = element_blank())+
  #theme(panel.spacing =unit(0.5, "lines")) +
  scale_fill_manual(values = cols) +
  theme(legend.position = "none") +
  #ylim(-1,1) +
  xlab("") + ylab("") + 
  ylim(-0.5,1) +
  theme(axis.ticks.y = element_line(size=1,color='#476b6b')) +
  theme(axis.ticks.x = element_blank()) + 
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_blank())

p

dev.off()

data$group = factor(data$group, levels = c("NRMSE","\u03c1", "Sign Acc"))
m1 = split(data, data$group)
m1 = sapply(1:length(m1), function(x) mean(m1[[x]]$value,na.rm=T))
data2$group = factor(data2$group, levels = c("NRMSE","\u03c1", "Sign Acc"))
m2 = split(data2, data2$group)
m2 = sapply(1:length(m2), function(x) mean(m2[[x]]$value,na.rm=T))

m1
m2
