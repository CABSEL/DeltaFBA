rm(list = ls(all=T))
graphics.off()
library(here)
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
conditions = colnames(wd)[-c(1,2)]
conditions[1] = "D = 0.1/h"
conditions[2] = "D = 0.4/h"
conditions[3] = "D = 0.5/h"
conditions[4] = "D = 0.7/h"
fluxes = as.character(wd[,1])
wtflux = wd[,2]
M = wd[,-c(1,2)]

setwd(paste0(here(),"/Simulations/Ishii/Sensitivity/Ishii_RT_PCR_Relaxed/"))
data.files = list.files(pattern = "*.xlsx")
eps = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2.5,5,10)
res = c()
for (p in 1:length(eps)) {
  
  Md = read.xlsx(paste("Results_", eps[[p]],".xlsx", sep="", collapse=""),1, header = F)
  Pd = read.xlsx(paste("Results_", eps[[p]],".xlsx", sep="", collapse=""),2, header = F)
  Rd = read.xlsx(paste("Results_", eps[[p]],".xlsx", sep="", collapse=""),3, header = F)
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
    condition=rep(conditions,4),
    group=c(rep('NRMSE', length(conditions)), rep('Rho', length(conditions)),rep('Acc', length(conditions)), rep('R^2', length(conditions))) ,
    value= c(m1, m2,m3, m6))
  data$Eps = eps[[p]]
  
  res = rbind(res,data)
  
}

data = res
data$group = factor(data$group, levels = c("NRMSE", "Rho","Acc",'R^2'))
data = data[order(data$group),]
data$Eps = factor(data$Eps, levels = eps)
data$Method = "dFBA"


rm(list= ls()[!(ls() %in% c('data', 'p_cond','eps'))])

p1 = data[data$group == "NRMSE",]
p1$value = log10(p1$value)
p2 = data[data$group == "Rho",]
p3 = data[data$group == "Acc",]
dummyp1 = p1[c(1,2),]
dummyp1$value = c(-2,1.2)
dummyp1 = droplevels(dummyp1)
dummyp1$group = factor(dummyp1$group, levels = c("NRMSE","Rho", "Acc"))

dummyp2 = p2[c(1,2),]
dummyp2$value = c(-0.5,1)
dummyp2 = droplevels(dummyp2)
dummyp2$group = factor(dummyp2$group, levels = c("NRMSE","Rho", "Acc"))

dummyp3 = p3[c(1,2),]
dummyp3$value = c(0.2,0.9)
dummyp3 = droplevels(dummyp3)
dummyp3$group = factor(dummyp3$group, levels = c("NRMSE","Rho", "Acc"))


dat = rbind(p1,p2,p3)
dat$group = factor(dat$group, levels = c("NRMSE","Rho", "Acc"))
cols = rep("#0b84a5", times  = length(eps))

setwd(paste0(here(),"/Plotting/Ishii/"))

tiff('Sensitivity_Panel.tiff', units="px", width=(1200*3), height=(850*3), res=300)

p <- ggplot(dat, aes(x = Eps, y = value, fill= Eps)) + 
  geom_boxplot(outlier.size = 2) + 
  stat_boxplot(geom = "errorbar",width = 0.4,size = 0.75) + 
  geom_boxplot(lwd=0.55) + 
  theme_bw(base_size = 12) +
  geom_blank(data = dummyp1) +
  geom_blank(data = dummyp2) +
  geom_blank(data = dummyp3) +
  facet_grid(group~.,scales="free")+
  theme(strip.text.y = element_blank())+
  theme(panel.spacing =unit(1, "lines"),
        panel.border = element_rect( fill = NA, size = 1)) +
  scale_fill_manual(values = cols) +
  theme(legend.position = "none") +
  #ylim(-1,1) +
  xlab("Epsilon") + ylab("") +
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold",angle = 90, hjust = 1, vjust = 0.4))

p

dev.off()
