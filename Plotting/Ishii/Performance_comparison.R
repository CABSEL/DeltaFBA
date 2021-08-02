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

setwd(paste0(here(),"/OtherTools/transcript2flux-master/results/"))

fil = list.files()
fil = fil[grep("ishii_sim_all", fil)]

ncon = 28

nrmse = matrix(0, nrow = ncon, ncol = length(fil))
colnames(nrmse) = str_extract(fil, "[^_]+")

rho = nrmse
acc = nrmse
acc_med = nrmse
acc_75 = nrmse
zero_l = list()

for (i in 1:length(fil)) {
  df = readMat(fil[i])
  df = df[[1]]
  f_exp = df[[7]]
  f_p = df[[8]]
  
  m1 = c()
  m2 = c()
  m3 = c()
  m4 = c()
  m5 = c()
  
  zero = matrix(0, nrow = 3, ncol = (length(f_exp)-1))
  
  for (j in 2:length(f_exp)) {
    sim = tryCatch(f_p[[j]][[1]] - f_p[[1]][[1]], error = function(e) NA)
    exp = f_exp[[j]][[1]] - f_exp[[1]][[1]]
    
    m1[j-1] = tryCatch(tdStats(sim, exp, functions= c("nrmse")) , error  = function(e) NA)
    
    sim = tryCatch(sim[,1], error = function(e) NA)
    exp = tryCatch(exp[,1], error = function(e) NA)
    m2[j-1] =  tryCatch(sim %*%  exp / (norm(as.matrix(sim), type = "2") * norm(as.matrix(exp), type = "2")), error = function(e) NA)
    
    
    s_Md = as.numeric(exp)
    sel = which(abs(s_Md)>median(abs(s_Md)))
    sel2 = which(abs(s_Md)>quantile(abs(s_Md),0.75))
    s_Md[s_Md>0] = 1;
    s_Md[s_Md<0] = -1;
    
    s_Pd = as.numeric(sim)
    s_Pd[s_Pd>0] = 1;
    s_Pd[s_Pd<0] = -1;
    
    m3[j-1] = accuracy(s_Md,s_Pd)
    m4[j-1] = accuracy(s_Md[sel], s_Pd[sel])
    m5[j-1] = accuracy(s_Md[sel2], s_Pd[sel2])
    
    zero[1,j-1] = length(which(exp==0))
    zero[2,j-1] = length(which(sim==0))
    zero[3,j-1] = length(intersect(which(exp==0), which(sim==0)))
    
  }
  
  nrmse[,i] = m1
  rho[,i] = m2
  acc[,i] = m3
  acc_med[,i] = m4
  acc_75[,i] = m5
  zero_l[[i]] = zero
}
zero_dat = data.frame(Method = str_extract(fil, "[^_]+"),
                      Exp = sapply(1:length(zero_l), function(x) mean(zero_l[[x]][1,])),
                      Sim = sapply(1:length(zero_l), function(x) mean(zero_l[[x]][2,])),
                      Int = sapply(1:length(zero_l), function(x) mean(zero_l[[x]][3,])))


rm(list= ls()[!(ls() %in% c('nrmse','rho','acc','acc_med','acc_75', 'zero_dat'))])

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
m7 = c()
zero = matrix(0, nrow = 3, ncol = length(conditions))

for (i in 1:length(conditions)) {
  sim = Pd[,i]
  exp = Md[,i]
  
  m1[i] = tryCatch(tdStats(sim, exp, functions= c("nrmse")) , error  = function(e) NA)
  
  m2[i] =  tryCatch(sim %*%  exp / (norm(as.matrix(sim), type = "2") * norm(as.matrix(exp), type = "2")), error = function(e) NA)
  
  s_Md = as.numeric(exp)
  sel = which(abs(s_Md)>median(abs(s_Md)))
  sel2 = which(abs(s_Md)>quantile(abs(s_Md),0.75))
  sel3 = which(Rd[,i]!=0)
  s_Md[s_Md>0] = 1;
  s_Md[s_Md<0] = -1;
  
  s_Pd = as.numeric(sim)
  s_Pd[s_Pd>0] = 1;
  s_Pd[s_Pd<0] = -1;
  
  m3[i] = accuracy(s_Md,s_Pd)
  m4[i] = accuracy(s_Md[sel], s_Pd[sel])
  m5[i] = accuracy(s_Md[sel2], s_Pd[sel2])
  
  m6[i] = summary(lm(MRatio[,i]~Rd[,i]))$r.squared
  m7[i] = accuracy(s_Md[sel3], s_Pd[sel3])
  
  zero[1,i] = length(which(exp==0))
  zero[2,i] = length(which(sim==0))
  zero[3,i] = length(intersect(which(exp==0), which(sim==0)))
  
}

nrmse = cbind(nrmse,m1)
rho = cbind(rho,m2)
acc = cbind(acc,m3)
acc_med = cbind(acc_med,m4)
acc_75 = cbind(acc_75,m5)
dfbam6 = m6
adddat = data.frame(Method = "dFBA", Exp = mean(zero[1,]), Sim = mean(zero[2,]), Int = mean(zero[3,]))
zero_dat = rbind(zero_dat, adddat)

rm(list= ls()[!(ls() %in% c('nrmse','rho','acc','acc_med','acc_75','dfbam6','conditions','zero_dat','m7'))])

setwd(paste0(here(),"/Data/Ishii/"))
wd = read.xlsx("Fluxomics.xlsx", 1)
conditions = colnames(wd)[-c(1,2)]
conditions[1] = "WT_0.1h-1"
conditions[2] = "WT_0.4h-1"
conditions[3] = "WT_0.5h-1"
conditions[4] = "WT_0.7h-1"
fluxes = as.character(wd[,1])
wtflux = wd[,2]
M = wd[,-c(1,2)]

setwd(paste0(here(),"/Simulations/Ishii/RT_PCR_REMI/"))
Md = matrix(0, nrow= length(fluxes), ncol = length(conditions))
Pd = Md
Pmut = Pd
Rd = Pd
MRatio = Md
for (j in 1:length(conditions)) {
  M_wt = read.xlsx(paste(conditions[j], ".xlsx", sep = ""), 1, header = F)
  M_exp = read.xlsx(paste(conditions[j], ".xlsx", sep = ""), 2, header = F)
  P_wt = read.xlsx(paste(conditions[j], ".xlsx", sep = ""), 3, header = F)
  P_exp = read.xlsx(paste(conditions[j], ".xlsx", sep = ""), 4, header = F)
  R_d = read.xlsx(paste(conditions[j], ".xlsx", sep = ""), 5, header = F)
  sol = read.xlsx(paste(conditions[j], ".xlsx", sep = ""), 9, header = F)
  sol = as.numeric(sol)
  ind = which(sol == min(sol))
  Md[,j] = M_exp[,ind] - M_wt[,ind]
  MRatio[,j] = M_exp[,ind] / M_wt[,ind]
  Pd[,j] = P_exp[,ind] - P_wt[,ind]
  Pmut[,j] = P_exp[,ind]
  Rd[,j] = as.numeric(R_d$X1)
}
MRatio[is.nan(MRatio)] = 0
MRatio[which(MRatio==Inf)] = 100
MRatio[which(MRatio==-Inf)] = 0.001

m1 = c()
m2 = c()
m3 = c()
m4 = c()
m5 = c()
m6 = c()
zero = matrix(0, nrow = 3, ncol = length(conditions))


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
  zero[1,i] = length(which(exp==0))
  zero[2,i] = length(which(sim==0))
  zero[3,i] = length(intersect(which(exp==0), which(sim==0)))
  
}

nrmse = cbind(nrmse,m1)
rho = cbind(rho,m2)
acc = cbind(acc,m3)
acc_med = cbind(acc_med,m4)
acc_75 = cbind(acc_75,m5)
remim6 = m6
adddat = data.frame(Method = "REMI", Exp = mean(zero[1,]), Sim = mean(zero[2,]), Int = mean(zero[3,]))
zero_dat = rbind(zero_dat, adddat)

rm(list= ls()[!(ls() %in% c('nrmse','rho','acc','acc_med','acc_75','dfbam6','remim6','conditions','zero_dat','m7'))])

cols =  c("#4379a7", rep("#d3d3d3",9))
setwd(paste0(here(),"/Plotting/Ishii/"))


colnames(nrmse) = c(colnames(nrmse)[1:8],"dFBA","REMI")
dat = melt(nrmse)
dat$Var2 = factor(dat$Var2, levels = c("dFBA","REMI","pFBA","GIMME","iMAT","MADE", "Lee-12","RELATCH", "GX-FBA","E-Flux"))
p1 = dat
dat = dat[!is.nan(dat$value),]

dat$Var2 = as.character(dat$Var2)
dat$Var2[which(dat$Var2=="dFBA")] = paste0("\u0394", "FBA")
dat$Var2 = factor(dat$Var2, levels = c(paste0("\u0394", "FBA"),"REMI","pFBA","GIMME","iMAT","MADE", "Lee-12","RELATCH", "GX-FBA","E-Flux"))

tiff('Log_NRMSE.tiff', units="px", width=(800*3), height=(575*3), res=300)
p <- ggplot(dat, aes(x = Var2, y = log10(value), fill= Var2)) + 
  geom_boxplot(outlier.size = 2) + 
  stat_boxplot(geom = "errorbar",width = 0.4,size = 0.75) + 
  geom_boxplot(lwd=0.55) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = cols) +
  theme(legend.position = "none") +
  ylim(-2,2) +
  xlab("Method") + ylab(expression(bold(Log[10]~(NRMSE)))) +
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold",angle = 90, hjust = 1, vjust = 0.4))

p
dev.off()

p1$value = log10(p1$value)
t1 = split(p1,p1$Var2)
nrmse_mean = sapply(1:length(t1), function(x) mean(10^(t1[[x]]$value),na.rm=T))
nrmse_p = c()
for (i in 1:(length(t1)-1)) {
  #p1p[i] = wilcox.test(t1[[1]]$value,t1[[i+1]]$value, alternative = "less")$p.value
  nrmse_p[i] = t.test(t1[[1]]$value,t1[[i+1]]$value, paired = TRUE, alternative = "two.sided")$p.value
}



colnames(rho) = c(colnames(rho)[1:8],"dFBA","REMI")
dat = melt(rho)
dat$Var2 = factor(dat$Var2, levels = c("dFBA","REMI","pFBA","GIMME","iMAT","MADE", "Lee-12","RELATCH", "GX-FBA","E-Flux"))
p2 = dat
dat = dat[!is.na(dat$value),]

dat$Var2 = as.character(dat$Var2)
dat$Var2[which(dat$Var2=="dFBA")] = paste0("\u0394", "FBA")
dat$Var2 = factor(dat$Var2, levels = c(paste0("\u0394", "FBA"),"REMI","pFBA","GIMME","iMAT","MADE", "Lee-12","RELATCH", "GX-FBA","E-Flux"))

tiff('Rho.tiff', units="px", width=(800*3), height=(575*3), res=300)
p <- ggplot(dat, aes(x = Var2, y = value, fill= Var2)) + 
  geom_boxplot(outlier.size = 2) + 
  stat_boxplot(geom = "errorbar",width = 0.4,size = 0.75) + 
  geom_boxplot(lwd=0.55) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = cols) +
  theme(legend.position = "none") +
  ylim(-1,1) +
  xlab("Method") + ylab(expression(bold(rho))) +
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold",angle = 90, hjust = 1, vjust = 0.4))

p
dev.off()

p2$value = p2$value
t1 = split(p2,p2$Var2)
rho_mean = sapply(1:length(t1), function(x) mean(t1[[x]]$value,na.rm=T))
rho_p = c()
for (i in 1:(length(t1)-1)) {
  #p1p[i] = wilcox.test(t1[[1]]$value,t1[[i+1]]$value, alternative = "less")$p.value
  rho_p[i] = t.test(t1[[1]]$value,t1[[i+1]]$value, paired = TRUE, alternative = "two.sided")$p.value
}


colnames(acc) = c(colnames(acc)[1:8],"dFBA","REMI")
dat = melt(acc)
dat$Var2 = factor(dat$Var2, levels = c("dFBA","REMI","pFBA","GIMME","iMAT","MADE", "Lee-12","RELATCH", "GX-FBA","E-Flux"))
p3 = dat
dat = dat[!is.na(dat$value),]

dat$Var2 = as.character(dat$Var2)
dat$Var2[which(dat$Var2=="dFBA")] = paste0("\u0394", "FBA")
dat$Var2 = factor(dat$Var2, levels = c(paste0("\u0394", "FBA"),"REMI","pFBA","GIMME","iMAT","MADE", "Lee-12","RELATCH", "GX-FBA","E-Flux"))

tiff('Acc.tiff', units="px", width=(800*3), height=(575*3), res=300)
p <- ggplot(dat, aes(x = Var2, y = value, fill= Var2)) + 
  geom_boxplot(outlier.size = 2) + 
  stat_boxplot(geom = "errorbar",width = 0.4,size = 0.75) + 
  geom_boxplot(lwd=0.55) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = cols) +
  theme(legend.position = "none") +
  ylim(0,1) +
  xlab("Method") + ylab("Sign Acc") +
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold",angle = 90, hjust = 1, vjust = 0.4))

p
dev.off()

p3$value = p3$value
t1 = split(p3,p3$Var2)
acc_mean = sapply(1:length(t1), function(x) mean(t1[[x]]$value,na.rm=T))

acc_p = c()
for (i in 1:(length(t1)-1)) {
  #p1p[i] = wilcox.test(t1[[1]]$value,t1[[i+1]]$value, alternative = "less")$p.value
  acc_p[i] = t.test(t1[[1]]$value,t1[[i+1]]$value, paired = TRUE, alternative = "two.sided")$p.value
}


dat$value[which(dat$Var2==paste0("\u0394", "FBA"))] = m7
tiff('Acc_Mapped.tiff', units="px", width=(800*3), height=(575*3), res=300)
p <- ggplot(dat, aes(x = Var2, y = value, fill= Var2)) + 
  geom_boxplot(outlier.size = 2) + 
  stat_boxplot(geom = "errorbar",width = 0.4,size = 0.75) + 
  geom_boxplot(lwd=0.55) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = cols) +
  theme(legend.position = "none") +
  ylim(0,1) +
  xlab("Method") + ylab("Sign Acc") +
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold",angle = 90, hjust = 1, vjust = 0.4))

p
dev.off()



colnames(acc_med) = c(colnames(acc_med)[1:8],"dFBA","REMI")
dat = melt(acc_med)
dat$Var2 = factor(dat$Var2, levels = c("dFBA","REMI","pFBA","GIMME","iMAT","MADE", "Lee-12","RELATCH", "GX-FBA","E-Flux"))
p4 = dat
dat = dat[!is.na(dat$value),]

dat$Var2 = as.character(dat$Var2)
dat$Var2[which(dat$Var2=="dFBA")] = paste0("\u0394", "FBA")
dat$Var2 = factor(dat$Var2, levels = c(paste0("\u0394", "FBA"),"REMI","pFBA","GIMME","iMAT","MADE", "Lee-12","RELATCH", "GX-FBA","E-Flux"))

tiff('Acc_med.tiff', units="px", width=(800*3), height=(575*3), res=300)
p <- ggplot(dat, aes(x = Var2, y = value, fill= Var2)) + 
  geom_boxplot(outlier.size = 2) + 
  stat_boxplot(geom = "errorbar",width = 0.4,size = 0.75) + 
  geom_boxplot(lwd=0.55) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = cols) +
  theme(legend.position = "none") +
  ylim(0,1.2) +
  xlab("Method") + ylab("Sign Acc") +
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold",angle = 90, hjust = 1, vjust = 0.4))

p
dev.off()

p4$value = p4$value
t1 = split(p4,p4$Var2)
acc_med_mean = sapply(1:length(t1), function(x) mean(t1[[x]]$value,na.rm=T))

acc_med_p = c()
for (i in 1:(length(t1)-1)) {
  #p1p[i] = wilcox.test(t1[[1]]$value,t1[[i+1]]$value, alternative = "less")$p.value
  acc_med_p[i] = t.test(t1[[1]]$value,t1[[i+1]]$value, paired = TRUE, alternative = "two.sided")$p.value
}


colnames(acc_75) = c(colnames(acc_75)[1:8],"dFBA","REMI")
dat = melt(acc_75)
dat$Var2 = factor(dat$Var2, levels = c("dFBA","REMI","pFBA","GIMME","iMAT","MADE", "Lee-12","RELATCH", "GX-FBA","E-Flux"))
p5 = dat
dat = dat[!is.na(dat$value),]

dat$Var2 = as.character(dat$Var2)
dat$Var2[which(dat$Var2=="dFBA")] = paste0("\u0394", "FBA")
dat$Var2 = factor(dat$Var2, levels = c(paste0("\u0394", "FBA"),"REMI","pFBA","GIMME","iMAT","MADE", "Lee-12","RELATCH", "GX-FBA","E-Flux"))

tiff('Acc_75.tiff', units="px", width=(800*3), height=(575*3), res=300)
p <- ggplot(dat, aes(x = Var2, y = value, fill= Var2)) + 
  geom_boxplot(outlier.size = 2) + 
  stat_boxplot(geom = "errorbar",width = 0.4,size = 0.75) + 
  geom_boxplot(lwd=0.55) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = cols) +
  theme(legend.position = "none") +
  ylim(0,1.2) +
  xlab("Method") + ylab("Sign Acc") +
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold",angle = 90, hjust = 1, vjust = 0.4))

p
dev.off()

p5$value = p5$value
t1 = split(p5,p5$Var2)
acc_75_mean = sapply(1:length(t1), function(x) mean(t1[[x]]$value,na.rm=T))

acc_75_p = c()
for (i in 1:(length(t1)-1)) {
  #p1p[i] = wilcox.test(t1[[1]]$value,t1[[i+1]]$value, alternative = "less")$p.value
  acc_75_p[i] = t.test(t1[[1]]$value,t1[[i+1]]$value, paired = TRUE, alternative = "two.sided")$p.value
}


p1$Measure = "NRMSE"
p2$Measure = "Rho"
p3$Measure = "Acc"
dummyp1 = data.frame(Var1 = p1$Var1[c(1,2)], Var2 = p1$Var2[c(1,2)], value = c(-2,2), Measure = rep("NRMSE",2))
dummyp1 = droplevels(dummyp1)
dummyp1$Measure = factor(dummyp1$Measure, levels = c("NRMSE","Rho", "Acc"))

dummyp2 = data.frame(Var1 = p2$Var1[c(1,2)], Var2 = p2$Var2[c(1,2)], value = c(-1,1.3), Measure = rep("Rho",2))
dummyp2 = droplevels(dummyp2)
dummyp2$Measure = factor(dummyp2$Measure, levels = c("NRMSE","Rho", "Acc"))

dummyp3 = data.frame(Var1 = p3$Var1[c(1,2)], Var2 = p3$Var2[c(1,2)], value = c(0,1), Measure = rep("Acc",2))
dummyp3 = droplevels(dummyp3)
dummyp3$Measure = factor(dummyp3$Measure, levels = c("NRMSE","Rho", "Acc"))

dat = rbind(p1,p2,p3)
dat$Measure = factor(dat$Measure, levels = c("NRMSE","Rho", "Acc"))
dat$Var2 = as.character(dat$Var2)
dat = dat[!is.na(dat$value),]
dat$Var2[which(dat$Var2=="dFBA")] = paste0("\u0394", "FBA")
dat$Var2 = factor(dat$Var2, levels = c(paste0("\u0394", "FBA"),"REMI","pFBA","GIMME","iMAT","MADE", "Lee-12","RELATCH", "GX-FBA","E-Flux"))

tiff('Fig1.tiff', units="px", width=(720*3), height=(900*3), res=300)

p <- ggplot(dat, aes(x = Var2, y = value, fill= Var2)) + 
  geom_boxplot(outlier.size = 2) + 
  stat_boxplot(geom = "errorbar",width = 0.4,size = 0.75) + 
  geom_boxplot(lwd=0.55) + 
  theme_bw(base_size = 12) +
  geom_blank(data = dummyp1) +
  geom_blank(data = dummyp2) +
  geom_blank(data = dummyp3) +
  facet_grid(Measure~.,scales="free")+
  theme(strip.text.y = element_blank())+
  theme(panel.spacing =unit(1, "lines")) +
  scale_fill_manual(values = cols) +
  theme(legend.position = "none") +
  #ylim(-1,1) +
  xlab("Method") + ylab("") +
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold",angle = 90, hjust = 1, vjust = 0.4))

p

dev.off()

mean(dfbam6)
sd(dfbam6)
nrmse_p
rho_p
acc_p
acc_med_p
acc_75_p

nrmse_mean
rho_mean
acc_mean
acc_med_mean
acc_75_mean
