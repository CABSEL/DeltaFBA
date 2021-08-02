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
library(reshape2)

setwd(paste0(here(),"/Data/Gerosa/"))
wd = read.xlsx("Fluxomics.xlsx", 1)
conditions = colnames(wd)[-c(1)] # Carbon sources
fluxes = as.character(wd[,1])

setwd(paste0(here(),"/Simulations/Gerosa/"))

Cor = matrix(0, nrow = length(conditions), ncol = length(conditions))
rownames(Cor) = conditions
colnames(Cor) = conditions
d_cor = Cor
d_err = Cor
accuracy = Cor
r.sq = Cor
r.sqratio = Cor


for (j in 1:length(conditions)) {
  
  wtflux = wd[,j+1]
  M = wd[,-c(1)]
  
  Md = read.xlsx(paste("Results",j,".xlsx", sep = "", collapse = ""),1, header = F)
  Pd = read.xlsx(paste("Results",j,".xlsx", sep = "", collapse = ""),2, header = F)
  Rd = read.xlsx(paste("Results",j,".xlsx", sep = "", collapse = ""),3, header = F)
  
  MRatio = M/wtflux
  MRatio = as.matrix(MRatio)
  MRatio[is.nan(MRatio)] = 0
  MRatio[which(MRatio==Inf)] = 100
  MRatio[which(MRatio==-Inf)] = 0.001
  
  for (i in 1:length(conditions)) {
    
    P = Pd[,i] + wtflux
    sel = intersect(which(P!=0), which(M[,i]!=0))
    Cor[j,i] =  P[sel] %*%  M[sel,i] / (norm(as.matrix(P[sel]), type = "2") * norm(as.matrix(M[sel,i]), type = "2"))
    sel = intersect(which(Pd[,i]!=0), which(Md[,i]!=0))
    
    if (length(sel)>0)  {
      d_cor[j,i] =  Pd[sel,i] %*%  Md[sel,i] / (norm(as.matrix(Pd[sel,i]), type = "2") * norm(as.matrix(Md[sel,i]), type = "2"))
      d_err[j,i] = tdStats(as.numeric(Pd[,i]), as.numeric(Md[,i]), functions= c("nrmse"))
    }
    
    s_Md = as.numeric(Md[,i])
    s_Md[s_Md>0] = 1;
    s_Md[s_Md<0] = -1;
    
    s_Pd = as.numeric(Pd[,i])
    s_Pd[s_Pd>0] = 1;
    s_Pd[s_Pd<0] = -1;
    accuracy[j,i] = accuracy(s_Pd, s_Md)
    
    r.sq[j,i] = summary(lm(Rd[,i]~Md[,i]))$r.squared
    r.sqratio[j,i] = summary(lm(Rd[,i]~MRatio[,i]))$r.squared
  }
}


df = melt(d_cor)
df2 = melt(d_err)
df3 = melt(accuracy)
v = c()
t = df
t$value[t$value==0] = NA
v[1] = mean(t$value, na.rm = T)
t = df2
t$value[t$value==0] = NA
v[2] = mean(t$value, na.rm = T)
t = df3
t$value[t$value==1] = NA
v[3] =mean(t$value, na.rm = T)
v = as.matrix(v)
rownames(v) = c("Rho", "NRMSE","Accuracy")
setwd(paste0(here(),"/Plotting/Gerosa/"))
write.table(v, "Gerosa_Avg.txt", row.names = T, col.names = F, sep = "\t")

df$value[df$value==0] = NA
df$value2 = -log(df2$value)
df$value2[is.infinite(df$value2)] = NA
df$value2[is.na(df$value)] = NA
df$value3 = df3$value
df$text = round(df$value3,2)
df$hjust = 0.5
df$vjust = 0.5
df2 = df
df3 = df
df2 = df2[-which(df2$text==1),]
df3 = df3[which(df3$text==1),]


tiff('Gerosa.tiff', units="px", width=(925*3), height=(650*3), res=300)
p <- ggplot(df,
            aes_(x = ~Var1, y = ~Var2, size = ~value2)) +
  
  geom_point(shape = 16, alpha = 0.9) +
  aes_string(color="value") +
  
  theme_bw(base_size = 12) +
  scale_color_gradient2(midpoint = 0,low = "#784809", mid = "#e5e5e5", high = "#004742", name = expression(bold(rho)), breaks = c(-1,0,1),limits=c(-1,1), labels = c("-1", "0", "1")) +
  scale_size(range = c(5,10),name = "- Log (NRMSE)", breaks = c(min(df$value2,na.rm = T),mean(df$value2, na.rm= T), max(df$value2, na.rm = T)), labels = c("1.5", "2.0", "2.75")) + 
  theme(axis.ticks.y = element_line(size=1,color='#476b6b')) +
  theme(axis.ticks.x = element_line(size=1,color='#476b6b')) +
  xlab("Reference Condition") + ylab("Experimental Condition") + ggtitle("") + 
  theme(axis.title = element_text(size = 16, face = "bold", family = "Calibri")) + 
  theme(axis.text.y = element_text(size=12, family = "Calibri"), axis.text.x = element_text(size=12, family = "Calibri", angle = 90, vjust = 0.2)) + 
  theme(legend.text=element_text(size=12, family = "Calibri")) + 
  theme(legend.title = element_text(size = 14, family = "Calibri", face = "bold")) +
  theme(panel.border = element_blank()) +
  geom_text(data=df2, aes(x=Var1, y=Var2, label=as.character(text), hjust=hjust, vjust = vjust), color="peachpuff",family = "Calibri", fontface="bold",alpha=1, size=3, inherit.aes = FALSE ) +
  geom_text(data=df3, aes(x=Var1, y=Var2, label=as.character(text), hjust=hjust, vjust = vjust), color="black",family = "Calibri", fontface="bold",alpha=1, size=3, inherit.aes = FALSE )

p    
dev.off()