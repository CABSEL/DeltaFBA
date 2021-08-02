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

setwd(paste0(here(),"/Simulations/Myocyte/"))

df2 = read.xlsx("Results_5_2_GSE19420.xlsx", 1, header = F)
df3 = read.xlsx("Results_5_2_GSE25462.xlsx", 1, header = F)
delta = data.frame(GSE19420 = df2$X1, GSE25462 = df3$X1)


df2 = read.xlsx("Results_5_2_GSE19420.xlsx", 2, header = F)
df3 = read.xlsx("Results_5_2_GSE25462.xlsx", 2, header = F)
rxnexpr = data.frame(GSE19420 = df2$X1, GSE25462 = df3$X1)
Cond = c("van Tienen et al.", "Jin et al.")
setwd(paste0(here(),"/Simulations/Myocyte/"))
ReactionData = readMat("model.mat")
ReactionData = lapply(ReactionData$model[,,1], unlist)
ReactionData = data.frame(ID = ReactionData$rxns, Subsystems = ReactionData$subSystems, stringsAsFactors = F)

Tres = c()
eps = 10

Data = as.data.frame(table(ReactionData$Subsystems), stringsAsFactors = F)
colnames(Data) = c("Subsystems", "Total_Reactions")

for (j in 1:length(Cond)) {
  flux = data.frame(ID = ReactionData$ID, Flux = delta[,j], stringsAsFactors = F)
  
  up_flux = as.character(flux$ID[flux$Flux>=eps])
  select = ReactionData[ReactionData$ID %in% up_flux,]
  selectData = as.data.frame(table(select$Subsystems), stringsAsFactors = F)
  Data$Up = 0
  for (i in 1:length(Data$Subsystems)) {
    ind = which(selectData$Var1==Data$Subsystems[i])
    if (length(ind)>0) {
      Data$Up[i] = selectData$Freq[ind]
    }
  }
  Res = data.frame(Subsystems=Data$Subsystems, stringsAsFactors = F)
  Res$pval = 0
  Res$FDR = 0
  Res$Odds = 0
  for (i in 1:length(Res$Subsystems)) {
    f1 = fisher.test(matrix(c(Data$Up[i], sum(Data$Up)-Data$Up[i], Data$Total_Reactions[i] - Data$Up[i], sum(Data$Total_Reactions)-Data$Total_Reactions[i]-sum(Data$Up)+Data$Up[i]),2,2), alternative ="greater")
    Res$pval[i] = f1$p.value
    Res$Odds[i] = f1$estimate
  }
  Res$FDR = p.adjust(Res$pval, method = "BH")
  Res = Res[Res$FDR<0.05,]
  Res = Res[Res$Odds>1,]
  if (length(Res$Subsystems)>0) {
    Res$Cluster = Cond[j]
    Res$Reg = "Up-reactions"
    Tres = rbind(Tres,Res)
  }
  
  down_flux = as.character(flux$ID[flux$Flux<=(-eps)])
  select = ReactionData[ReactionData$ID %in% down_flux,]
  selectData = as.data.frame(table(select$Subsystems), stringsAsFactors = F)
  Data$Down = 0
  for (i in 1:length(Data$Subsystems)) {
    ind = which(selectData$Var1==Data$Subsystems[i])
    if (length(ind)>0) {
      Data$Down[i] = selectData$Freq[ind]
    }
  }
  Res = data.frame(Subsystems=Data$Subsystems, stringsAsFactors = F)
  Res$pval = 0
  Res$FDR = 0
  Res$Odds = 0
  
  for (i in 1:length(Res$Subsystems)) {
    f2 = fisher.test(matrix(c(Data$Down[i], sum(Data$Down)-Data$Down[i], Data$Total_Reactions[i] - Data$Down[i], sum(Data$Total_Reactions)-Data$Total_Reactions[i]-sum(Data$Down)+Data$Down[i]),2,2), alternative ="two.sided")
    Res$pval[i] = f2$p.value
    Res$Odds[i] = f2$estimate
  }
  Res$FDR = p.adjust(Res$pval, method = "BH")
  Res = Res[Res$FDR<0.05,]
  Res = Res[Res$Odds>1,]
  if (length(Res$Subsystems)>0) {
    Res$Cluster = Cond[j]
    Res$Reg = "Down-reactions"
    Tres = rbind(Tres,Res)
  }
}

Tres$Neg = -log(Tres$FDR)
Tres$Odds[which(Tres$Odds>100)] = 100
Tres = Tres[order(Tres$Neg),]
Tres$Cluster = factor(Tres$Cluster, levels = Cond)
Tres$Subsystems = str_wrap(Tres$Subsystems, width = 80)
Tres$Subsystems = factor(Tres$Subsystems, levels = unique(Tres$Subsystems))

library(extrafont)
library(viridis)
#font_import()
loadfonts(device = "win")
setwd(paste0(here(),"/Plotting/Myocyte/"))
tiff('T2D_5_2.tiff', units="px", width=(1000*3), height=(625*3), res=300)
p <- ggplot(Tres,
            aes_(x = ~Cluster, y = ~Subsystems, size = ~Neg))

p <- p +
  geom_point(shape = 16, alpha = 0.9) +
  aes_string(color="Odds") +
  
  theme_bw(base_size = 12) +
  scale_x_discrete(drop=FALSE)+
  # scale_colour_gradient(limits=c(0, 0.05), low="red", na.value = "black")
  scale_color_viridis(alpha = 0.7, begin = 0, end = 0.9,breaks = c(min(Tres$Odds), 50, max(Tres$Odds)), labels = c("10", "50", "100"), name = "Odds Ratio\n")+
  #scale_fill_gradientn(colours = topo.colors(7), name = "Fraction (Down/Total)")+
  #scale_color_distiller(palette = "RdYlGn", limits = c(0,90), direction = 1, name = "Fraction of Downregulated genes\n", guide = F)+
  
  scale_size(range = c(4,7),breaks = c(min(Tres$Neg), 40
                                       , max(Tres$Neg)), labels  = c("4", "40", "80"),name = "- Log (Adj.P)") + 
  
  facet_grid(~Reg,scales="free",space="free")+
  theme(strip.text.x = element_text(family="Calibri", face = "bold", size = 18)) +
  theme(panel.spacing =unit(0.5, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, size = 1.5), 
        strip.background = element_rect(color = "#476b6b", size = 1.5, fill = "#d6dce4"))+
  theme(axis.ticks.y = element_line(size=1,color='#476b6b'))+
  theme(axis.ticks.x = element_blank())

p <- p + xlab("") + ylab("") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70"))
p<- p+theme(axis.title = element_text(size = 14, face = "bold", family = "Calibri"))
p<-p+theme(axis.text.y = element_text(size=16,face = "bold", family = "Calibri"), axis.text.x = element_text(size=16,face = "bold", family = "Calibri", angle = 90, hjust = 1, vjust = 0.5))
p<- p+ theme(legend.text=element_text(size=15, family = "Calibri"))
p<- p + theme(legend.title = element_text(size = 16, family = "Calibri", face = "bold"))
p = p+theme(legend.position="bottom")
#axis.text.x=element_blank(),
#axis.ticks.x=element_blank()

p
dev.off()






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

setwd(paste0(here(),"/Simulations/Myocyte/"))

df2 = read.xlsx("Results_25_1_GSE19420.xlsx", 1, header = F)
df3 = read.xlsx("Results_25_1_GSE25462.xlsx", 1, header = F)
delta = data.frame(GSE19420 = df2$X1, GSE25462 = df3$X1)


df2 = read.xlsx("Results_25_1_GSE19420.xlsx", 2, header = F)
df3 = read.xlsx("Results_25_1_GSE25462.xlsx", 2, header = F)
rxnexpr = data.frame(GSE19420 = df2$X1, GSE25462 = df3$X1)
Cond = c("van Tienen et al.", "Jin et al.")
setwd(paste0(here(),"/Simulations/Myocyte/"))
ReactionData = readMat("model.mat")
ReactionData = lapply(ReactionData$model[,,1], unlist)
ReactionData = data.frame(ID = ReactionData$rxns, Subsystems = ReactionData$subSystems, stringsAsFactors = F)

Tres = c()
eps = 10

Data = as.data.frame(table(ReactionData$Subsystems), stringsAsFactors = F)
colnames(Data) = c("Subsystems", "Total_Reactions")

for (j in 1:length(Cond)) {
  flux = data.frame(ID = ReactionData$ID, Flux = delta[,j], stringsAsFactors = F)
  
  up_flux = as.character(flux$ID[flux$Flux>=eps])
  select = ReactionData[ReactionData$ID %in% up_flux,]
  selectData = as.data.frame(table(select$Subsystems), stringsAsFactors = F)
  Data$Up = 0
  for (i in 1:length(Data$Subsystems)) {
    ind = which(selectData$Var1==Data$Subsystems[i])
    if (length(ind)>0) {
      Data$Up[i] = selectData$Freq[ind]
    }
  }
  Res = data.frame(Subsystems=Data$Subsystems, stringsAsFactors = F)
  Res$pval = 0
  Res$FDR = 0
  Res$Odds = 0
  for (i in 1:length(Res$Subsystems)) {
    f1 = fisher.test(matrix(c(Data$Up[i], sum(Data$Up)-Data$Up[i], Data$Total_Reactions[i] - Data$Up[i], sum(Data$Total_Reactions)-Data$Total_Reactions[i]-sum(Data$Up)+Data$Up[i]),2,2), alternative ="greater")
    Res$pval[i] = f1$p.value
    Res$Odds[i] = f1$estimate
  }
  Res$FDR = p.adjust(Res$pval, method = "BH")
  Res = Res[Res$FDR<0.05,]
  Res = Res[Res$Odds>1,]
  if (length(Res$Subsystems)>0) {
    Res$Cluster = Cond[j]
    Res$Reg = "Up-reactions"
    Tres = rbind(Tres,Res)
  }
  
  down_flux = as.character(flux$ID[flux$Flux<=(-eps)])
  select = ReactionData[ReactionData$ID %in% down_flux,]
  selectData = as.data.frame(table(select$Subsystems), stringsAsFactors = F)
  Data$Down = 0
  for (i in 1:length(Data$Subsystems)) {
    ind = which(selectData$Var1==Data$Subsystems[i])
    if (length(ind)>0) {
      Data$Down[i] = selectData$Freq[ind]
    }
  }
  Res = data.frame(Subsystems=Data$Subsystems, stringsAsFactors = F)
  Res$pval = 0
  Res$FDR = 0
  Res$Odds = 0
  
  for (i in 1:length(Res$Subsystems)) {
    f2 = fisher.test(matrix(c(Data$Down[i], sum(Data$Down)-Data$Down[i], Data$Total_Reactions[i] - Data$Down[i], sum(Data$Total_Reactions)-Data$Total_Reactions[i]-sum(Data$Down)+Data$Down[i]),2,2), alternative ="two.sided")
    Res$pval[i] = f2$p.value
    Res$Odds[i] = f2$estimate
  }
  Res$FDR = p.adjust(Res$pval, method = "BH")
  Res = Res[Res$FDR<0.05,]
  Res = Res[Res$Odds>1,]
  if (length(Res$Subsystems)>0) {
    Res$Cluster = Cond[j]
    Res$Reg = "Down-reactions"
    Tres = rbind(Tres,Res)
  }
}

Tres$Neg = -log(Tres$FDR)
Tres$Odds[which(Tres$Odds>100)] = 100
Tres = Tres[order(Tres$Neg),]
Tres$Cluster = factor(Tres$Cluster, levels = Cond)
Tres$Subsystems = str_wrap(Tres$Subsystems, width = 80)
Tres$Subsystems = factor(Tres$Subsystems, levels = unique(Tres$Subsystems))

library(extrafont)
library(viridis)
#font_import()
loadfonts(device = "win")
setwd(paste0(here(),"/Plotting/Myocyte/"))
tiff('T2D_25_1.tiff', units="px", width=(1300*3), height=(820*3), res=300)
p <- ggplot(Tres,
            aes_(x = ~Cluster, y = ~Subsystems, size = ~Neg))

p <- p +
  geom_point(shape = 16, alpha = 0.9) +
  aes_string(color="Odds") +
  
  theme_bw(base_size = 12) +
  scale_x_discrete(drop=FALSE)+
  # scale_colour_gradient(limits=c(0, 0.05), low="red", na.value = "black")
  scale_color_viridis(alpha = 0.7, begin = 0, end = 0.9,breaks = c(min(Tres$Odds), 50, max(Tres$Odds)), labels = c("4", "50", "100"), name = "Odds Ratio\n")+
  #scale_fill_gradientn(colours = topo.colors(7), name = "Fraction (Down/Total)")+
  #scale_color_distiller(palette = "RdYlGn", limits = c(0,90), direction = 1, name = "Fraction of Downregulated genes\n", guide = F)+
  
  scale_size(range = c(4,7),breaks = c(min(Tres$Neg), 40, max(Tres$Neg)), labels  = c("4", "50", ">100"),name = "- Log (Adj.P)") + 
  
  facet_grid(~Reg,scales="free",space="free")+
  theme(strip.text.x = element_text(family="Calibri", face = "bold", size = 18)) +
  theme(panel.spacing =unit(0.5, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, size = 1.5), 
        strip.background = element_rect(color = "#476b6b", size = 1.5, fill = "#d6dce4"))+
  theme(axis.ticks.y = element_line(size=1,color='#476b6b'))+
  theme(axis.ticks.x = element_blank())

p <- p + xlab("") + ylab("") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70"))
p<- p+theme(axis.title = element_text(size = 14, face = "bold", family = "Calibri"))
p<-p+theme(axis.text.y = element_text(size=16,face = "bold", family = "Calibri"), axis.text.x = element_text(size=16,face = "bold", family = "Calibri", angle = 90, hjust = 1, vjust = 0.5))
p<- p+ theme(legend.text=element_text(size=15, family = "Calibri"))
p<- p + theme(legend.title = element_text(size = 16, family = "Calibri", face = "bold"))
p = p+theme(legend.position="bottom")
#axis.text.x=element_blank(),
#axis.ticks.x=element_blank()

p
dev.off()
