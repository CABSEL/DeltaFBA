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
delta = read.xlsx("Results_25_1_GSE19420.xlsx", 1, header = F)
setwd(paste0(here(),"/Simulations/Myocyte//"))
ReactionData = readMat("model.mat")
ReactionData = lapply(ReactionData$model[,,1], unlist)

S = as.matrix(ReactionData$S)
mets =  ReactionData$mets
Names = ReactionData$metNames
Comp = as.numeric(ReactionData$metComps)

Compartments = data.frame(Names = ReactionData$compNames, Abbr = ReactionData$comps, stringsAsFactors = F)
Comp = as.character(Compartments$Abbr[Comp])

MetNames = paste(Names, "[", Comp, "]", sep ="")

ReactionData = data.frame(ID = ReactionData$rxns, Subsystems = ReactionData$subSystems, stringsAsFactors = F)
rownames(S) = MetNames
colnames(S) = ReactionData$ID

Data = as.data.frame(table(ReactionData$Subsystems))
Data$Var1 = as.character(Data$Var1)
colnames(Data) =  c("SubSytems","Number of DE reactions")


testmetNames = unique(Names)
nc = c()
np = c()
MetP = matrix(0, nrow = length(Data$SubSytems), ncol = length(testmetNames))
MetC = matrix(0, nrow = length(Data$SubSytems), ncol = length(testmetNames))
creactions = list()
preactions = list()


for (i in 1:length(testmetNames)) {
  allIdx = which(Names==testmetNames[[i]])
  test = as.matrix(S[allIdx,])
  allIdx1 = as.numeric(unlist(lapply(1:nrow(test), function(x) which(test[x,]>0))))
  allIdx2 = as.numeric(unlist(lapply(1:nrow(test), function(x) which(test[x,]<0))))
  nc[i] = length(allIdx2)
  np[i] = length(allIdx1)
  test1 = setdiff(allIdx1, allIdx2)
  
  if (length(test1)>0) {
    preactions[[i]] = cbind(rep(testmetNames[[i]], times = length(test1)), ReactionData[test1,], delta[test1,])
    colnames(preactions[[i]]) = c("MetNames",colnames(preactions[[i]])[2:3], "Flux")
    testP = aggregate(Flux~Subsystems, data = preactions[[i]], FUN = sum)
    MetP[match(testP$Subsystems, Data$SubSytems),i] = testP$Flux
  }
  
  
  test2 = setdiff(allIdx2, allIdx1)
  
  if (length(test2)>0) {
    creactions[[i]] = cbind(rep(testmetNames[[i]], times = length(test2)), ReactionData[test2,], delta[test2,])
    colnames(creactions[[i]]) = c("MetNames",colnames(creactions[[i]])[2:3], "Flux")
    testC = aggregate(Flux~Subsystems, data = creactions[[i]], FUN = sum)
    MetC[match(testC$Subsystems, Data$SubSytems),i] = testC$Flux
  }
}


colnames(MetP) = testmetNames
colnames(MetC) = testmetNames


rownames(MetP) = Data$SubSytems
rownames(MetC) = Data$SubSytems

cat = do.call("rbind", creactions)
prod = do.call("rbind",preactions)
# MetP_del[which(MetP_del==0)] = ""
# MetC_del[which(MetC_del==0)] = ""
#write.table(cat,"AnchoredTrp_MinEx_Mets_Cons.txt", row.names = F, col.names = T, sep = "\t")
#write.table(prod,"AnchoredTrp_MinEx_Mets_Prod.txt", row.names = F, col.names = T, sep = "\t")

library(plotly)
library(RColorBrewer)
library(BBmisc)
library(ggplot2)
library(hrbrthemes)
library(plotly)
library(reshape2)
library(plyr)
test = data.frame(Met = testmetNames, Value = colSums(MetP), N = np, stringsAsFactors = F)
p_fil = test[order(-abs(test$Value)),]


test = data.frame(Met = testmetNames, Value = colSums(MetC), N = nc, stringsAsFactors = F)
c_fil = test[order(-abs(test$Value)),]


rm(list= ls()[!(ls() %in% c('p_fil','c_fil'))])
graphics.off()
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

setwd(paste0(here(),"/Simulations/Myocyte//"))
delta = read.xlsx("Results_25_1_GSE25462.xlsx", 1, header = F)

setwd(paste0(here(),"/Simulations/Myocyte/"))
ReactionData = readMat("model.mat")
ReactionData = lapply(ReactionData$model[,,1], unlist)

S = as.matrix(ReactionData$S)
mets =  ReactionData$mets
Names = ReactionData$metNames
Comp = as.numeric(ReactionData$metComps)

Compartments = data.frame(Names = ReactionData$compNames, Abbr = ReactionData$comps, stringsAsFactors = F)
Comp = as.character(Compartments$Abbr[Comp])

MetNames = paste(Names, "[", Comp, "]", sep ="")

ReactionData = data.frame(ID = ReactionData$rxns, Subsystems = ReactionData$subSystems, stringsAsFactors = F)
rownames(S) = MetNames
colnames(S) = ReactionData$ID

Data = as.data.frame(table(ReactionData$Subsystems))
Data$Var1 = as.character(Data$Var1)
colnames(Data) =  c("SubSytems","Number of DE reactions")


testmetNames = unique(Names)
nc = c()
np = c()
MetP = matrix(0, nrow = length(Data$SubSytems), ncol = length(testmetNames))
MetC = matrix(0, nrow = length(Data$SubSytems), ncol = length(testmetNames))
creactions = list()
preactions = list()


for (i in 1:length(testmetNames)) {
  allIdx = which(Names==testmetNames[[i]])
  test = as.matrix(S[allIdx,])
  allIdx1 = as.numeric(unlist(lapply(1:nrow(test), function(x) which(test[x,]>0))))
  allIdx2 = as.numeric(unlist(lapply(1:nrow(test), function(x) which(test[x,]<0))))
  nc[i] = length(allIdx2)
  np[i] = length(allIdx1)
  test1 = setdiff(allIdx1, allIdx2)
  
  if (length(test1)>0) {
    preactions[[i]] = cbind(rep(testmetNames[[i]], times = length(test1)), ReactionData[test1,], delta[test1,])
    colnames(preactions[[i]]) = c("MetNames",colnames(preactions[[i]])[2:3], "Flux")
    testP = aggregate(Flux~Subsystems, data = preactions[[i]], FUN = sum)
    MetP[match(testP$Subsystems, Data$SubSytems),i] = testP$Flux
  }
  
  
  test2 = setdiff(allIdx2, allIdx1)
  
  if (length(test2)>0) {
    creactions[[i]] = cbind(rep(testmetNames[[i]], times = length(test2)), ReactionData[test2,], delta[test2,])
    colnames(creactions[[i]]) = c("MetNames",colnames(creactions[[i]])[2:3], "Flux")
    testC = aggregate(Flux~Subsystems, data = creactions[[i]], FUN = sum)
    MetC[match(testC$Subsystems, Data$SubSytems),i] = testC$Flux
  }
}


colnames(MetP) = testmetNames
colnames(MetC) = testmetNames


rownames(MetP) = Data$SubSytems
rownames(MetC) = Data$SubSytems

cat = do.call("rbind", creactions)
prod = do.call("rbind",preactions)
# MetP_del[which(MetP_del==0)] = ""
# MetC_del[which(MetC_del==0)] = ""
#write.table(cat,"AnchoredTrp_MinEx_Mets_Cons.txt", row.names = F, col.names = T, sep = "\t")
#write.table(prod,"AnchoredTrp_MinEx_Mets_Prod.txt", row.names = F, col.names = T, sep = "\t")

library(plotly)
library(RColorBrewer)
library(BBmisc)
library(ggplot2)
library(hrbrthemes)
library(plotly)
library(reshape2)
library(ggrepel)

test = data.frame(Met = testmetNames, Value = colSums(MetP), N = np, stringsAsFactors = F)
p_fil2 = test[order(-abs(test$Value)),]


test = data.frame(Met = testmetNames, Value = colSums(MetC), N = nc, stringsAsFactors = F)
c_fil2 = test[order(-abs(test$Value)),]


p_fil$Study = "van Tienen et al."
c_fil$Study = "van Tienen et al."
p_fil2$Study = "Jin et al."
c_fil2$Study = "Jin et al."

p_fil = rbind(p_fil , p_fil2)
c_fil = rbind(c_fil, c_fil2)

rm(list= ls()[!(ls() %in% c('p_fil','c_fil'))])
cf = c_fil
pf = p_fil

p_fil = p_fil[p_fil$N>12,]
p_fil = p_fil[abs(p_fil$Value)>10,]
p_fil$Reg = paste0("Positive", "  \u0394\U1D4E5")
p_fil$Reg[p_fil$Value<0] = paste0("Negative", "  \u0394\U1D4E5")


setwd(paste0(here(),"/Plotting/Myocyte/"))
lim = round_any(max(c(abs(min(p_fil$Value)), max(p_fil$Value))),100,ceiling)
#tiff('T2D_25_1_P.tiff', units="px", width=(3600), height=(3600), res=300)
p = ggplot(data = p_fil, aes(x=N, y=Value, shape = Study, color = Reg)) +
  
  geom_point(alpha=0.65, size = 3.5) +
  geom_text_repel(data=p_fil, aes(x=N, y=Value, label=Met, color = factor(Reg)),arrow = arrow(
    length = unit(0.01, "npc"),  ends = "first"
  ),family = "Calibri", fontface="bold",alpha=0.9, size=6,force = 10,inherit.aes = FALSE, show.legend = F ) +
  scale_color_manual(values = c("#008000", "#ff0000"), name = "Regulation")+
  
  theme_bw(base_size = 12) +
  geom_segment(aes(x = 0, y = 0, xend = 500, yend = 0), colour = "black", alpha=0.5, size=0.25 , inherit.aes = FALSE) +
  xlim(0,500) + 
  geom_segment(aes(x = 0, y = lim, xend = 0, yend = -lim), colour = "black", alpha=0.5, size=0.25 , inherit.aes = FALSE) +
  ylim(-lim,lim) +
  theme(panel.border = element_rect(color = "#476b6b", fill = NA, size = 1.5), 
        strip.background = element_rect(color = "#476b6b", size = 1.5, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  xlab("Number of Reactions") + ylab("Production Flux Throughput") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold")) +
  theme(legend.text=element_text(size=16, family = "Calibri", face = "bold")) +
  theme(legend.title = element_text(size = 18, family = "Calibri", face = "bold")) +
  theme(legend.position="bottom") 
p
#dev.off()

p_fil$v = 0
p_fil$v[which(p_fil$Value>0)] = log10(p_fil$Value[which(p_fil$Value>0)])
p_fil$v[which(p_fil$Value<0)] = log10(abs(p_fil$Value[which(p_fil$Value<0)]))

#tiff('T2D_25_1_P2.tiff', units="px", width=(3600), height=(3600), res=300)
p = ggplot(data = p_fil, aes(x=N, y=v, shape = Study, color = Reg)) +
  
  geom_point(alpha=0.7, size = 4) +
  geom_text_repel(data=p_fil, aes(x=N, y=v, label=Met, color = factor(Reg)),arrow = arrow(
    length = unit(0.01, "npc"),  ends = "first"
  ),family = "Calibri", fontface="bold",alpha=0.9, size=6,force = 10,inherit.aes = FALSE, show.legend = F ) +
  scale_color_manual(values = c("#008000", "#ff0000"), name = "Regulation")+
  
  theme_bw(base_size = 12) +
  geom_segment(aes(x = 0, y = 1, xend = 500, yend = 1), colour = "black", alpha=0.5, size=0.25 , inherit.aes = FALSE) +
  xlim(0,500) + 
  geom_segment(aes(x = 0, y = 3.25, xend = 0, yend = 1), colour = "black", alpha=0.5, size=0.25 , inherit.aes = FALSE) +
  
  theme(panel.border = element_rect(color = "#476b6b", fill = NA, size = 1.5), 
        strip.background = element_rect(color = "#476b6b", size = 1.5, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  xlab("Number of Production Reactions") + ylab("Log(10) Predicted Flux Difference") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold")) +
  theme(legend.text=element_text(size=16, family = "Calibri", face = "bold")) +
  theme(legend.title = element_text(size = 18, family = "Calibri", face = "bold")) +
  theme(legend.position="bottom") 
p
#dev.off()


p_fil  = pf
p_fil = p_fil[p_fil$N>12,]
p_fil = p_fil[abs(p_fil$Value)>10,]
p_fil$Reg = paste0("Positive", "  \u0394\U1D4E5")
p_fil$Reg[p_fil$Value<0] = paste0("Negative", "  \u0394\U1D4E5")
p_fil$v = 0
p_fil$v[which(p_fil$Value>0)] = log10(p_fil$Value[which(p_fil$Value>0)])
p_fil$v[which(p_fil$Value<0)] = log10(abs(p_fil$Value[which(p_fil$Value<0)]))
p_fil = p_fil[order(-p_fil$Value),]
rownames(p_fil) = c(1:length(p_fil$Met))
p_fil = p_fil[-c(4,5,9,10,15,16,26,34,48,51),]
add = c()
cond = unique(p_fil$Study)
for (i in 1:length(p_fil$Met)) {
  k = length(which(p_fil$Met==p_fil$Met[i]))
  if (k==1) {
    t = p_fil[i,]
    t$Value = 0
    t$v = 0
    t$Study = setdiff(cond, t$Study)
    add = rbind(add, t)
  }
}
final = rbind(p_fil,add)
final = final[order(-final$Value),]
final$Met = factor(final$Met, levels = unique(final$Met))

tiff('T2D_25_1_P3.tiff', units="px", width=(3900), height=(2250), res=300)

p = ggplot(final, aes(x = Met, y = v, fill = Study)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.9) +
  scale_fill_manual(values = c("#0b84a5", "#f6c85f")) +
  theme_bw(base_size = 12) +
  theme(panel.border = element_rect(color = "#476b6b", fill = NA, size = 1.5), 
        strip.background = element_rect(color = "#476b6b", size = 1.5, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  xlab("") + ylab("Log(10) Predicted Flux Difference") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold",angle = 90, hjust = 1, vjust = 0.4)) +
  theme(legend.text=element_text(size=16, family = "Calibri", face = "bold")) +
  theme(legend.title = element_text(size = 18, family = "Calibri", face = "bold")) +
  theme(legend.position="bottom") 
p
dev.off()

p_fil  = pf
p_fil = p_fil[abs(p_fil$Value)>10,]
p_fil$Reg = paste0("Positive", "  \u0394\U1D4E5")
p_fil$Reg[p_fil$Value<0] = paste0("Negative", "  \u0394\U1D4E5")
p_fil$v = 0
p_fil$v[which(p_fil$Value>0)] = log10(p_fil$Value[which(p_fil$Value>0)])
p_fil$v[which(p_fil$Value<0)] = log10(abs(p_fil$Value[which(p_fil$Value<0)]))
p_fil = p_fil[order(-p_fil$Value),]
p_fil$r = round(p_fil$Value, digits = 1)
add = c()
for (i in 1:length(p_fil$Met)) {
  k = length(intersect(which(p_fil$r==p_fil$r[i]), which(p_fil$Study==p_fil$Study[i])))
  if (k>1) {
    add = c(add,i)
  }
}
p_fil = p_fil[-add,]
rownames(p_fil) = c(1:length(p_fil$Met))
p_fil = p_fil[-c(21,23,25,26,38,56,70,81,86),]
add = c()
cond = unique(p_fil$Study)
for (i in 1:length(p_fil$Met)) {
  k = length(which(p_fil$Met==p_fil$Met[i]))
  if (k==1) {
    t = p_fil[i,]
    t$Value = 0
    t$v = 0
    t$Study = setdiff(cond, t$Study)
    add = rbind(add, t)
  }
}
final = rbind(p_fil,add)
final = final[order(-final$Value),]
final$Met = factor(final$Met, levels = unique(final$Met))

#tiff('T2D_25_1_P4.tiff', units="px", width=(1800*3), height=(1000*3), res=300)

p = ggplot(final, aes(x = Met, y = v, fill = Study)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.9) +
  scale_fill_manual(values = c("#0b84a5", "#f6c85f")) +
  theme_bw(base_size = 12) +
  theme(panel.border = element_rect(color = "#476b6b", fill = NA, size = 1.5), 
        strip.background = element_rect(color = "#476b6b", size = 1.5, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(size=1,color='#476b6b')) +
  xlab("") + ylab("Log(10) Predicted Flux Difference") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 20, family = "Calibri", face = "bold")) +
  theme(axis.text.y = element_text(size=18, family = "Calibri", face = "bold"), axis.text.x = element_text(size = 18, family = "Calibri", face = "bold",angle = 90, hjust = 1, vjust = 0.4)) +
  theme(legend.text=element_text(size=16, family = "Calibri", face = "bold")) +
  theme(legend.title = element_text(size = 18, family = "Calibri", face = "bold")) +
  theme(legend.position="bottom") 
p
#dev.off()

