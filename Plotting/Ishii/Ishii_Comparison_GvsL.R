rm(list = ls(all=T))
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

setwd("D:/GitHub/D_FBA/Data/Ishii/")
wd = read.xlsx("Fluxomics.xlsx", 1)
conditions = colnames(wd)[c(9,10,15,21,24,5,6)]
conditions[6] = "D = 0.5/h"
conditions[7] = "D = 0.7/h"
fluxes = as.character(wd[,1])
wtflux = wd[,2]
M = wd[,c(9,10,15,21,24,5,6)]

setwd("D:/GitHub/D_FBA/Simulations/Ishii/Array_Relaxed_DFBA//")
Md = read.xlsx("Results.xlsx",1, header = F)
Pd = read.xlsx("Results.xlsx",2, header = F)
Rd = read.xlsx("Results.xlsx",3, header = F)
MRatio = M/wtflux
MRatio = as.matrix(MRatio)
MRatio[is.nan(MRatio)] = 0
MRatio[which(MRatio==Inf)] = 100
MRatio[which(MRatio==-Inf)] = 0.001


Cor = c()
d_cor = c()
d_err = c()
accuracy = c()
r.sq = c()
r.sqratio = c()

for (i in 1:length(conditions)) {
  P = Pd[,i] + wtflux
  sel = intersect(which(P!=0), which(M[,i]!=0))
  Cor[i] =  P[sel] %*%  M[sel,i] / (norm(as.matrix(P[sel]), type = "2") * norm(as.matrix(M[sel,i]), type = "2"))
  sel = intersect(which(Pd[,i]!=0), which(Md[,i]!=0))
  d_cor[i] =  Pd[sel,i] %*%  Md[sel,i] / (norm(as.matrix(Pd[sel,i]), type = "2") * norm(as.matrix(Md[sel,i]), type = "2"))
  d_err[i] = tdStats(as.numeric(Pd[,i]), as.numeric(Md[,i]), functions= c("nrmse"))
  s_Md = as.numeric(Md[,i])
  s_Md[s_Md>0] = 1;
  s_Md[s_Md<0] = -1;
  
  s_Pd = as.numeric(Pd[,i])
  s_Pd[s_Pd>0] = 1;
  s_Pd[s_Pd<0] = -1;
  accuracy[i] = accuracy(s_Pd, s_Md)
  
  r.sq[i] = summary(lm(Rd[,i]~Md[,i]))$r.squared
  r.sqratio[i] = summary(lm(Rd[,i]~MRatio[,i]))$r.squared
  
}

p_cond = conditions
p_dcor1 = d_cor
p_acc1 = accuracy
data <- data.frame(
  individual=rep(p_cond,2),
  group=c(rep('Acc', length(p_cond)), rep('dpCC', length(p_cond))) ,
  value= c(p_acc1, p_dcor1))
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

setwd("D:/GitHub/D_FBA/Data/Ishii/")
wd = read.xlsx("Fluxomics.xlsx", 1)
conditions = colnames(wd)[c(9,10,15,21,24,5,6)]
conditions[6] = "D = 0.5/h"
conditions[7] = "D = 0.7/h"
fluxes = as.character(wd[,1])
wtflux = wd[,2]
M = wd[,c(9,10,15,21,24,5,6)]

setwd("D:/GitHub/D_FBA/Simulations/Ishii/RT_PCR_Relaxed_DFBA//")
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


Cor = c()
d_cor = c()
d_err = c()
accuracy = c()
r.sq = c()
r.sqratio = c()

for (i in 1:length(conditions)) {
  P = Pd[,i] + wtflux
  sel = intersect(which(P!=0), which(M[,i]!=0))
  Cor[i] =  P[sel] %*%  M[sel,i] / (norm(as.matrix(P[sel]), type = "2") * norm(as.matrix(M[sel,i]), type = "2"))
  sel = intersect(which(Pd[,i]!=0), which(Md[,i]!=0))
  d_cor[i] =  Pd[sel,i] %*%  Md[sel,i] / (norm(as.matrix(Pd[sel,i]), type = "2") * norm(as.matrix(Md[sel,i]), type = "2"))
  d_err[i] = tdStats(as.numeric(Pd[,i]), as.numeric(Md[,i]), functions= c("nrmse"))
  s_Md = as.numeric(Md[,i])
  s_Md[s_Md>0] = 1;
  s_Md[s_Md<0] = -1;
  
  s_Pd = as.numeric(Pd[,i])
  s_Pd[s_Pd>0] = 1;
  s_Pd[s_Pd<0] = -1;
  accuracy[i] = accuracy(s_Pd, s_Md)
  
  r.sq[i] = summary(lm(Rd[,i]~Md[,i]))$r.squared
  r.sqratio[i] = summary(lm(Rd[,i]~MRatio[,i]))$r.squared
  
}

p_cond = conditions
p_dcor2 = d_cor
p_acc2 = accuracy
data2 <- data.frame(
  individual=rep(p_cond,2),
  group=c(rep('Acc', length(p_cond)), rep('dpCC', length(p_cond))) ,
  value= c(p_acc2, p_dcor2))
data2$test = "Local"
rm(list= ls()[!(ls() %in% c('data', 'data2', 'p_cond'))])
all = c(mean(data$value[data$group=="Acc"]), mean(data2$value[data2$group=="Acc"]), mean(data$value[data$group=="dpCC"]),mean(data2$value[data2$group=="dpCC"]))
s1 = all

data$group = factor(data$group, levels = c("dpCC", "Acc"))
data2$group = factor(data2$group, levels = c("dpCC", "Acc"))
data = data[order(data$group),]
data2 = data2[order(data2$group),]

data$id = seq(1, nrow(data))
data2$id = seq(1, nrow(data2))
data = rbind(data,data2)
data = data[order(data$id),]


# Set a number of 'empty bar' to add at the end of each group (Spacing between groups)
empty_bar <- 1
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)




data$id[15] = 8
data$id[data$group=="Acc"] = data$id[data$group=="Acc"]+1
data$id[length(data$individual)] = 16

data$newgroup = paste(as.character(data$group), "-", as.character(data$test), sep = "")
# Get the name and the y position of each label
label_data <- data
label_data$value[which(label_data$value<0)] = 0
number_of_bar <- nrow(label_data)
label_data$labelid = seq(1, nrow(label_data))
angle <- 90 - 360 * (label_data$labelid-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

keys <- colnames(label_data)[c(1,2,5)]
X <- as.data.table(label_data)
label_data = as.data.frame(X[,list(angle= mean(angle), hjust = mean(hjust), value = max(value)),keys])
label_data = label_data[complete.cases(label_data),]


#label_data$value[label_data$value>1] = 0.9

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))


# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end + 1
grid_data$start <- grid_data$start - 1 



## Adjust
base_data$start = base_data$start-0.5
base_data$end = base_data$end+0.5
base_data$title = rowMeans(cbind(base_data$start,base_data$end))

label_data$angle[1] = 65
label_data$angle[2] = 45

grid_data$start = grid_data$start+0.5
grid_data$end = grid_data$end-0.5


# Make the plot
loadfonts(device = "win")
setwd("D:/GitHub/D_FBA/Plotting/Ishii/")
par(family = "Calibri", font = 2)
tiff('Ishii_Comparison_GlobalvsLocal.tiff', units="px", width=(1200*3), height=(1000*3), res=300)
p <- ggplot(data, aes(x=as.factor(id), y=value, fill = newgroup)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(aes(x=as.factor(id), y=value, fill = newgroup),position = "dodge", stat="identity", alpha=0.75) +
  #scale_fill_manual(values = c("#e15759", "white", "#f0abac", "#4379a7","white", "#a1bcd3")) +
  scale_fill_manual(values = c("#2d506f", "#4983b5", "white", "#2d506f", "#4983b5", "white")) +
  #scale_fill_viridis(alpha = 0.8, begin = 0, end = 1, direction = 1,discrete =T, option = "D") +
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 1, xend = start, yend = 1), colour = "grey", alpha=0.5, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.5, xend = start, yend = 0.5), colour = "grey", alpha=0.5, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = -0.5, xend = start, yend = -0.5), colour = "grey", alpha=0.5, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "black", alpha=0.5, size=1 , inherit.aes = FALSE ) +
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(0,3), y = c(-0.5,0.5,1), label = c("-.5",".5", "1") , color="black", size= 3.75 , angle=0, fontface="bold", hjust=0.5) +
  annotate("text", x = rep(0,1), y = c(0), label = c("0") , color="black", size= 3.75 , angle=0, fontface="bold", hjust=0.5) +
  geom_bar(aes(x=as.factor(id), y=value, fill = newgroup),position = "dodge", stat="identity", alpha=0.75) +
  ylim(-1.250,1.50) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+0.1, label=individual, hjust=hjust), color="black",family = "Calibri", fontface="bold",alpha=0.9, size=4.75, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_rect(data=base_data, aes(xmin = start, ymin = -0.75, xmax = end, ymax = 1.42), fill = "lightsalmon", color= "darkolivegreen", alpha=0.1, size=0.5 , inherit.aes = FALSE )  
#geom_text(data=base_data, aes(x = title, y = 1.45, label=c("bold(rho)","bold(Accuracy)")), hjust=c(0,1), vjust=c(0.5,0.5), colour = "darkolivegreen", alpha=0.9, size=8, family = "Calibri", fontface="bold", inherit.aes = FALSE,parse = T) 


p
dev.off()

tiff('Ishii_Comparison_GlobalvsLocalC.tiff', units="px", width=(900*3), height=(900*3), res=300)
p <- ggplot(data, aes(x=as.factor(id), y=value, fill = newgroup)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(aes(x=as.factor(id), y=value, fill = newgroup),position = "dodge", stat="identity", alpha=0.75) +
  #scale_fill_manual(values = c("#e15759", "white", "#f0abac", "#4379a7","white", "#a1bcd3")) +
  scale_fill_manual(values = c("#2d506f", "#4983b5", "white", "#2d506f", "#4983b5", "white")) +
  #scale_fill_viridis(alpha = 0.8, begin = 0, end = 1, direction = 1,discrete =T, option = "D") +
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 1, xend = start, yend = 1), colour = "grey", alpha=0.5, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.5, xend = start, yend = 0.5), colour = "grey", alpha=0.5, size=0.5 , inherit.aes = FALSE ) +
  #geom_segment(data=grid_data, aes(x = end, y = -0.5, xend = start, yend = -0.5), colour = "grey", alpha=0.5, size=0.5 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "black", alpha=0.5, size=1 , inherit.aes = FALSE ) +
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(0,2), y = c(0.5,1), label = c(".5", "1") , color="black", size= 3.75 , angle=0, fontface="bold", hjust=0.5) +
  annotate("text", x = rep(0,1), y = c(0), label = c("0") , color="black", size= 3.75 , angle=0, fontface="bold", hjust=0.5) +
  geom_bar(aes(x=as.factor(id), y=value, fill = newgroup),position = "dodge", stat="identity", alpha=0.75) +
  ylim(-0.50,1.30) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+0.1, label=individual, hjust=hjust), color="black",family = "Calibri", fontface="bold",alpha=0.9, size=4.75, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_rect(data=base_data, aes(xmin = start, ymin = -0.25, xmax = end, ymax = 1.22), fill = "lightsalmon", color= "darkolivegreen", alpha=0.1, size=0.5 , inherit.aes = FALSE )  
#geom_text(data=base_data, aes(x = title, y = 1.45, label=c("bold(rho)","bold(Accuracy)")), hjust=c(0,1), vjust=c(0.5,0.5), colour = "darkolivegreen", alpha=0.9, size=8, family = "Calibri", fontface="bold", inherit.aes = FALSE,parse = T) 


p
dev.off()

rm(list= ls()[!(ls() %in% c('s1'))])

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

setwd("D:/GitHub/D_FBA/Data/Ishii/")
wd = read.xlsx("Fluxomics.xlsx", 1)
conditions = colnames(wd)[c(9,10,15,21,24,5,6)]
conditions[6] = "D = 0.5/h"
conditions[7] = "D = 0.7/h"
fluxes = as.character(wd[,1])
wtflux = wd[,2]
M = wd[,c(9,10,15,21,24,5,6)]

setwd("D:/GitHub/D_FBA/Simulations/Ishii/Array_Relaxed_DFBA/")
Md = read.xlsx("Results.xlsx",1, header = F)
Pd = read.xlsx("Results.xlsx",2, header = F)
Rd = read.xlsx("Results.xlsx",3, header = F)
MRatio = M/wtflux
MRatio = as.matrix(MRatio)
MRatio[is.nan(MRatio)] = 0
MRatio[which(MRatio==Inf)] = 100
MRatio[which(MRatio==-Inf)] = 0.001


Cor = c()
d_cor = c()
d_err = c()
accuracy = c()
r.sq = c()
r.sqratio = c()

for (i in 1:length(conditions)) {
  P = Pd[,i] + wtflux
  sel = intersect(which(P!=0), which(M[,i]!=0))
  Cor[i] =  P[sel] %*%  M[sel,i] / (norm(as.matrix(P[sel]), type = "2") * norm(as.matrix(M[sel,i]), type = "2"))
  sel = intersect(which(Pd[,i]!=0), which(Md[,i]!=0))
  d_cor[i] =  Pd[sel,i] %*%  Md[sel,i] / (norm(as.matrix(Pd[sel,i]), type = "2") * norm(as.matrix(Md[sel,i]), type = "2"))
  d_err[i] = tdStats(as.numeric(Pd[,i]), as.numeric(Md[,i]), functions= c("nrmse"))
  s_Md = as.numeric(Md[,i])
  s_Md[s_Md>0] = 1;
  s_Md[s_Md<0] = -1;
  
  s_Pd = as.numeric(Pd[,i])
  s_Pd[s_Pd>0] = 1;
  s_Pd[s_Pd<0] = -1;
  accuracy[i] = accuracy(s_Pd, s_Md)
  
  r.sq[i] = summary(lm(Rd[,i]~Md[,i]))$r.squared
  r.sqratio[i] = summary(lm(Rd[,i]~MRatio[,i]))$r.squared
  
}

p_cond = conditions
p_dcor1 = d_err
p_acc1 = r.sqratio
data <- data.frame(
  individual=rep(p_cond,2),
  group=c(rep('R^2', length(p_cond)), rep('NRMSE', length(p_cond))) ,
  value= c(p_acc1, p_dcor1))
data$test = "Global"

rm(list= ls()[!(ls() %in% c('data','s1'))])
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

setwd("D:/GitHub/D_FBA/Data/Ishii/")
wd = read.xlsx("Fluxomics.xlsx", 1)
conditions = colnames(wd)[c(9,10,15,21,24,5,6)]
conditions[6] = "D = 0.5/h"
conditions[7] = "D = 0.7/h"
fluxes = as.character(wd[,1])
wtflux = wd[,2]
M = wd[,c(9,10,15,21,24,5,6)]

setwd("D:/GitHub/D_FBA/Simulations/Ishii/RT_PCR_Relaxed_DFBA/")
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


Cor = c()
d_cor = c()
d_err = c()
accuracy = c()
r.sq = c()
r.sqratio = c()

for (i in 1:length(conditions)) {
  P = Pd[,i] + wtflux
  sel = intersect(which(P!=0), which(M[,i]!=0))
  Cor[i] =  P[sel] %*%  M[sel,i] / (norm(as.matrix(P[sel]), type = "2") * norm(as.matrix(M[sel,i]), type = "2"))
  sel = intersect(which(Pd[,i]!=0), which(Md[,i]!=0))
  d_cor[i] =  Pd[sel,i] %*%  Md[sel,i] / (norm(as.matrix(Pd[sel,i]), type = "2") * norm(as.matrix(Md[sel,i]), type = "2"))
  d_err[i] = tdStats(as.numeric(Pd[,i]), as.numeric(Md[,i]), functions= c("nrmse"))
  s_Md = as.numeric(Md[,i])
  s_Md[s_Md>0] = 1;
  s_Md[s_Md<0] = -1;
  
  s_Pd = as.numeric(Pd[,i])
  s_Pd[s_Pd>0] = 1;
  s_Pd[s_Pd<0] = -1;
  accuracy[i] = accuracy(s_Pd, s_Md)
  
  r.sq[i] = summary(lm(Rd[,i]~Md[,i]))$r.squared
  r.sqratio[i] = summary(lm(Rd[,i]~MRatio[,i]))$r.squared
  
}

p_cond = conditions
p_dcor2 = d_err
p_acc2 = r.sqratio
data2 <- data.frame(
  individual=rep(p_cond,2),
  group=c(rep('R^2', length(p_cond)), rep('NRMSE', length(p_cond))) ,
  value= c(p_acc2, p_dcor2))
data2$test = "Local"
rm(list= ls()[!(ls() %in% c('data', 'data2', 'p_cond','s1'))])
all = c(mean(data$value[data$group=="NRMSE"]), mean(data2$value[data2$group=="NRMSE"]))
s2 = all

data$group = factor(data$group, levels = c("NRMSE", "R^2"))
data2$group = factor(data2$group, levels = c("NRMSE", "R^2"))
data = data[order(data$group),]
data2 = data2[order(data2$group),]

data$id = seq(1, nrow(data))
data2$id = seq(1, nrow(data2))
data = rbind(data,data2)
data = data[order(data$id),]

library(fmsb)

t1 = data[data$group=="NRMSE",]
t1 = t1[complete.cases(t1),]

t2 = matrix(0, nrow = 2, ncol = length(p_cond))
t2[1,] = t1$value[t1$test=="Global"]
t2[2,] = t1$value[t1$test=="Local"]
rownames(t2) = c("Global", "Local")
colnames(t2) = p_cond
t2 = as.data.frame(t2)
# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
t2 <- rbind(rep(0.25,length(p_cond)) , rep(0,length(p_cond)) , t2)

# Color vector
colors_border=c("#2d506f", "#4983b5")
colors_in=c(rgb(0.6,0.2,0.2,0.3),rgb(0.2,0.3,0.5,0.3))
par(family = "Calibri", font = 2)

setwd("D:/GitHub/D_FBA/Plotting/Ishii/")
par(family = "Calibri", font = 2)
tiff('Ishii_Comparison_GlobalvsLocal_B.tiff', units="px", width=(1000*3), height=(1000*3), res=300)
par(family = "Calibri", font = 2)
# plot with default options:
radarchart(t2 , axistype=1 , 
           #custom polygon
           pcol=colors_border ,pfcol=colors_in , plwd=4 , plty=1, seg = 5,pty = 32,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(0,0.25,0.05), cglwd=0.8,
           #custom labels
           vlcex=1.35, calcex = 1.35
)
dev.off()

s1
s2

s3 =  matrix(s1, 2,2)
s3 = cbind(s3, s2)
colnames(s3) = c("Accuracy","Rho", "NRMSE")
rownames(s3) = c("Global", "Local")

write.table(s3,"Ishii_Comparison_GvsL.txt", row.names = T, col.names = T, sep = "\t")