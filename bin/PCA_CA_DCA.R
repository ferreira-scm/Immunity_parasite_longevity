# PCA and CA and DCA
# data analysis
library("vegan")
library(MASS)
library(ggplot2)
setwd("/home/susana/Documents/HyenaMan/ImmuneHyena")
complete=read.csv("data/hyenadataset.csv")

#Hypothesis: adults have distinct parasite and immune profiles than juveniles
# Ordination in reduced space
# Detrended Correspondance Analysis - DCA
# 18 individuals do not have any parasites. We need to remove those for the
# next steps.
nrow(na.omit(complete[,c(6:13)])[rowSums(na.omit(complete[,c(6:13)]))==0,])
# I consider that rare taxa are informative, so I don't use
#the option to down-weight them: "iweigh=1"
Pdca=decorana(na.omit(complete[,c(6:13)])[rowSums(na.omit(complete[,c(6:13)]))>0,])
#correspondence analysis
Pca=decorana(na.omit(complete[,c(6:13)])[rowSums(na.omit(complete[,c(6:13)]))>0,], ira=1)
plot(Pca)

# let's check axis length
summary(Pdca)
# sample(species) scores for the first axis
PDC1= scores(Pdca, display = c("sites"), choices = c(1))
PC1= scores(Pca, display = c("sites"), choices = c(1))
#Immune DCA
Idca=decorana(na.omit(complete[,c(17:19)]))
Ica=decorana(na.omit(complete[,c(17:19)]), ira=1)
plot(Idca)
plot(Pdca)
# sample(species) scores for the first axis
IDC1= scores(Idca, display = c("sites"), choices = c(1))
IC1= scores(Ica, display = c("sites"), choices = c(1))
summary(Idca)

# immune and parasite PCA
#first pca follwed by randomization test
Ppca = prcomp((na.omit(complete[,c(6:13)])), scale. = T, center = T)
Ipca = prcomp((na.omit(complete[,c(17:19)])), scale. = T, center = T)
sign.pc<-function(x,R=1000,s=ncol(x)){
  # run PCA
  pc.out<-prcomp(x,scale.=T, center=T)
  # the proportion of variance of each PC
  pve=(pc.out$sdev^2/sum(pc.out$sdev^2))[1:s]
  # a matrix with R rows and s columns that contains
  # the proportion of variance explained by each pc
  # for each randomization replicate.
  pve.perm<-matrix(NA,ncol=s,nrow=R)
  for(i in 1:R){
    # permutation each column
    x.perm<-apply(x,2,sample)
    # run PCA
    pc.perm.out<-prcomp(x.perm, scale. = T, center = T)
    # the proportion of variance of each PC.perm
    pve.perm[i,]=(pc.perm.out$sdev^2/sum(pc.perm.out$sdev^2))[1:s]
  }
  # calcalute the p-values
  pval<-apply(t(pve.perm)>pve,1,sum)/R
  return(list(pve=pve,pval=pval))
}
sign.pc(na.omit(complete[,c(6:13)]))
sign.pc(na.omit(complete[,c(17:19)]))

#extracting scores from first axis and adding to dataset
PDC1<- as.data.frame(PDC1)
PDC1$row <- rownames(PDC1)
complete$row <- rownames(complete)
complete <- merge(complete, PDC1, by = "row", all = T)
PC1<- as.data.frame(PC1)
PC1$row <- rownames(PC1)
complete <- merge(complete, PC1, by = "row", all = T)

complete$PDC1
IDC1<- as.data.frame(IDC1)
IDC1$row <- rownames(IDC1)
complete <- merge(complete, IDC1, by = "row", all = T)
IC1<- as.data.frame(IC1)
IC1$row <- rownames(IC1)
complete <- merge(complete, IC1, by = "row", all = T)


complete$IC1


PPC2 <- Ppca$x[,2]
PPC2<- as.data.frame(PPC2)
PPC2$row <- rownames(PPC2)
complete <- merge(complete, PPC2, by = "row", all = T)
complete$PPC2
IPC1 <- Ipca$x[,1]
IPC1<- as.data.frame(IPC1)
IPC1$row <- rownames(IPC1)
complete <- merge(complete, IPC1, by = "row", all=T)
complete$IPC1
complete <- complete[, -1]



### anosim
## ANOSIM
dataP$Sex=droplevels(dataP$Sex)
dataP <- complete[complete.cases(complete$PAxis1),]
dataI <- complete[complete.cases(complete$IAxis1),]

table(dataI$Age)

P.dist = vegdist(na.omit(complete[,c(6:13)])[rowSums(na.omit(complete[,c(6:13)]))>0,])
P.ano_Agesex=anosim(P.dist,dataP$Age_sex, permutations = 99999, distance = "bray")
summary(P.ano_Agesex)
P.ano_Age=anosim(P.dist,dataP$Age, permutations = 10000, distance = "bray")
P.ano_Sex=anosim(P.dist,dataP$Sex, permutations = 10000, distance = "bray")
summary(P.ano_Agesex)
I.dist = vegdist(na.omit(complete[,c(17:19)]))
I.ano.AgeSex=anosim(I.dist,dataI$Age_sex, permutations = 10000, distance = "bray")
I.ano.Age=anosim(I.dist,dataI$Age, permutations = 10000, distance = "bray")
I.ano.Sex=anosim(I.dist,dataI$Sex, permutations = 10000, distance = "bray")

saveRDS(P.ano_Agesex, "tmp/P.ano_AgeSex")
P.ano_Agesex=readRDS("tmp/P.ano_AgeSex")
P.ano_sex=readRDS("tmp/P.ano_sex")
P.ano_age=readRDS("tmp/P.ano.Age")
I.ano_sex=readRDS("tmp/I.ano.Sex")
I.ano_age=readRDS("tmp/I.ano.Age")



Imodel <- with(data, lm(IPC1 ~ PPC2 + Age + Sex + PPC2:Age + Sex:Age))
Imodel <- with(data, lm(IDC1 ~ PDC1 + Age + Sex + PDC1:Age + Sex:Age))
summary(Imodel)
Imodel_P <- with(data, lm(((IDC1) ~ Age + Sex + Sex:Age)))
Imodel_Age <- with(data, lm(((IDC1) ~ PDC1 + Sex )))
Imodel_Sex <- with(data, lm(((IDC1) ~ PDC1 + Age + PDC1:Age)))
Imodel_PA <- with(data, lm(((IDC1) ~ PDC1 + Age + Sex + Sex:Age)))
Imodel_SA <- with(data, lm(((IDC1) ~ PDC1 + Age + Sex + PDC1:Age)))
Imodel_int <- with(data, lm(((IDC1) ~ 1)))
# plot(Imodel)
anova(Imodel, test = "LRT")
anova(Imodel, test = "Chisq")
