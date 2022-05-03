# data analysis
#setwd("/home/susana/Documents/Ph_D/Methods_paper/")
#source("table_mod.R")
#master <- read.csv("master_methods_20200207.csv", sep = ",")
#Paired <- na.omit(subset(master, select = c(LabCode, IgA, IgA_serum,
#                                            IgG, IgG_serum, Lysozyme_serum,
#                                            Lysozyme_faeces, Age)))
#write.csv(Paired, "/home/susana/Documents/HyenaMan/ImmuneHyena/data/pairedBloodFaeces.csv", row.names=FALSE)
#write.csv(complete, "/home/susana/Documents/HyenaMan/ImmuneHyena/data/hyenadataset.csv", row.names=FALSE)
#write.csv(zoos, "/home/susana/Documents/HyenaMan/ImmuneHyena/data/zoodataset.csv", row.names=FALSE)
library("vegan")
library(MASS)
library(ggplot2)
library(ggpubr)
setwd("/home/susana/Documents/HyenaMan/ImmuneHyena")
complete=read.csv("data/hyenadataset.csv")
zoos=read.csv("data/zoodataset.csv")

#Hypotheses: 
#1) adults have distinct parasite and immune profiles than juveniles
#2) sex differences in immune and parasites in adults

#### preparing dataset
Adult <- complete[complete$Age == 3,]
Juv <- complete[complete$Age == 1,]

for (i in 1:nrow(Juv)) {
    Juv$survAd[i] <- if (Juv$age_death[i] < 730) 1 else 0
}
complete$rich_a <- as.numeric(complete$rich) - as.numeric(complete$ancy)
complete$Age <- as.factor(complete$Age)
complete$Sex <- as.factor(complete$Sex)
complete$Tri=as.numeric(complete$Tri)
complete$Tani=as.numeric(complete$Tani)
complete$Dipy=as.numeric(complete$Dipy)

# add year
complete$year <- substr(complete$S_date, 0, 4)
complete$year <- as.factor(complete$year)
# add clan
complete$Clan <- as.factor(substr(complete$ID, 0, 1))

### nMDS and permutation
# 18 individuals do not have any parasites. We need to remove those for the
# next steps.
nrow(na.omit(complete[,c(6:13)])[rowSums(na.omit(complete[,c(6:13)]))==0,])
# Non parametric multidimentional scaling
Pmds=metaMDS(na.omit(complete[,c(6:13)])[rowSums(na.omit(complete[,c(6:13)]))>0,])
Imds=metaMDS(na.omit(complete[,c(17:19)]))


PAxis=Pmds$points[,1:2]
PAxis=as.data.frame(PAxis)
PAxis$row=rownames(PAxis)
IAxis=Imds$points[,1:2]
IAxis=as.data.frame(IAxis)
IAxis$row=rownames(IAxis)
names(IAxis)=c("IMDS1", "IMDS2", "row")
names(PAxis)=c("PMDS1", "PMDS2", "row")
complete$row <- rownames(complete)
complete <- merge(complete, IAxis, by = "row", all = T)
complete <- merge(complete, PAxis, by = "row", all = T)
complete <- complete[, -1]
# Alice is not intersex, she is a female
complete$Sex[complete$ID=="P 657"]=2
complete$Sex=droplevels(complete$Sex)

## permutation tests for quantitative and qualitative variables. A lot of immune
# measures are NA's. For this reason, I test seperatedly Sex and Age and immune
# measures, to optimize sample size.
# Immune parameters first (quantitative)
dataP <- complete[complete.cases(complete$PMDS1),]
dataP$Sex=droplevels(dataP$Sex)

# permanova
k=complete.cases(dataP[,17:19])
P.dist = vegdist(dataP[k,6:13])
# permanovas are sensitive to the order of the variables, so let's shuffle a bit
adonis2(P.dist~year+Clan+Age+Sex+IgA+IgG+mucin, data=dataP[k,], permutations = 10000, method = "bray")
adonis2(P.dist~Sex+IgA+IgG+mucin+Age+year+Clan, data=dataP[k,], permutations = 10000, method = "bray")
adonis2(P.dist~IgA+IgG+mucin+year+Clan+Age+Sex, data=dataP[k,], permutations = 10000, method = "bray")
adonis2(P.dist~IgG+mucin+year+Clan+Age+Sex+IgA, data=dataP[k,], permutations = 10000, method = "bray")
adonis2(P.dist~mucin+year+Clan+Age+Sex+IgA+IgG, data=dataP[k,], permutations = 10000, method = "bray")
# we report this one
adonis2(P.dist~Age+Sex+IgA+IgG+mucin + year+Clan, data=dataP[k,], permutations = 10000, method = "bray")
table(dataP[k,]$Age)

# Models for individual Immune measures 
complete <- within(complete, year<-relevel(year, ref="2010"))
Amodel_lm <- with(complete, glm((IgA) ~ year + Sex + Ancy + Age + rich_a + Ancy:Age + Sex:Age))
summary(Amodel_lm)

# using gaussian distribution gives substancial overdispersion.
Amodel <- with(complete, glm.nb((IgA) ~ year +Clan + Sex + Ancy + Age + rich_a + Ancy:Age + Sex:Age))
AmodelSexAnc<-with(complete, glm.nb((IgA) ~ year +Clan + Sex + Ancy + Age + rich_a + Ancy:Age + Sex:Age + Sex:Ancy))
# negative binomial looks pretty good
# not we do LRT for each predictor
Amodel_Sex <- with(complete, glm.nb((IgA) ~ year + Clan +Ancy + Age + rich_a + Ancy:Age))
Amodel_SA <- with(complete, glm.nb((IgA) ~ year + Clan +Sex + Ancy + Age + rich_a + Ancy:Age))
Amodel_Ancy <- with(complete, glm.nb((IgA) ~ year + Clan +Sex + Age + rich_a + Sex:Age))
Amodel_Age <- with(complete, glm.nb((IgA) ~ year + Clan +Sex + Ancy + rich_a))
Amodel_rich <- with(complete, glm.nb((IgA) ~ year + Clan +Sex + Ancy + Age + Ancy:Age + Sex:Age))
AmodelAA <- with(complete, glm.nb((IgA) ~ year + Clan +Sex + Ancy + Age + rich_a+ Sex:Age))
Amodelyear <- with(complete, glm.nb((IgA) ~ Sex + Clan +Ancy + Age + rich_a + Ancy:Age + Sex:Age))
Amodelclan <- with(complete, glm.nb((IgA) ~ year + Sex + Ancy + Age + rich_a + Ancy:Age + Sex:Age))
Amodel0 <- with(complete, glm.nb((IgA) ~ 1))
sink("tmp/IgA_summ.txt")
summary(Amodel)
sink()
anova(Amodel, Amodel0)
anova(Amodel, Amodelyear)
anova(Amodel, Amodelclan)
anova(Amodel, Amodel_Sex)
anova(Amodel, Amodel_Ancy)
anova(Amodel, Amodel_Age)
anova(Amodel, Amodel_rich)
anova(Amodel, AmodelAA)
anova(Amodel, Amodel_SA)

# Gaussian distribution also looks pretty bad
Gmodel_lm <- with(complete, glm((IgG) ~ year + Clan +Sex + Ancy + Age + rich_a + Ancy:Age + Sex:Age))
Gmodel <- with(complete, glm.nb((IgG) ~ Sex + Ancy + Age + rich_a + year + Clan + Ancy:Age + Sex:Age))
GmodelSexAncy <- with(complete, glm.nb((IgG) ~ Sex + Ancy + Age + rich_a + year + Clan + Ancy:Age + Sex:Age + Sex:Ancy))
Gmodel_Sex <- with(complete, glm.nb((IgG) ~ year + Clan +Ancy + Age + rich_a + Ancy:Age))
Gmodel_SA <- with(complete, glm.nb((IgG) ~ year + Clan +Sex + Ancy + Age + rich_a + Ancy:Age))
Gmodel_Ancy <- with(complete, glm.nb((IgG) ~ year + Clan +Sex + Age + rich_a + Sex:Age))
Gmodel_Age <- with(complete, glm.nb((IgG) ~ year + Clan +Sex + Ancy + rich_a))
Gmodel_rich <- with(complete, glm.nb((IgG) ~ year + Clan +Sex + Ancy + Age + Ancy:Age + Sex:Age))
GmodelAA <- with(complete, glm.nb((IgG) ~ year + Clan +Sex + Ancy + Age + rich_a+ Sex:Age))
Gmodelyear <- with(complete, glm.nb((IgG) ~ Clan +Sex + Ancy + Age + rich_a + Ancy:Age + Sex:Age))
Gmodelclan <- with(complete, glm.nb((IgG) ~ year +Sex + Ancy + Age + rich_a + Ancy:Age + Sex:Age))
Gmodel0 <- with(complete, glm.nb((IgG) ~ 1))
sink("tmp/IgG_summ.txt")
summary(Gmodel)
sink()
anova(Gmodel, Gmodel0)
anova(Gmodel, Gmodel_Sex)
anova(Gmodel, Gmodel_Ancy)
anova(Gmodel, Gmodel_Age)
anova(Gmodel, Gmodel_rich)
anova(Gmodel, GmodelAA)
anova(Gmodel, Gmodel_SA)
anova(Gmodel, Gmodelclan)
anova(Gmodel, Gmodelyear)

Mmodel <- with(complete, glm.nb((mucin) ~ Sex + Ancy + Age + rich_a + year + Clan + Ancy:Age + Sex:Age))
MmodelSexAncy <- with(complete, glm.nb((mucin) ~ Sex + Ancy + Age + rich_a + year + Clan + Ancy:Age + Sex:Age + Sex:Ancy))
Mmodel_lm <- with(complete, glm((mucin) ~ year + Clan +Sex + Ancy + Age + rich_a + Ancy:Age + Sex:Age))
plot(Mmodel_lm)

Mmodel_Sex <- with(complete, glm.nb((mucin) ~ year + Clan +Ancy + Age + rich_a + Ancy:Age))
Mmodel_SA <- with(complete, glm.nb((mucin) ~ year + Clan +Sex + Ancy + Age + rich_a + Ancy:Age))
Mmodel_Ancy <- with(complete, glm.nb((mucin) ~ year + Clan +Sex + Age + rich_a + Sex:Age))
Mmodel_Age <- with(complete, glm.nb((mucin) ~ year + Clan +Sex + Ancy + rich_a))
Mmodel_rich <- with(complete, glm.nb((mucin) ~ year + Clan +Sex + Ancy + Age + Ancy:Age + Sex:Age))
MmodelAA <- with(complete, glm.nb((mucin) ~ year + Clan +Sex + Ancy + Age + rich_a+ Sex:Age))
Mmodelyear <- with(complete, glm.nb((mucin) ~ Clan +Sex + Ancy + Age + rich_a + Ancy:Age + Sex:Age))
Mmodelclan <- with(complete, glm.nb((mucin) ~ year + Sex + Ancy + Age + rich_a + Ancy:Age + Sex:Age))
Mmodel0 <- with(complete, glm.nb((mucin) ~ 1))
sink("tmp/Muc_summ.txt")
summary(Mmodel)
sink()
anova(Mmodel, Mmodel0)
anova(Mmodel, Mmodel_Sex)
anova(Mmodel, Mmodel_Ancy)
anova(Mmodel, Mmodel_Age)
anova(Mmodel, Mmodel_rich)
anova(Mmodel, MmodelAA)
anova(Mmodel, Mmodel_SA)
anova(Mmodel, Mmodelclan)
anova(Mmodel, Mmodelyear)
# plot(Amodel)

#### longevity analysis for juveniles
library(survival)
Juv <- Juv[Juv$age < 365,]
surv_object <- Surv(time=Juv$age_death/365, event = Juv$censor)

Juv$ancy <- as.factor(Juv$ancy)
Juv$CatAncy3[Juv$Ancy >= mean(na.omit(Juv$Ancy))] <- "high Ancylostoma"
Juv$CatAncy3[Juv$Ancy < mean(na.omit(Juv$Ancy))] <- "medium Ancylostoma"
Juv$CatAncy3[Juv$Ancy < quantile(na.omit(Juv$Ancy, 0.25))[2]] <- "low Ancylostoma"

Juv$CatAncy[Juv$Ancy >= mean(na.omit(Juv$Ancy))] <- "high Ancylostoma"
Juv$CatAncy[Juv$Ancy < mean(na.omit(Juv$Ancy))] <- "low Ancylostoma"

Juv1$CatCysto[Juv$Cysto >= mean(na.omit(Juv$Cysto))] <- "high Cystoisospora"
Juv1$CatCysto[Juv$Cysto < mean(na.omit(Juv$Cysto))] <- "low Cystoisospora"
summary(Juv1$Ancy)
hist(Juv1$Ancy)

Juv$CatIgA[Juv$IgA >= median(na.omit(Juv$IgA))] <- "high IgA"
Juv$CatIgA[Juv$IgA < median(na.omit(Juv$IgA))] <- "low IgA"

Juv$CatIgA3[Juv$IgA >= median(na.omit(Juv$IgA))] <- "high IgA"
Juv$CatIgA3[Juv$IgA < median(na.omit(Juv$IgA))] <- "medium IgA"
Juv$CatIgA3[Juv$IgA < quantile(na.omit(Juv$IgA, 0.25))[2]] <- "low IgA"
summary(as.factor(Juv$CatIgA3))

IgA_Juv <- survfit(surv_object ~ CatIgA, data = Juv)

Juv$CatIgG[Juv$IgG >= median(na.omit(Juv$IgG))] <- "high IgG"
Juv$CatIgG[Juv$IgG < median(na.omit(Juv$IgG))] <- "low IgG"
Juv$CatIgG3[Juv$IgG >= median(na.omit(Juv$IgG))] <- "high IgG"
Juv$CatIgG3[Juv$IgG < median(na.omit(Juv$IgG))] <- "medium IgG"
Juv$CatIgG3[Juv$IgG < quantile(na.omit(Juv$IgG, 0.25))[2]] <- "low IgG"

IgG_Juv <- survfit(surv_object ~ CatIgG, data = Juv)

Juv$Catmucin[Juv$mucin >= median(na.omit(Juv$mucin))] <- "high mucin"
Juv$Catmucin[Juv$mucin < median(na.omit(Juv$mucin))] <- "low mucin"
Juv$Catmucin3[Juv$mucin >= median(na.omit(Juv$mucin))] <- "high mucin"
Juv$Catmucin3[Juv$mucin < median(na.omit(Juv$mucin))] <- "medium mucin"
Juv$Catmucin3[Juv$mucin < quantile(na.omit(Juv$mucin, 0.25))[2]] <- "low mucin"

mucin_Juv <- survfit(surv_object ~ Catmucin, data = Juv)
Juv$Sex <- as.factor(Juv$Sex)
Juv1<- Juv[complete.cases(Juv$age),]
Juv1<- Juv1[complete.cases(Juv1$Ancy),]
Juv1<- Juv1[complete.cases(Juv1$CatIgA),]

Juv$ancy <- as.factor(Juv$ancy)
fit.coxphA <- coxph(surv_object ~ Ancy + age + CatIgA, data = Juv)
fit.coxphA2 <- coxph(surv_object ~ Ancy + age + IgA, data = Juv)
fit.coxphA3 <- coxph(surv_object ~ Ancy + age + CatIgA3, data = Juv)
fit.coxphA4 <- coxph(surv_object ~ Ancy * IgA + age, data = Juv)
fit.coxphA5 <- coxph(surv_object ~ Ancy * CatIgA + age, data = Juv)

ggforest(fit.coxphA)
ggforest(fit.coxphA3)
ggforest(fit.coxphA4)
ggforest(fit.coxphA5)

fit.coxphA
fit.coxphG <- coxph(surv_object ~ Ancy+ age + CatIgG, data = Juv)
fit.coxphG3 <- coxph(surv_object ~ Ancy+ age + CatIgG3, data = Juv)
fit.coxphG4 <- coxph(surv_object ~ Ancy * CatIgG3 + age, data = Juv)
fit.coxphGq <- coxph(surv_object ~ Ancy+ age + IgG, data = Juv)
fit.coxphGI <- coxph(surv_object ~ Ancy * IgG +age, data = Juv)
# survG<- ggforest(fit.coxphG, data = Juv)
ggforest(fit.coxphG)
ggforest(fit.coxphG3)
ggforest(fit.coxphG4)
ggforest(fit.coxphGq)
fit.coxphM <- coxph(surv_object ~ Ancy + age + Catmucin, data = Juv)
fit.coxphM3 <- coxph(surv_object ~ Ancy + age + Catmucin3, data = Juv)
fit.coxphMq <- coxph(surv_object ~ Ancy + age + mucin, data = Juv)
fit.coxphMqI <- coxph(surv_object ~ Ancy * mucin+ age, data = Juv)
fit.coxphM4 <- coxph(surv_object ~ Ancy * Catmucin + age, data = Juv)

ggforest(fit.coxphM, data = Juv)
ggforest(fit.coxphM3, data = Juv)
ggforest(fit.coxphM4, data = Juv)
ggforest(fit.coxphMq, data = Juv)

azfit <- cox.zph(fit.coxphA)
# plot(azfit)
gzfit <- cox.zph(fit.coxphG)
# plot(gzfit)
mzfit <- cox.zph(fit.coxphM)
# plot(mzfit)

ADs <- subset(zoos, grepl("AD", zoos$LabCode))
ADs$days <- c(0, 6, 9, 16, 18)
ADs <- droplevels(ADs)
cor.test(ADs$IgA, ADs$days, method = "spearman")
cor.test(ADs$IgG, ADs$days, method = "spearman")
cor.test(ADs$muc, ADs$days, method = "spearman")

## spearman for paired blood and faeces immune measures
Paired<- read.csv("data/pairedBloodFaeces.csv")
cor.test(Paired$IgA, Paired$IgA_serum, test = "spearman")
cor.test(Paired$IgG, Paired$IgG_serum, test = "spearman")
## save dataset

names(complete)
keep=c("LabCode", "Sex", "Age", "Ancy", "Spiro", "Cysto", "Tri", "Tani", "Dipy", "Meso", "Spiru", "IgA", "IgG", "mucin", "censor", "age", "Clan", "age_death", "rich_a", "year", "IMDS1", "IMDS2", "PMDS1", "PMDS2")
dryadComplete <- complete[, colnames(complete)%in% keep]

dryadComplete$Diphyllobothrium = dryadComplete$Spiro
dryadComplete$Spiro = NULL

names(zoos)
keep2=c("LabCode", "Sex.x", "Age", "IgA", "IgG", "mucin")
dryadZoo <- zoos[, colnames(zoos)%in% keep2]
dryadZoo$Sex = dryadZoo$Sex.x
dryadZoo$Sex.x = NULL

dryadPaired=Paired
dryadPaired$LabCode=c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9")

write.table(dryadPaired, "/home/susana/Documents/HyenaMan/ImmuneHyena/data/pairedBloodFaeces_dryad.txt", row.names=FALSE, sep="\t")
write.table(dryadComplete, "/home/susana/Documents/HyenaMan/ImmuneHyena/data/hyenadataset_dryad.txt", row.names=FALSE, sep="\t")
write.table(dryadZoo, "/home/susana/Documents/HyenaMan/ImmuneHyena/data/zoodataset_dryad.txt", row.names=FALSE, sep="\t")
    