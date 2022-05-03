# Plotting_methods
setwd("/home/susana/Documents/HyenaMan/ImmuneHyena")
source("bin/Ord_Models.R")

library(ggplot2)
library(showtext)
library(pec)
library(survminer)
library(ggeffects)
showtext_auto()
theme_set(theme_light(base_size = 10, base_family = "Sans"))

#figure 1
Fig1b <- ggplot(ADs, aes(days, (IgG))) +
  geom_point(shape = 21, size=3, colour = "black", fill="gray", stroke = 1) +
  scale_y_continuous(breaks = seq(0, 120, 10), name = "Total faecal IgG") +
  scale_x_continuous(breaks = seq(0, 18, 2), name = "Days post anaesthesia") +
  coord_cartesian(ylim=c(0, 120)) +
  theme_classic()
Fig1a <- ggplot(ADs, aes(days, (IgA))) +
  geom_point(shape = 21, size=3, colour = "black", fill="gray", stroke = 1) +
  scale_y_continuous(breaks = seq(0, 100000, 10000), name = "Total faecal IgA") +
  scale_x_continuous(breaks = seq(0, 18, 2), name = "Days post anaesthesia") +
  coord_cartesian(ylim=c(0, 100000)) +
  theme_classic()
Fig1c <- ggplot(ADs, aes(days, (mucin))) +
  geom_point(shape = 21, size=3, colour = "black", fill="gray", stroke = 1) +
  scale_y_continuous(breaks = seq(0, 100, 10), name = "Mucin (μmol OE)") +
  scale_x_continuous(breaks = seq(0, 18, 2), name = "Days anaesthesia") +
  coord_cartesian(ylim=c(0, 100)) +
  theme_classic()

pdf("figures/figure1.pdf",
     width = 7, height = 3.5, pointsize = 12)
ggarrange(Fig1a, Fig1b, Fig1c,ncol=3, labels = c("a", "b", "c"))
dev.off()

# plotting figure Validation nMDS
pdf("figures/figureS1.pdf", 
    useDingbats=FALSE,
    width=10, 
    height=10)
par(mfrow=c(2,2))
stressplot(Pmds, main="Parasite nMDS Shepard plot")
gofp=goodness(Pmds)
rownames(Pmds$species) = c("Ancylostoma", "Diphyllobothrium", "Cystoisospora",
                           "Trichuris", "Taniidae", "Dipylidium", "Mesocestoides", "Spirurida")
plot(Pmds, display = c("sites"), main = "Parasite nMDS goodness of fit")
points(Pmds, display = "sites", cex=gofp*50, pch = 21, bg="gray")
text(Pmds, display=c("species"), col="blue", cex=0.5)
stressplot(Imds, main="Immune nMDS Shepard plot")
gofi=goodness(Imds)
rownames(Imds$species) = c("IgA", "IgG", "Mucin")
plot(Imds, display = c("sites"), main = "Immune nMDS goodness of fit")
points(Imds, display = "sites", cex=gofi*200, pch = 21, bg="gray")
text(Imds, display=c("species"), col="blue", cex=0.5)
dev.off()


# plotting survival
IgAsurv <- ggadjustedcurves(fit.coxphA, data = Juv, method = "conditional", variable = "CatIgA",
                            palette = "jco",
                            legend.title = "Groups:",
                            tables.y.text = FALSE,
                            legend.labs = c("high IgA", "low IgA"),
                            legend=c(0.8, 0.8),
                            xlab= "Longevity (years)",
                            xlim=c(0,8)) +
  theme_classic()+
  theme(legend.position=c(0.8, 0.8))
jpeg("IgA_surv.jpeg",
     width = 3.5, height = 3.5, units = "in", pointsize = 11,
     res = 200)
IgAsurv
dev.off() 

pdf("Figure4.pdf",
     width = 4, height = 4, pointsize = 11)
IgAsurv
dev.off() 


IgGsurv <- ggadjustedcurves(fit.coxphG, data = Juv, method = "conditional", variable = "CatIgG", palette = "jco")
Mucsurv <- ggadjustedcurves(fit.coxphM, data = Juv, method = "conditional", variable = "Catmucin", palette = "jco")

# plot error from survivorship
### annoying package so we need to remove NA's from the table
Juv$age_death <- as.numeric(Juv$age_death)

Juvv <- Juv[complete.cases(Juv$CatAncy),]
Juvv <- Juvv[complete.cases(Juvv$IgA),]
surv_object <- Surv(time=Juvv$age_death, event = Juvv$censor)


# fit some candidate Cox models and compute the Kaplan-Meier estimate
models <- list("CoX.AA"=fit.coxphA <- coxph(Surv(age_death/365, censor) ~ Ancy + age + CatIgA, data = Juvv,x=TRUE,y=TRUE),
               "CoX.IgA"=fit.coxphA <- coxph(Surv(age_death/365, censor) ~ Ancy + age, data = Juvv,x=TRUE,y=TRUE))

# compute the apparent prediction error
PredErrorA <- pec(object=models,
                  formula= Surv(age_death/365, censor)~Ancy+age+CatIgA,
                  data=Juvv,
                  exact=TRUE,
                  cens.model="marginal",
                  splitMethod="Bootcv",
                  B=100,
                  verbose=TRUE)

print(PredErrorA)
summary(PredErrorA)
PredA <- plot(PredErrorA, what="BootCvErr", smooth = T, legend=F)

### annoying package so we need to remove NA's

Juv$age_death <- as.numeric(Juv$age_death)

Juvv <- Juv[complete.cases(Juv$Ancy),]
Juvv <- Juvv[complete.cases(Juvv$IgG),]

surv_object <- Surv(time=Juvv$age_death, event = Juvv$censor)
Juvv$CatAncy
# fit some candidate Cox models and compute the Kaplan-Meier estimate
models <- list("CoX.GA"=fit.coxphA <- coxph(Surv(age_death/365, censor) ~ Ancy + age + CatIgG, data = Juvv,x=TRUE,y=TRUE),
  "Cox.IgG"=fit.coxphA <- coxph(Surv(age_death/365, censor) ~ Ancy + age, data = Juvv,x=TRUE,y=TRUE))
# compute the apparent prediction error
GPredError <- pec(object=models,
                  formula= Surv(age_death/365, censor)~Ancy+age +CatIgG,
                  data=Juvv,
                  exact=TRUE,
                  cens.model="marginal",
                  splitMethod="Bootcv",
                  B=100,
                  verbose=TRUE)

print(GPredError)
summary(GPredError)
plot(GPredError, what="BootCvErr", smooth = T, legend=F)

### annoying package so we need to remove NA's
Juv$age_death <- as.numeric(Juv$age_death)

Juvv <- Juv[complete.cases(Juv$Ancy),]
Juvv <- Juvv[complete.cases(Juvv$Catmucin),]
Juvv$CatAncy
surv_object <- Surv(time=Juvv$age_death, event = Juvv$censor)
Juvv$Catmucin
# fit some candidate Cox models and compute the Kaplan-Meier estimate
models <- list(
               "CoX.Amuc"=fit.coxphA <- coxph(Surv(age_death/365, censor) ~ Ancy + age + Catmucin, data = Juvv,x=TRUE,y=TRUE),
               "Cox.muc"=fit.coxphA <- coxph(Surv(age_death/365, censor) ~ Ancy + age, data = Juvv,x=TRUE,y=TRUE))
# compute the apparent prediction error
mPredError <- pec(object=models,
                  formula= Surv(age_death/365, censor)~Ancy+age +Catmucin,
                  data=Juvv,
                  exact=TRUE,
                  cens.model="marginal",
                  splitMethod="Bootcv",
                  B=100,
                  verbose=TRUE)

print(mPredError)
summary(PredError)
mucerr<-plot(mPredError, what="BootCvErr", smooth = T, legend=F)
pdf("A_surverr.pdf",
    width = 3.5, height = 3.5, pointsize = 12)
plot(PredErrorA, what="BootCvErr", smooth = T, legend=T)
dev.off()

pdf("G_surverr.pdf",
    width = 3.5, height = 3.5)
plot(GPredError, what="BootCvErr", smooth = T, legend=F)
dev.off()

pdf("M_surverr.pdf",
    width = 3.5, height = 3.5)
plot(mPredError, what="BootCvErr", smooth = T, legend=F)
dev.off()

jpeg("A_surverr.jpeg",
     width = 3.5, height = 3.5, pointsize = 12,
     units = "in", res = 200)
plot(PredErrorA, what="BootCvErr", smooth = T, legend=F)
dev.off()

jpeg("G_surverr.jpeg",
     width = 3.5, height = 3.5, pointsize = 12,
     units = "in", res = 200)
plot(GPredError, what="BootCvErr", smooth = T, legend=F)
dev.off()

jpeg("m_surverr.jpeg",
     width = 3.5, height = 3.5, pointsize = 12,
     units = "in", res = 200)
plot(mPredError, what="BootCvErr", smooth = T, legend=F)
dev.off()


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
### predict IgA for Ancy-Age
Adata <- complete[complete.cases(complete$IgA),]
Adata <- Adata[complete.cases(Adata$Ancy),]
newdata1 <- with(Adata, data.frame
                 (Sex = factor('1'), Ancy = seq(min(Ancy), max(Ancy), 110),
                   Age = factor("1"), rich_a = median(rich_a),
                   year = factor("2010"),
                   Clan = factor("I")))
newdata2 <- with(Adata, data.frame
                 (Sex = factor('1'), Ancy = seq(min(Ancy), 1575, 110),
                   Age = factor("3"), rich_a = median(rich_a),
                   year = factor("2010"),
                   Clan = factor("I")))
newdata1 <- cbind(newdata1, predict(Amodel, newdata1, type = "link", se.fit=TRUE))
newdata2 <- cbind(newdata2, predict(Amodel, newdata2, type = "link", se.fit=TRUE))
newdata1 <- within(newdata1, {
  IgA <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
newdata2 <- within(newdata2, {
  IgA <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
new1 <- rbind(newdata1, newdata2)

Age_Ancy_A <-ggplot(new1, aes(x = Ancy, y = IgA, colour = Age)) +
  geom_line(aes(color=Age)) +
  geom_point(data=Adata, aes(x=Ancy, y=IgA, colour=Age), alpha=0.4)+
  scale_x_continuous("Ancylostoma egg load") +
  scale_y_continuous("Predicted IgA", breaks = seq(0, 100000, 10000)) +
  coord_cartesian(ylim=c(0,100000), xlim = c(0, 10000)) +
  geom_ribbon(aes(ymin = LL, ymax = UL),
              linetype=2,
              show.legend = TRUE,
              alpha=0.1) +
  scale_color_manual(name = "Age class", labels =c("juvenile", "adult"),
                     values=c('#999999','#E69F00'))+
  theme_classic()

# Predicted values of IgG counts in relation to Age, for males and females, holding Ancylostoma and richness at mean.
# library(sjPlot)
Gdata <- complete[complete.cases(complete$IgG),]
Gdata <- Gdata[complete.cases(Gdata$Ancy),]
summary(Gmodel)
#Ggmodel <- with(Gdata, glm.nb((IgG) ~ Sex + Ancy + Age + rich_a + Ancy:Age + Sex:Age))
Gdata$Sex
newdata1 <- with(Gdata, data.frame
                 (Sex = rep(c('1', '2'), 57), Ancy = median(Ancy),
                                      Age = rep(c("1"), 57), rich_a = median(rich_a),
                   year = rep(c("2010"), 57),
                   Clan = rep(c("I"), 57)))
newdata2 <- with(Gdata, data.frame
                 (Sex = rep(c('1', '2'), 57), Ancy = median(Ancy),
                   Age = rep(c("3"), 57), rich_a = median(rich_a),
                   year = rep(c("2010"), 57),
                   Clan = rep(c("I"), 57)))

# newdata1$juvenile <- predict.glm(Gmodel, newdata1, type = "response")
newdata1 <- cbind(newdata1, predict(Gmodel, newdata1, type = "link", se.fit=TRUE))
# newdata2$adult <- predict.glm(Gmodel, newdata2, type = "response")
newdata2 <- cbind(newdata2, predict(Gmodel, newdata2, type = "link", se.fit=TRUE))
newdata1 <- within(newdata1, {
  IgG <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
newdata2 <- within(newdata2, {
  IgG <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
newdata1 <- newdata1[-c(1:57),]
newdata2 <- newdata2[-c(2:58),]

new1 <- rbind(newdata1, newdata2)

Age_Sex_G <- ggplot(new1, aes(x = Age, y = IgG, colour = Sex)) +
  geom_point(position="dodge",size=3, shape=15) +
  geom_jitter(data=Gdata, aes(x=Age, y=IgG, colour=Sex), alpha=0.4, width = 0.3, height = 0) +
  scale_x_discrete("Age class", labels=c("juvenile", "adult")) +
  scale_y_continuous("Predicted IgG", breaks = seq(0, 600, 100)) +
  coord_cartesian(ylim=c(0,600)) +
  geom_errorbar(aes(ymin = LL, ymax = UL, width=.1), alpha=0.1,
                position=position_dodge(0.02)) +
  scale_color_manual(name = "Sex", labels =c("male", "female"),
  values=c('#999999','#E69F00'))+
  theme_classic()
# Predicted IgG for Ancy-Age
newdata1 <- with(Gdata, data.frame
                 (Sex = factor('1'), Ancy = seq(min(Ancy), max(Ancy), 110),
                  Age = factor("1"), rich_a = median(rich_a),
                  year = factor("2010"),
                  Clan = factor("I")))
newdata2 <- with(Gdata, data.frame
                 (Sex = factor('1'), Ancy = seq(min(Ancy), 1575, 110),
                  Age = factor("3"), rich_a = median(rich_a),
                  year = factor("2010"),
                  Clan = factor("I")))
newdata1 <- cbind(newdata1, predict(Gmodel, newdata1, type = "link", se.fit=TRUE))
newdata2 <- cbind(newdata2, predict(Gmodel, newdata2, type = "link", se.fit=TRUE))
newdata1 <- within(newdata1, {
  IgG <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
newdata2 <- within(newdata2, {
  IgG <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
new1 <- rbind(newdata1, newdata2)
Age_Ancy_G <- ggplot(new1, aes(x = Ancy, y = IgG, colour = Age)) +
  geom_line(aes(color=Age)) +
  geom_point(data=Gdata, aes(x=Ancy, y=IgG, colour=Age), alpha=0.4)+
  scale_x_continuous("Ancylostoma egg load") +
  scale_y_continuous("Predicted IgG") +
  coord_cartesian(ylim=c(0,3000), xlim = c(0, 10000)) +
  geom_ribbon(aes(ymin = LL, ymax = UL),
              linetype=2,
              show.legend = TRUE,
              alpha=0.1) +
  scale_color_manual(name = "Age class", labels =c("juvenile", "adult"),
                     values=c('#999999','#E69F00'))+
  theme_classic()
# Predicted values of mucin in relation to Age and Ancy,females, richness at mean.
# library(sjPlot)
mdata <- complete[complete.cases(complete$mucin),]
mdata <- mdata[complete.cases(mdata$Ancy),]
summary(Mmodel)
mdata$Sex
newdata1 <- with(mdata, data.frame
                 (Sex = factor('1'), Ancy = seq(min(Ancy), max(Ancy), 151),
                   Age = factor("1"), rich_a = median(rich_a),
                   year = factor("2010"),
                   Clan = factor("I")))
newdata2 <- with(mdata, data.frame
                 (Sex = factor('1'), Ancy = seq(min(Ancy), 1575, 151),
                   Age = factor("3"), rich_a = median(rich_a),
                   year = factor("2010"),
                   Clan = factor("I")))
newdata1 <- cbind(newdata1, predict(Mmodel, newdata1, type = "link", se.fit=TRUE))
newdata2 <- cbind(newdata2, predict(Mmodel, newdata2, type = "link", se.fit=TRUE))
newdata1 <- within(newdata1, {
  muc <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
newdata2 <- within(newdata2, {
  muc <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
new1 <- rbind(newdata1, newdata2)
Age_Ancy_M <-ggplot(new1, aes(x = Ancy, y = (muc), colour = Age)) +
  geom_line(aes(color=Age)) +
  geom_point(data=mdata, aes(x=Ancy, y=mucin, colour=Age), alpha=0.4)+
  scale_x_continuous("Ancylostoma egg load") +
  scale_y_continuous("Predicted mucin", breaks = seq(0, 500, 100)) +
  coord_cartesian(ylim=c(0,500), xlim = c(0, 10000)) +
  geom_ribbon(aes(ymin = LL, ymax = UL),
              linetype=2,
              show.legend = TRUE,
              alpha=0.1) +
  scale_color_manual(name = "Age class", labels =c("juvenile", "adult"),
                     values=c('#999999','#E69F00'))+
  theme_classic()
# Predicted values of muci in relation to Age, for males and females, holding Ancylostoma and richness at mean.
newdata1 <- with(mdata, data.frame
                 (Sex = rep(c('1', '2'), 57), Ancy = median(Ancy),
                   Age = rep(c("1"), 57), rich_a = median(rich_a),
                   year = rep(c("2010"), 57),
                   Clan = rep(c("I"), 57)))
newdata2 <- with(mdata, data.frame
                 (Sex = rep(c('1', '2'), 57), Ancy = median(Ancy),
                   Age = rep(c("3"), 57), rich_a = median(rich_a),
                   year = rep(c("2010"), 57),
                   Clan = rep(c("I"), 57)))

newdata1 <- cbind(newdata1, predict(Mmodel, newdata1, type = "link", se.fit=TRUE))
newdata2 <- cbind(newdata2, predict(Mmodel, newdata2, type = "link", se.fit=TRUE))
newdata1 <- within(newdata1, {
  muc <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
newdata2 <- within(newdata2, {
  muc <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
newdata1 <- newdata1[-c(1:57),]
newdata2 <- newdata2[-c(2:58),]

new1 <- rbind(newdata1, newdata2)
names(mdata)
Age_Sex_M <- ggplot(new1, aes(x = Age, y = muc, colour = Sex)) +
  geom_point(position="dodge",size=3, shape=15) +
  geom_jitter(data=mdata, aes(x=Age, y=mucin, colour=Sex), alpha=0.4, width = 0.3, height = 0) +
  scale_x_discrete("Age class", labels=c("juvenile", "adult")) +
  scale_y_continuous("Predicted mucin", breaks = seq(0, 600, 100)) +
  coord_cartesian(ylim=c(0,600)) +
  geom_errorbar(aes(ymin = LL, ymax = UL, width=.1), alpha=0.1,
                position=position_dodge(0.02)) +
  scale_color_manual(name = "Sex", labels =c("male", "female"),
                     values=c('#999999','#E69F00'))+
  theme_classic()
# ploting mucin ~ richness
newdata1 <- with(mdata, data.frame(
  Sex = factor('1'),
  Ancy = median(Ancy),
  Age = factor("1"), rich_a = seq(0, 5, 1),
  year = factor("2010"),
  Clan = factor("I")))
newdata2 <- with(mdata, data.frame(
  Sex = factor('1'),
  Ancy = median(Ancy),
  Age = factor("3"), rich_a = seq(0, 5, 1),
  year = factor("2010"),
  Clan = factor("I")))
newdata1 <- cbind(newdata1, predict(Mmodel, newdata1, type = "link", se.fit=TRUE))
newdata2 <- cbind(newdata2, predict(Mmodel, newdata2, type = "link", se.fit=TRUE))
newdata1 <- within(newdata1, {
  muc <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
newdata2 <- within(newdata2, {
  muc <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
new1 <- rbind(newdata1, newdata2)
Age_Rich_M <- ggplot(new1, aes(x = rich_a, y = muc, colour = Age)) +
  geom_line(aes(color=Age)) +
  geom_point(data=mdata, aes(x=rich_a, y=mucin, colour=Age), alpha=0.4) +
  scale_x_continuous("Parasite richness") +
  scale_y_continuous("Predicted mucin", breaks = seq(0, 700, 100)) +
  coord_cartesian(ylim=c(0,700)) +
  geom_ribbon(aes(ymin = LL, ymax = UL),
              linetype=2,
              show.legend = TRUE,
              alpha=0.1) +
  scale_color_manual(name = "Age class", labels =c("juvenile", "adult"),
                     values=c('#999999','#E69F00'))+
  theme_classic()

# plotting only a window from the plot before
newdata1 <- with(mdata, data.frame
                 (Sex = factor('1'), Ancy = seq(min(Ancy), 1575, 151),
                   Age = factor("1"), rich_a = median(rich_a)))
newdata2 <- with(mdata, data.frame
                 (Sex = factor('1'), Ancy = seq(min(Ancy), 1575, 151),
                   Age = factor("3"), rich_a = median(rich_a)))
newdata1 <- cbind(newdata1, predict(Mmodel, newdata1, type = "link", se.fit=TRUE))
newdata2 <- cbind(newdata2, predict(Mmodel, newdata2, type = "link", se.fit=TRUE))
newdata1 <- within(newdata1, {
  muc <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
newdata2 <- within(newdata2, {
  muc <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})
new1 <- rbind(newdata1, newdata2)
Age_Ancy_Mm <- ggplot(new1, aes(x = Ancy, y = (muc), colour = Age)) +
  geom_line(aes(color=Age)) +
  scale_x_continuous("Ancylostoma egg load") +
  scale_y_continuous("Predicted mucin", breaks = seq(0, 500, 100)) +
  coord_cartesian(ylim=c(0,500)) +
  geom_ribbon(aes(ymin = LL, ymax = UL),
              linetype=2,
              show.legend = TRUE,
              alpha=0.1) +
  scale_color_manual(name = "Age class", labels =c("juvenile", "adult"),
                     values=c('#999999','#E69F00'))+
  theme_classic()
Age_Ancy_Mm


pdf("figures/Revfigure2.pdf",
     width = 5, height = 7, pointsize = 11)
ggarrange(Age_Ancy_A, Age_Ancy_G, Age_Ancy_M, ncol = 1, labels=c("a", "b", "c"))
dev.off()

pdf("figures/Revfigure3.pdf",
    width = 5, height = 5, pointsize = 11)
ggarrange(Age_Rich_M,
          ggarrange(Age_Sex_G, Age_Sex_M, labels = c("b","c")),
          nrow=2,labels="a")
dev.off()

complete$Sex_Age[complete$Age == "1" & complete$Sex == "1"] <- "juvenile male"
complete$Sex_Age[complete$Age == "1" & complete$Sex == "2"] <- "juvenile female"
complete$Sex_Age[complete$Age == "3" & complete$Sex == "1"] <- "adult male"
complete$Sex_Age[complete$Age == "3" & complete$Sex == "2"] <- "adult female"

### nMDS
# parasite

ordiplot(Pmds, type="text")
ordiplot(Pmds, display = "species")
ordiplot(Pmds, display="sites")
PmdsDF<- complete[complete.cases(complete$PMDS1),]
ImdsDF<- complete[complete.cases(complete$IMDS1),]

rownames(PmdsDF)==rownames(Pmds$points)
plot(Imds)
##
spp <- as.data.frame(Pmds$species)
spp$species <- rownames(spp)

library(ggrepel)
PmdsDF$group=PmdsDF$Age
PmdsDF.mean=aggregate(PmdsDF[,c("PMDS1", "PMDS2")], list(group=PmdsDF$group), mean)
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(PmdsDF$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(PmdsDF[PmdsDF$group==g,],
                                                   veganCovEllipse(cov.wt(cbind(PMDS1,PMDS2),wt=rep(1/length(PMDS1),length(PMDS1)))$cov,center=c(mean(PMDS1),mean(PMDS2)))))
                                ,group=g))
}
df_ell

Pmds <-ggplot(complete, aes(x = PMDS1, y = PMDS2)) +
  geom_point(size = 1, alpha=0.8, aes(color = Age)) +
  #  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  labs(x = "nMDS1", y = "nMDS2", color="age class") +
  scale_colour_manual(values = c("#999999", "#E69F00"), labels=c("juvenile", "adult"))+
#  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  geom_segment(data = spp, aes(x = 0, xend=MDS1, y=0, yend=MDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey10", lwd=0.3)+
  ggrepel::geom_text_repel(data = spp, aes(x=MDS1, y=MDS2, label = species),
                           cex = 3, direction = "both",
                           segment.size = 0.25)+
  geom_path(data = df_ell, aes(x = PMDS1, y = PMDS2, group = group, colour=group),
            alpha=0.8)

ImdsDF$group=ImdsDF$Age
ImdsDF.mean=aggregate(ImdsDF[,c("IMDS1", "IMDS2")], list(group=ImdsDF$group), mean)

df_ellI <- data.frame()
for(g in levels(ImdsDF$group)){
    df_ellI <- rbind(df_ellI, cbind(as.data.frame(with(ImdsDF[ImdsDF$group==g,],
                  veganCovEllipse(cov.wt(cbind(IMDS1,IMDS2),
                                  wt=rep(1/length(IMDS1),length(IMDS1)))$cov,
                                  center=c(mean(IMDS1),mean(IMDS2))))),group=g))
  }
df_ellI
imm <- as.data.frame(Imds$species)
imm$imm <- rownames(imm)

Imds <-ggplot(complete, aes(x = IMDS1, y = IMDS2)) +
    geom_point(size = 1, alpha=0.8, aes(color = Age)) +
    #  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
    labs(x = "nMDS1", y = "nMDS2", color="age class") +
    scale_colour_manual(values = c("#999999", "#E69F00"), labels=c("juvenile", "adult"))+
#    coord_fixed()+
    theme_classic()+ 
    theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
    geom_segment(data = imm, aes(x = 0, xend=MDS1, y=0, yend=MDS2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "grey10", lwd=0.3)+
    ggrepel::geom_text_repel(data = imm, aes(x=MDS1, y=MDS2, label = imm),
                             cex = 3, direction = "both",
                             segment.size = 0.25)+
    geom_path(data = df_ellI, aes(x = IMDS1, y = IMDS2, group = group, colour=group),
              alpha=0.8)
head(df_ellI)
pdf("figures/Revfigure4.pdf",
    width = 9, height = 4, pointsize = 11)
ggarrange(Pmds, Imds, nrow=1,
          ncol=2, labels=c("a", "b"))
dev.off()

