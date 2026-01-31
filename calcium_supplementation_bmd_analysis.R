library(nlme)

BoneWide = read.csv("./BoneWide.csv")
BoneLong = read.csv("./BoneLong.csv")

BoneWide$Trt = as.factor(BoneWide$Trt)
BoneLong$Trt = as.factor(BoneLong$Trt)


BoneWide$Trt = relevel(BoneWide$Trt, ref="P")
BoneLong$Trt = relevel(BoneLong$Trt, ref="P")

## Complete cases:
ID.comp = BoneWide$ID[complete.cases(BoneWide)]
BoneLong.comp = BoneLong[BoneLong$ID %in% ID.comp,]
BoneWide.comp = BoneWide[BoneWide$ID %in% ID.comp,]


### Missing data
BoneWide[!complete.cases(BoneWide),]

na.count = function(x){ return( sum(is.na(x)) ) }
NumberOfObservations = 5-apply(BoneWide, 1, na.count)/2
table(NumberOfObservations)

### Summary statistics

# Summary for Trt, BMICat
table(BoneWide$Trt)
table(BoneWide$BMICat)

# Mean, SD, and number missing by visit for Age and BMD:

# Summary statistics for age:
tapply(BoneLong$Age, list(BoneLong$Trt,BoneLong$Visit), mean, na.rm=TRUE)
tapply(BoneLong$Age, list(BoneLong$Trt,BoneLong$Visit), sd, na.rm=TRUE)
# Summary statistics for BMD:
tapply(BoneLong$BMD, list(BoneLong$Trt,BoneLong$Visit), mean, na.rm=TRUE)
tapply(BoneLong$BMD, list(BoneLong$Trt,BoneLong$Visit), sd, na.rm=TRUE)
# Number missing:
tapply(BoneLong$BMD, list(BoneLong$Trt,BoneLong$Visit), na.count)

boxplot(Age ~ Visit, data=BoneLong)

### Plots

# Mean BMD by treatment and visit number:
means = tapply(BoneLong$BMD, list(BoneLong$Trt,BoneLong$Visit), mean, na.rm=TRUE)
plot(c(1:5),means[1,],type="o",xlab="Visit Number", ylab="Mean Bone Mineral Density", ylim=c(0.8,1), main="Mean BMD by Visit and Treatment")
lines(c(1:5),means[2,],type="o",xlab="Visit Number", ylab="Mean Bone Mineral Density", lty=2,pch=2)
legend("bottomright",c("Placebo","Calcium"),lty=c(1,2),pch=c(1,2))

# BMD by treatment and age:
plot(BMD ~ Age, data=BoneLong, pch=as.numeric(BoneLong$Trt), col=as.numeric(BoneLong$Trt)+1, main="BMD versus Age by Treatment",ylim=c(0.6,1.2))
lines(lowess(BoneLong.comp$Age[BoneLong.comp$Trt=="P"], BoneLong.comp$BMD[BoneLong.comp$Trt=="P"]),col=2)
lines(lowess(BoneLong.comp$Age[BoneLong.comp$Trt=="C"], BoneLong.comp$BMD[BoneLong.comp$Trt=="C"]),col=3, lty=2)
legend("bottomright", c("Placebo","Calcium"), lty=c(1,2), col=c(2,3), pch=c(1,2))

# Additional plots
plot(BMD ~ relevel(as.factor(BMICat),ref="Under"), data=BoneLong, main="BMD versus BMI Category",col="lightblue",xlab="BMI Category",ylab="Bone Mineral Density")
tapply(BoneLong$BMD, BoneLong$BMICat, median, na.rm=TRUE)

mosaicplot( ~ factor(Trt, levels=c("P","C"), labels=c("Placebo","Calcium")) + relevel(as.factor(BMICat),ref="Under"), data=BoneWide, xlab="Treatment",ylab="BMI Category Proportions", main="Mosaic Plot of Treatment by BMI Category", las=1, col=c("red","blue","green"))


tab = table(relevel(as.factor(BoneWide$BMICat),ref="Under"),factor(BoneWide$Trt, levels=c("P","C"), labels=c("Placebo","Calcium")))
barplot(tab,xlab="Treatment",ylab="Frequency",legend=TRUE,ylim=c(0,50),beside=TRUE, main="Bar Plot of BMI Category by Treatment")

plot(Age ~ relevel(as.factor(BMICat),ref="Under"), data=BoneLong, main="Age versus BMI Category",col="lightblue",xlab="BMI Category",ylab="Age")


### Models
library(nlme)

mod1 = gls(BMD ~ Trt+I(Visit-1)+Trt*I(Visit-1), correlation=corCompSymm(, form=~Visit | ID), method="REML", data=BoneLong.comp)
summary(mod1)

mod1.1 = gls(BMD ~ Trt+I(Visit-1)+Trt*I(Visit-1), correlation=corAR1(, form=~Visit | ID), method="REML", data=BoneLong.comp)
summary(mod1.1)


mod1.ML = gls(BMD ~ Trt+I(Visit-1)+Trt*I(Visit-1), correlation=corCompSymm(, form=~Visit | ID), method="ML", data=BoneLong.comp)

mod1.1.ML = gls(BMD ~ Trt+I(Visit-1)+Trt*I(Visit-1), correlation=corAR1(, form=~Visit | ID), method="ML", data=BoneLong.comp)


mod0.ML = gls(BMD ~ Trt+I(Visit-1), correlation=corCompSymm(, form=~Visit | ID), method="ML", data=BoneLong.comp)

mod0.1.ML = gls(BMD ~ Trt+I(Visit-1), correlation= corAR1(, form=~Visit | ID), method="ML", data=BoneLong.comp)

anova(mod0.ML, mod1.ML)

anova(mod0.1.ML, mod1.1.ML)

anova(mod1.1, mod1)

mod2 = gls(BMD ~ Trt+I(Visit-1)+Trt*I(Visit-1), correlation=corCompSymm(, form=~Visit | ID), method="REML", weights=varIdent(form = ~1|Visit), data=BoneLong.comp)

mod2.1 = gls(BMD ~ Trt+I(Visit-1)+Trt*I(Visit-1), correlation=corAR1(, form=~Visit | ID), method="REML", weights=varIdent(form = ~1|Visit), data=BoneLong.comp)

anova(mod1, mod2)

anova(mod1.1, mod2.1)

getVarCov(mod1)

getVarCov(mod1.1)



mod3 = lme(BMD ~ Trt+Age+Trt*Age, random = ~ 1|ID, data=BoneLong, na.action=na.omit, method="ML")
summary(mod3)

mod3.1 = lme(BMD ~ Trt , random = ~ 1|ID, data=BoneLong, na.action=na.omit, method="ML")
summary(mod3.1)

mod3s = lme(BMD ~ Trt+Age+Trt*Age, random = ~ 1+Age|ID, data=BoneLong, na.action=na.omit, method="ML")
summary(mod3s)

anova(mod3, mod3s)


mod3.REML = lme(BMD ~ Trt+Age+Trt*Age, random = ~ 1|ID, data=BoneLong, na.action=na.omit, method="REML")

mod3s.REML = lme(BMD ~ Trt+Age+Trt*Age, random = ~ 1+Age|ID, data=BoneLong, na.action=na.omit, method="REML")

anova(mod3.REML, mod3s.REML)


mod4.0 = lme(BMD ~ Trt+Age+Trt*Age + BMICat, random = ~ 1 + Age|ID, data=BoneLong, na.action=na.omit, method="ML")

mod4.1 = lme(BMD ~ Trt+Age+Trt*Age, random = ~ 1 + Age|ID, data=BoneLong, na.action=na.omit, method="ML")

anova(mod4.0, mod4.1)
anova(mod4.1, mod4.0)


mod4.0.1 = lme(BMD ~ Trt+Age+Trt*Age + BMICat, random = ~ 1 |ID, data=BoneLong, na.action=na.omit, method="ML")

mod4.1.1 = lme(BMD ~ Trt+Age+Trt*Age, random = ~ 1 |ID, data=BoneLong, na.action=na.omit, method="ML")

anova(mod4.0.1, mod4.1.1)


mod4.0.REML = lme(BMD ~ Trt+Age+Trt*Age + BMICat, random = ~ 1 + Age|ID, data=BoneLong, na.action=na.omit, method="REML")

mod4.1.REML = lme(BMD ~ Trt+Age+Trt*Age + BMICat, random = ~ 1  |ID, data=BoneLong, na.action=na.omit, method="REML")

anova(mod4.0.REML, mod4.1.REML)

mod4.0 = lme(BMD ~ Trt+Age+Trt*Age + BMICat, random = ~ 1 + Age|ID, data=BoneLong, na.action=na.omit, method="ML")
summary(mod4.0)



