

ldata = read.delim("Part_2/data/lipid_data.txt", header=T, as.is=T, check.names=F)
sdata = read.delim("Part_2/data/sample_data.txt", header=T, as.is=T, check.names=F) 

colnames(ldata) = gsub("\\(-H20\\)", "", colnames(ldata))

### First match the sample IDs between the two data sets
all(ldata$ID %in% sdata$ID)
all(sdata$ID %in% ldata$ID)
all(ldata$ID == sdata$ID) ## All sample IDs are aligned! 

################################
### Plotting sample meta data
### while learning R syntax
################################

#install.packages("scales")
library(scales)

attach(sdata)
nsamples = nrow(sdata)
  
hist(BMI, breaks=50)
hist(BMI[DM == 1], breaks=50, add=TRUE, col=2)

plot(SBP ~ Age, pch=19)
plot(SBP ~ Age, cex=.3, pch=19, col=alpha(2, 0.4))  ### red, opaque dots

par(mfrow=c(2,2))
plot(SBP ~ Age, cex=.3, pch=19, col=alpha(2, 0.4))
plot(HbA1c ~ Age, cex=.3, pch=19, col=alpha(2, 0.4))
plot(TG ~ Age, cex=.3, pch=19, col=alpha(2, 0.4))
plot(LDL ~ Age, cex=.3, pch=19, col=alpha(2, 0.4))
dev.off()

## Boxplots
boxplot(BMI ~ Gender)
boxplot(HbA1c ~ DM, boxwex=0.3)

################################
### Checking categorial data
################################

## Chi-squared test for categorical data
## and see if there is any correlation / association
tab = table(DM, Gender)
chisq.test(tab)  ### not significant

## Discretizing a continuously scaled variable 
bmi.new = rep(NA, nsamples)
bmi.new[BMI <= 18.5] = 1
bmi.new[BMI > 18.5 & BMI <= 23] = 2
bmi.new[BMI > 23 & BMI <= 27.5] = 3
bmi.new[BMI > 27.5] = 4
tab = table(DM, bmi.new)
chisq.test(tab)  
### significant

detach(sdata)



################################
### Projection of lipidomics data
################################
rownames(ldata) = ldata$ID
ldata = ldata[,-1]

any(ldata <= 0)  ### check before log transform
sum(ldata <= 0)  ### turns out only one value

######### In case there is zero, plug-in with a minimum value for each analyte
for(k in 1:ncol(ldata)) {
  xvec = ldata[,k]
  xpos = xvec[xvec > 0]
  ximpute = 0.9 * min(xpos)
  zid = xvec <= 0
  if(any(zid)) {
    cat(colnames(ldata)[k])
    xvec[zid] = ximpute 
    ldata[,k] = xvec
  }
}

tmp = log2(ldata)
#tmp = ldata   ### If you prefer working with raw data, do this instead. 

tmp.pca = prcomp(tmp)
vv = tmp.pca$sdev^2
vv = vv / sum(vv) * 100
vv = round(vv, 2)
print(vv)

ccc = rep("gray", nsamples)
ccc[sdata$DM == 1] = "red"

Thres = max(abs(tmp.pca$x[,1]))
XLAB = paste("PC1 (", vv[1], "%)", sep="")
YLAB = paste("PC2 (", vv[2], "%)", sep="")
library(scales)

pdf("Part_2/output/PCAplot.pdf", height=5.5, width=5, useDingbats = FALSE)
plot(tmp.pca$x[,1], tmp.pca$x[,2], col=alpha(ccc,0.2), pch=19,
     xlim=c(-Thres,Thres), ylim=c(-Thres,Thres), 
     xlab=XLAB, ylab=YLAB, 
     main="MEC cohort (N=2,299)", cex=0.8)
legend("bottomright", c("No event","Incident DM"), pch=19, col=c("gray","red"), cex=0.8)
abline(v=0, lty=2)
abline(h=0, lty=2)
dev.off()

############### Looks promising along PC2 axis
############### Draw the entire data (on a relative scale)
tmp.ctr = sweep(tmp, 2, apply(tmp, 2, median))
### We are normalizing each lipid by its own median value
### so that all lipids are on a comparable scale. 

#install.packages("gplots")
library(gplots)
pdf("Part_2/output/heatmap.pdf", height=20, width=25, useDingbats = FALSE)
heatmap.2(as.matrix(t(tmp.ctr)), trace="n", col=bluered(20), breaks=seq(-2,2,by=0.2), 
          distfun=function(x) as.dist(1-cor(t(x))), 
          hclustfun=function(x) hclust(x, method="average"), 
          ColSideColors = ccc, mar=c(10,10))
dev.off()

############### PLS-DA analysis (for fun)
############### Explain why PLS-DA is a "supervised" analysis
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("mixOmics")

library(mixOmics)
sample.class = ifelse(sdata$DM == 1, "Case", "Control")
X = as.matrix(tmp)
Y = sample.class
tmp.out = plsda(X=X, Y=Y, ncomp=2)
pdf("Part_2/output/PLSDA_projection.pdf")
plotIndiv(tmp.out, ind.names = TRUE, ellipse = TRUE, legend = TRUE)
dev.off()


##########################################
### Correlation with clinical parameters
##########################################

### We first re-arrange the columns of the sample data matrix
### so that all continous variables are shifted to the right
sdata = sdata[,c(1:3,5,4,6:14)]
cp = colnames(sdata)[5:14]

### Now we compute correlations between clinical data and lipids in a one liner
### Recall that ``tmp'' object holds lipidomic data in log2 scale
cormat = cor(sdata[,5:14], tmp, use="pairwise.complete.obs")

### Plot it in a heatmap
pdf("Part_2/output/cor_heatmap.pdf", height=20, width=6, useDingbats = FALSE)
heatmap.2(as.matrix((cormat)), trace="n", main="At baseline",
          col=bluered(20), breaks=seq(-1,1,by=0.1), 
          #distfun=function(x) as.dist(1-cor(t(x))), 
          hclustfun=function(x) hclust(x, method="average"), 
          mar=c(10,10))
dev.off()


##################################
### Hypothesis testing
### with multiple testing correction
##################################
nlipid = ncol(tmp)
nsample = nrow(tmp)
lipid.name = colnames(tmp)
tmp2 = tmp  
### We will standardize the data by mean centering 
### and setting standard deviation to 1 for each lipid

for(k in 1:nlipid) {
  mm = mean(tmp[,k], na.rm=TRUE)
  ss = sd(tmp[,k], na.rm=TRUE)
  tmp2[,k] = (tmp[,k] - mm) / ss
}

### Causal pathway: Lipid --> DM incidence, 
### controlling for HbA1c at baseline, age, gender
### Start the model with no lipid, "base model"
### Then run through a for loop to test contribution of one lipid

### Step 1: Build model (logistic)
coef.logit = rep(NA, nlipid)
pval.logit = rep(NA, nlipid)

sdata$Gender = factor(sdata$Gender, levels=c(1,2))
baseModel = glm(DM ~ Age + Gender + BMI + HbA1c + SBP + HDL + LDL + TG, family=binomial, data=sdata)
summary(baseModel)

for(k in 1:nlipid) {
  tmpModel = glm(DM ~ Age + Gender + BMI + HbA1c + SBP + HDL + LDL + TG + tmp2[,k], family=binomial, data=sdata)
  coef.logit[k] = summary(tmpModel)$coef[10,1]
  pval.logit[k] = summary(tmpModel)$coef[10,4]
  if(k %% 10 == 0) print(k)
}
hist(pval.logit, breaks=50)

# qvalue or BH
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("qvalue")
library(qvalue)
qval.logit = qvalue(pval.logit)$qvalues
plot(pval.logit, qval.logit, cex=.3, pch=19) 
abline(0,1,lty=2) ### properly modeled, it seems

tab.logit = data.frame(lipid=lipid.name, coefficient=round(coef.logit, 2),
                 pval=pval.logit, qval=qval.logit,
                 stringsAsFactors=FALSE, check.names=FALSE)
tab.logit = tab.logit[order(tab.logit$pval), ]


# volcano plot with labels for significant ones
plot(coef.logit, -log10(pval.logit), pch=19, cex=.8, col=alpha("gray", 0.5))
sid = qval.logit <= 0.2
points(coef.logit[sid], -log10(pval.logit[sid]), pch=19, cex=.8, col=alpha("red", 0.5))

tab.logit[tab.logit$qval <= 0.2, ]

#########################################################
### Step 2: Cox model with time info (time-to-event)
#install.packages("survival")
library(survival)
mid = is.na(sdata$Days)
sdata$Days[mid] = 10000
baseModelCox = coxph(Surv(Days, DM) ~ Age + Gender + BMI + HbA1c + SBP + HDL + LDL + TG, data=sdata)
summary(baseModelCox)

coef.cox = rep(NA, nlipid)
pval.cox = rep(NA, nlipid)
for(k in 1:nlipid) {
  tmpModel = coxph(Surv(Days, DM) ~ Age + Gender + BMI + HbA1c + SBP + HDL + LDL + TG + tmp2[,k], data=sdata)
  coef.cox[k] = summary(tmpModel)$coef[9,1]  ### In Cox, they don't report intercept in summary()
  pval.cox[k] = summary(tmpModel)$coef[9,5]  ### Also, need to track the columns
  if(k %% 10 == 0) print(k)
}
hist(pval.cox, breaks=100)

# qvalue or BH
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("qvalue")
library(qvalue)
qval.cox = qvalue(pval.cox)$qvalues
plot(pval.cox, qval.cox, cex=.3, pch=19) 
abline(0,1,lty=2) ### properly modeled, it seems

tab.cox = data.frame(lipid=lipid.name, coefficient=round(coef.cox, 2),
                 pval=pval.cox, qval=qval.cox,
                 stringsAsFactors=FALSE, check.names=FALSE)
tab.cox = tab.cox[order(tab.cox$pval), ]


# volcano plot with labels for significant ones
plot(coef.cox, -log10(pval.cox), pch=19, cex=.8, col=alpha("gray", 0.5))
sid = qval.cox <= 0.2
points(coef.cox[sid], -log10(pval.cox[sid]), pch=19, cex=.8, col=alpha("red", 0.5))

tab.cox[tab.cox$qval <= 0.2, ]

##### Kaplan-Meier curve

sdata$`SM d16:1/C18:0` = tmp2$`SM d16:1/C18:0`
#km1 = survfit(Surv(Days, DM) ~ (`SM d16:1/C18:0` > 0), data=sdata)
#km1.diff = survdiff(Surv(Days, DM) ~ (`SM d16:1/C18:0` > 0), data=sdata)

km1 = survfit(Surv(Days, DM) ~ (`SM d16:1/C18:0` > 0), data=sdata, subset=Age > 50)
km1.diff = survdiff(Surv(Days, DM) ~ (`SM d16:1/C18:0` > 0), data=sdata, subset=Age > 50)

plot(km1, lty=1, col=c("red","blue"), mark.time=TRUE, xlab="Days", ylab="Survival", main="SM d16:1/18:0", xlim=c(0,3000))
legend("bottomleft", c("Below zero", "Above zero"), lty=1, col=c("red","blue"))
p.val.1 <- 1 - pchisq(km1.diff$chisq, length(km1.diff$n) - 1)
text(2000, 0.8, paste("p=", round(p.val.1, 4), sep=""))


##################################
### Multivariable modeling
##################################



















