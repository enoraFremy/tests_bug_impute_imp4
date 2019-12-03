library(DAPAR)
library(DAPARdata)
library(imp4p)

#### Function to test the egality between the original (ordered) dataset and
#### the mixed one which has been imputed
testSpecialDatasets <- function(qData.original, qData.mixed){
  
  original.order <- colnames(qData.original)
  qData.mixed.reordered <- qData.mixed[,original.order]
  
  testthat::expect_equal(qData.mixed.reordered, qData.original,tolerance=1)
}



#### data loading ####
data("Exp1_R25_prot")
dat <- mvFilter(Exp1_R25_prot, type="allCond", th = 1)
qData <- Biobase::exprs(dat)

#dim(qData) # 2384 6
#head(qData)

#### data pre ####
qData <- cbind(qData, qData[,c(1,3,5)]) # add a third condition
#head(qData)
colnames(qData) <- c("A1","A2","A3","B1","B2","B3","C1","C2","C3")
#head(qData)

qData.NA <- is.na(qData) # keep NA for imputation

qData[,1:3] <- apply(qData[,1:3], 2, function(x) rnorm(x, mean = 1, sd = 1))
qData[,4:6] <- apply(qData[,4:6], 2, function(x) rnorm(x, mean = 500, sd = 1))
qData[,7:9] <- apply(qData[,7:9], 2, function(x) rnorm(x, mean = 1000, sd = 1))

qData[qData.NA] <- NA
rm(qData.NA)
# different rnorm per condition

#summary(qData)
#boxplot(qData)
# three different rnorm per condition OK

#### imputation ####

# wrapper.impute.slsa
conds <- as.factor(c("A","A","A","B","B","B","C","C","C"))
res.slsa <- impute.slsa(qData, conditions=conds, nknn=15, selec="all", weight=1, ind.comp=1)

boxplot(qData)
boxplot(res.slsa)

# wrapper.dapar.impute.mi
conditions <- as.factor(c("A","A","A","B","B","B","C","C","C"))
repbio <- as.factor(c("1","2","3","4","5","6","7","8","9"))
reptech <- as.factor(NULL)

tab <- qData

dat.slsa = impute.rand(tab = tab, conditions = conditions)
res = estim.mix(tab = tab, tab.imp = dat.slsa, conditions = conditions, 
                x.step.mod = 300, 
                x.step.pi = 300, nb.rei = 100)
# Error in optim(init, fr, gr = NULL, x = (abs[, j] - xmin), pi_est = pi_init,  : 
# L-BFGS-B nécessite des valeurs finies de 'fn'


res.mi <- impute.mi(qData, conditions=conds, nknn=15, selec="all", weight=1, ind.comp=1)
# Error in optim(init, fr, gr = NULL, x = (abs[, j] - xmin), pi_est = pi_init,  : 
# L-BFGS-B nécessite des valeurs finies de 'fn'

boxplot(qData)
boxplot(res.mi)
