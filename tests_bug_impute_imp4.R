library(DAPAR)
library(DAPARdata)


#### data loading ####
data("Exp1_R25_prot")
qData <- Biobase::exprs(Exp1_R25_prot)
#dim(qData) # 2384 6
#head(qData)

#### data pre ####
qData <- cbind(qData, qData[,1:3]) # add a third condition
#head(qData)
colnames(qData) <- c("A1","A2","A3","B1","B2","B3","C1","C2","C3")
#head(qData)

qData[,1:3] <- apply(qData[,1:3], 2, function(x) rnorm(x, mean = 1, sd = 1))
qData[,4:6] <- apply(qData[,4:6], 2, function(x) rnorm(x, mean = 500, sd = 1))
qData[,7:9] <- apply(qData[,7:9], 2, function(x) rnorm(x, mean = 1000, sd = 1))
# different rnorm per condition

#summary(qData)
#boxplot(qData)
# three different rnorm per condition OK

#### imputation ####
# wrapper.dapar.impute.mi
# wrapper.impute.slsa

