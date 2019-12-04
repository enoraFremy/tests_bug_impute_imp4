library(DAPAR)
library(DAPARdata)
library(imp4p)
library(testthat)

df_generation <- function(pData, qData, nCond, nRep, mismatch.nRep) {
  
  # pData for ncol and colnames
  pData.plus <- data.frame(matrix(nrow = nRep*nCond, ncol = ncol(pData)))
  colnames(pData.plus) <- names(pData)
  pData <- pData.plus ; rm(pData.plus)
  pData$Bio.Rep <- c( 1:nrow(pData) )
  
  # replace Sample.Name by A1,...,A(nRep),B1,...,B(nRep), C1,... 
  # and Condition by A,B,C,D...
  Sample.name <- vector()
  for (i in 1:nCond) {
    #i=4
    Sample.name <- c(Sample.name, paste0(LETTERS[i],1:nRep))
  }
  pData$Sample.name <- Sample.name
  pData$Condition <- substr(pData$Sample.name,1,1)
  rownames(pData) <- pData$Sample.name

  # qData
  
  
}

data("Exp1_R25_pept")
nCond = 3
nRep = 3
qData <- Biobase::exprs(Exp1_R25_pept)
pData <- Biobase::pData(Exp1_R25_pept)

#------------------------------------------------------------
df_generation(qData, pData, 3, 3, FALSE)


qData <- cbind(qData, qData[,c(1,3,5)]) # Third condition
colnames(qData) <- c("A1","A2","A3","B1","B2","B3","C1","C2","C3")
pData <- rbind(pData, pData[1:nbCond,])
pData$Condition <- c(rep("A",nbCond), rep("B", nbCond), rep("C", nbCond))
pData$Sample.name <- paste0(pData$Condition, 1:nbCond)
rownames(pData) <- pData$Sample.name
