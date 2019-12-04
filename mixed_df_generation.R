library(DAPAR)
library(DAPARdata)
library(imp4p)
library(testthat)

df_generation <- function(qData, pData, nCond, nRep, mismatch.nRep = FALSE) { 
  
  # Order of Sample.name of pData and qdata must be the same
  if ( table(colnames(qData) == pData$Sample.name) ) {
    
    
    # pData for ncol and colnames
    pData.plus <- data.frame(matrix(nrow = nRep*nCond, ncol = ncol(pData)))
    colnames(pData.plus) <- colnames(pData)
    pData <- pData.plus ; rm(pData.plus)
    pData$Bio.Rep <- c( 1:nrow(pData) )
    
    # set Sample.Name to A1,...,A(nRep),B1,...,B(nRep), C1,... 
    # and Condition to A,B,C,D...
    Sample.name <- vector()
    for (i in 1:nCond) {
      #i=4
      Sample.name <- c(Sample.name, paste0(LETTERS[i],1:nRep))
    }
    pData$Sample.name <- Sample.name
    pData$Condition <- substr(pData$Sample.name,1,1)
    rownames(pData) <- pData$Sample.name
    
    # qData
    qData.plus <- data.frame(matrix(nrow = nrow(qData), ncol = nRep*nCond))
    rownames(qData.plus) <- rownames(qData) ## ! rownames(qData) Begin by 0 ! ##
    qData.plus[,1:ncol(qData)] <- qData
    
    
    
    
    
    
    colnames(qData.plus) <- pData$Sample.name
  }
  res <- list(pData=pData,qData=qData.plus)
  return(res)
  
}

#------------------------------------------------------------
data("Exp1_R25_pept")
#nCond = 3
#nRep = 3
qData <- Biobase::exprs(Exp1_R25_pept)
pData <- Biobase::pData(Exp1_R25_pept)

#------------------------------------------------------------
res <- df_generation(qData, pData, 3, 3)
