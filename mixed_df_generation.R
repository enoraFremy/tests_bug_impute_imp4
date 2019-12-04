library(DAPAR)
library(DAPARdata)
library(imp4p)
library(testthat)

df_generation <- function(qData, pData, nCond, nRep, mismatch.nRep = FALSE, interC = 0, intraC = 0, fullRandom = 0) { 
  
  # Order of Sample.name of pData and qdata must be the same
  if ( table(colnames(qData) == pData$Sample.name) ) {
    
    #### dataset base ####
    # pData 
    pData.plus <- data.frame(matrix(nrow = nRep*nCond, ncol = ncol(pData)))
    colnames(pData.plus) <- colnames(pData)
    pData <- pData.plus ; rm(pData.plus)
    pData$Bio.Rep <- c( 1:nrow(pData) )
    
    Sample.name <- vector()# set Sample.Name to A1,...,A(nRep),B1,...,B(nRep), C1,... and Condition to A,B,C,D...
    for (i in 1:nCond) {
      #i=4
      Sample.name <- c(Sample.name, paste0(LETTERS[i],1:nRep))
    }
    pData$Sample.name <- Sample.name
    pData$Condition <- substr(pData$Sample.name,1,1)
    rownames(pData) <- pData$Sample.name
    
    # qData
    qData.plus <- data.frame(matrix(nrow = nrow(qData), ncol = nRep*nCond))
    rownames(qData.plus) <- rownames(qData) ## rownames(qData) must begin by 0 ! ##
    colnames(qData.plus) <- pData$Sample.name
    if (nRep == 3) {
      qData.plus[,1:ncol(qData)] <- qData
    }
    else { print("gerer les replicats != 3 quand qData pasted dans qData.plus")  }
    
    for (i in (ncol(qData)+1):ncol(qData.plus)) {
      random_col_qData <- sample(ncol(qData),1)
      print(paste0("Column qData random: ", random_col_qData))
      qData.plus[,i] <- qData[,random_col_qData]
    }
    
    colnames(qData.plus) <- pData$Sample.name
    qData <- qData.plus
    
    
    #### Mix columns qData ####
    
    if (fullRandom == 0) {
      
      if (interC == 1 && intraC == 0) { 
        
        print("conditions shuffled, replicates unchanged")
        # aller par multiples de nRep
        interC.list <- list()
        
        
      }
      
      if (interC == 0 && intraC == 1) { 
        
        print("conditions unchanged, replicates shuffled")
      }
      
      if (interC == 1 && intraC == 1) { 
        
        print("conditions and replicates shuffled")
      }
      
    }
    
    else {
      
      print("full random")
      random_col_indices <- sample(ncol(qData),ncol(qData))
      #print(paste0("random_col_indices: ",list(random_col_indices)))
      qData <- qData[,random_col_indices]
    }
    
  }
  
  res <- list(pData=pData,qData=qData)
  return(res)
  
}

#------------------------------------------------------------
data("Exp1_R25_pept")
nCond = 3
nRep = 3
interC = 1
intraC = 0
fullRandom = 0
qData <- Biobase::exprs(Exp1_R25_pept)
pData <- Biobase::pData(Exp1_R25_pept)

#------------------------------------------------------------
res <- df_generation(qData, pData, nCond = 3, nRep = 3)
View(res$pData)
View(res$qData)
