library(DAPAR)
library(DAPARdata)
library(imp4p)
library(testthat)

CreateMinimalistMSnset <- function(qData,pData){
  Intensity <- matrix(as.numeric(gsub(",", ".",as.matrix(qData )))
                      , ncol=ncol(qData)
                      , byrow=FALSE)
  
  colnames(Intensity) <- colnames(qData)
  rownames(Intensity) <- rownames(qData)
  
  ##building fData of MSnSet file
  fd <- data.frame( qData,stringsAsFactors = FALSE)
  
  pd <- as.data.frame(pData,stringsAsFactors = FALSE)
  
  ##Integrity tests
  if(identical(rownames(Intensity), rownames(fd))==FALSE)
    stop("Problem consistency between
             row names expression data and featureData")
  
  if(identical(colnames(Intensity), rownames(pd))==FALSE) 
    stop("Problem consistency between column names 
             in expression data and row names in phenoData")
  
  obj <- MSnSet(exprs = Intensity, fData = fd, pData = pd)
  
  return(obj)
}



df_generation <- function(qData, pData, nCond, nRep, mismatch.nRep = FALSE, interC = 0, intraC = 0, fullRandom = 0) { 
  
  # Order of Sample.name of pData and qdata must be the same
  if (table(colnames(qData) == pData$Sample.name)) {
    
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
    qData.plus[,1:ncol(qData)] <- qData
    
    for (i in (ncol(qData)+1):ncol(qData.plus)) {
      random_col_qData <- sample(ncol(qData),1)
      print(paste0("Column qData random: ", random_col_qData))
      qData.plus[,i] <- qData[,random_col_qData]
    }
    
    colnames(qData.plus) <- pData$Sample.name
    qData <- as.data.frame(qData.plus)
    
    
    #### Mix columns qData ####
    
    if (fullRandom == 0) {
      
      if (interC == 1 && intraC == 0) { 
        
        print("conditions shuffled, replicates unchanged")
        interC.list <- c(1:ncol(qData))
        interC.list <- split(interC.list, sort(interC.list%%nCond))
        new.interC.list <- unlist(sample(interC.list,nCond), use.names = F)
        qData <- qData[,new.interC.list]
        
      }
      
      if (interC == 0 && intraC == 1) { 
        
        print("conditions unchanged, replicates shuffled")
        intraC.list <- c(1:ncol(qData))
        intraC.list <- split(intraC.list, sort(intraC.list%%nCond))
        
        new.order <- vector()
        
        for (i in (1:nCond) ) {
          #print(i)
          #i=1
          new.order <- c(new.order,sample(intraC.list[[i]]))
          
          
        }
        qData <- qData[,new.order]
      }
      
      if (interC == 1 && intraC == 1) { 
        
        print("conditions and replicates shuffled")
        er.ac <- c(1:ncol(qData))
        er.ac <- split(er.ac, sort(er.ac%%nCond))
        new.order <- vector()
        
        for (i in (1:nCond) ) {
          #print(i)
          #i=1
          new.order <- c(new.order,sample(er.ac[[i]]))
          
        }
        qData <- qData[,new.order]
        new.interC.list <- unlist(sample(er.ac,nCond), use.names = F)
        qData <- qData[,new.interC.list]
        
        
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
# nCond = 3
# nRep = 3
# interC = 0
# intraC = 1
# fullRandom = 0
qData <- Biobase::exprs(Exp1_R25_pept)
pData <- Biobase::pData(Exp1_R25_pept)

#------------------------------------------------------------
res <- df_generation(qData, pData, nCond = 3, nRep = 3, mismatch.nRep = FALSE, interC = 0, intraC = 1, fullRandom = 0)
#View(res$pData)
View(res$qData)
new_MSnset <- CreateMinimalistMSnset(qData, pData) # pour les tests (d'imputation imp4p)