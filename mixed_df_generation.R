library(DAPAR)
library(DAPARdata)
library(imp4p)
library(testthat)


##########################################################
##
mix_dataset <- function(ll=NULL){
  
  if(is.null(ll)){warning('The list in empty')
    return(NULL)}
  
  qData <- ll$qData
  pData <- ll$pData
  
  #### Mix columns qData ####
  random_col_indices <- NULL
  if (fullRandom == FALSE) {
    
    if (interC == TRUE && intraC == FALSE) { 
      
      print("conditions shuffled, replicates unchanged")
      # aller par multiples de nRep
      interC.list <- list()
      ll <- lapply(unique(conditions), function(x) which(conditions==x))
      ll <- ll[sample(nbCond,nbCond)]
    }
    
    if (interC == FALSE && intraC == TRUE) { 
      
      print("conditions unchanged, replicates shuffled")
      ll <- lapply(unique(conditions), function(x) which(conditions==x))
      ll <- lapply(ll, function(x) sample(x,length(x)))
    }
    
    if (interC == TRUE && intraC == TRUE) { 
      
      print("conditions and replicates shuffled")
      ll <- lapply(unique(conditions), function(x) which(conditions==x))
      ll <- lapply(ll, function(x) sample(x,length(x)))
      ll <- ll[sample(nbCond,nbCond)]
    }
    
  }
  
  else {
    
    print("full random")
    random_col_indices <- sample(ncol(qData),ncol(qData))
    #print(paste0("random_col_indices: ",list(random_col_indices)))
  }
  qData <- qData[,random_col_indices]
  pData <- pData[random_col_indices,]

res <- list(pData=pData,qData=qData)
return(res)
}

##########################################################
##
df_generation <- function(qData, pData, nCond, nRep, mismatch.nRep = FALSE, 
                          interC = FALSE, intraC = FALSE, fullRandom = FALSE) { 
  
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
df.original <- df_generation(qData, pData, nCond = 3, nRep = 3)
df.mixed <- mix_dataset(df.original)
View(res$pData)
View(res$qData)





## Noms et parametres des fonctions d'imputation à tester
## utilises par do.call
GetListFuncs <- function(obj=NULL){
  
  ll <- NULL
  if (is.null(obj)) {
    ##liste des fonctions à tester
  ll <- c("wrapper.dapar.impute.mi",
                      "wrapper.impute.mle",
                      "wrapper.impute.slsa",
                      "wrapper.impute.detQuant",
                      "wrapper.impute.pa",
                      "wrapper.impute.fixedValue",
                      "wrapper.impute.KNN"
  ) }
  else {
  ll <- list(list(obj,nb.iter = 1),
             list(obj),
             list(obj),
             list(obj,qval=0.025, factor=1),
             list(obj,q.min = 0.025),
             list(obj,fixVal=0),
             list(obj,K=10)
             )
  }
  return(ll)
}


##Teste les fonctions d'imputation sur 2 datasets
test_impute_functions <- function(obj.original, obj.mixed){
  FUN <- GetListFuncs()
  # Pour chaque fonction d'imputation à tester
  for (i in 1:length(FUN)){
    print(paste0("test de la fonction : ",FUN[i]))
    #recupere les fonctions a tester et leurs parametres
    ll_params_original <- GetListFuncs(obj.original)
    ll_params_mixed <- GetListFuncs(obj.mixed)
    
    # execution de la fonction d'imputation a tester
    obj.original.imputed <- do.call(FUN[i],ll_params_original[[i]])
    obj.original.mixed <- do.call(FUN[i],ll_params_mixed[[i]])
    
    # tri des colonnes du dataset melange suivant l'ordre des 
    # colonnes du dataset original
    original.order <- colnames(exprs(obj.original.imputed))
    exprs(obj.original.mixed) <- exprs(obj.original.mixed)[,original.order]
    
    
    # test de comparaison
    expect_equal(exprs(obj.original.imputed), exprs(obj.original.mixed),tolerance=1)
  }
}
