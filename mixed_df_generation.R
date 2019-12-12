library(DAPAR)
library(DAPARdata)
library(imp4p)
library(testthat)
library(MSnbase)



IntroduceMEC <- function(qData, condition, nbMEC, verbose=FALSE){
  if (isTRUE(verbose)) {
    if (nbMEC == 0){
      warning("nbMEC = 0 : No MEC have been introduced")
      return(qData)
    }
  }
  n <- nbMEC
  conds <- unique(condition)
  conds <- sample(conds,length(conds))
  for (i in 1:length(conds)){
    nSamplesInCond <- length(which(condition==conds[i]))
    nMax.MEC <- floor(n/nSamplesInCond)
    if (nMax.MEC != 0){
      ind <- sample(nrow(qData), sample(nMax.MEC,1))
      qData[ind,which(condition==conds[i])] <- NA
      n <- n - length(ind)*nSamplesInCond
    }
  }
  return(qData)
}

IntroducePOV <- function(qData, condition, nbPOV, verbose=FALSE){
  if (isTRUE(verbose)) {
    if (nbPOV == 0){
      warning("nbPOV = 0 : No POV have been introduced")
      return(qData)
    }
  }
  n <- nbPOV
  condition <- sample(colnames(qData), ncol(qData))
  for (i in 1:length(condition)){
    if (n != 0){
      ind.autorises<- which(!is.na(qData)[,i])
      ind <- sample(ind.autorises, sample(n, 1))
      qData[ind,which(condition==colnames(qData)[i])] <- NA
      n <- n - length(ind)
    }
  }
  return(qData)
}




##' This function is xxxxxx
##'
##' @title xxxxxx
##' @param xxx.
##' @param sTab xxxx 
##' @param xxxxx
##' @param type xxxxx
##' @return A list of two items : xxxxx
##' @author Enora Fremy, Samuel Wieczorek
##' @examples
##' ll <- GenerateRandomDataset(nbCond=3, nRep=3, mismatch.nRep =TRUE, prop.MV=0)
GenerateRandomDataset <- function(params){
  
  qData <- pData <- NULL
  for(i in 1:length(params)){
    assign(names(params[i]), params[[i]])    
  }
  
  if (prop.MV==0 && (1 != (prop.MEC + prop.POV))){
    warning("The sum of probability of POV and MEC missing values must be equal to 1.")
    return(NULL)
  }
  
  if (isTRUE(fixedDesign)){
    nbCond <- nbCond
    nRep <- nRep
  } else {
    nbCond = sample(c(2:nbCond),1)
    nRep = sample(c(2:nRep),1)
  }
  
  
  #creation dataset sans missing values
  base <- LETTERS[1:nbCond]
  FC.factor = 100
  if (!mismatch.nRep) {
    sample.names <- unlist(lapply(base, function(x)paste0(x,"_",1:nRep)))
    
  } else {
    sample.names <- unlist(lapply(base, function(x)paste0(x,"_",1:sample(1:max(nRep),1 ))))
  }
  
  
  condition <- as.vector(sapply(sample.names, function(x) unlist(strsplit(x, "_"))[1]))
  nb.samples <- length(sample.names)
  qData <-  data.frame(matrix(rep(0,size * nb.samples), ncol = nb.samples))
  colnames(qData) <- sample.names
  
  for (i in 1:length(condition)){
    IndiceCurrentCond <- which (base==condition[i])
    qData[,i] <-rnorm(size, FC.factor*IndiceCurrentCond,0.1)
  }
  
  # creation du pData
  pData <- data.frame(Sample.name=sample.names,
                      Condition =condition,
                      Bio.Rep =1:ncol(qData))
  
  
  rownames(pData) <- pData$Sample.name
  
  
  nb.MV <- floor(size * prop.MV)
  print(nb.MV)
  nb.MEC <- floor(nb.MV * prop.MEC)
  nb.POV <- floor(nb.MV * prop.POV)
  
  # introduction des lignes vides
  # TODO
  
  #introduction des MEC
  print("Introdution of MEC MV")
  qData <- IntroduceMEC(qData, condition, nb.MEC)
  
  #introduction des POV
  print("Introdution of POV MV")
  qData <- IntroducePOV(qData, condition, nb.POV)
  
  return (list(qData=qData, pData=pData))
}



##' This function is xxxxxx
##'
##' @title xxxxxx
##' @param ll
##' @param nCond xxxx 
##' @param xxxxx
##' @param type xxxxx
##' @return A list of two items : xxxxx
##' @author Enora Fremy
##' @examples
##' ll <- GenerateRandomDataset(nbCond=3, nRep=3, mismatch.nRep =TRUE, prop.MV=0)
##' ll.mixed <- mix_dataset_Enora(ll, do.interC = TRUE, do.intraC = TRUE, do.fullRandom = TRUE)
mix_dataset_Enora <- function(ll, params){
  
  # print(params)
  # print(unlist(params))
  for(i in 1:length(params)){
    assign(names(params[i]), params[[i]])    
  }
  # print(do.interC)
  # print(do.intraC)
  # print(do.fullRandom)
  # 
  
  
  qData <- ll$qData
  pData <- ll$pData
  nCond <- length(unique(pData$Condition))
  interC <- intraC <- fullRandom <- 0
  if (isTRUE(do.interC)){interC <- sample(c(TRUE, FALSE), 1)}
  if (isTRUE(do.intraC)){intraC <- sample(c(TRUE, FALSE), 1)}
  if (isTRUE(do.fullRandom)){fullRandom <- sample(c(TRUE, FALSE), 1)}
  #### Mix columns qData ####
  
  if (fullRandom == 0) {
    
    if (interC == 1 && intraC == 0) { 
      
      print("conditions shuffled, replicates unchanged")
      interC.list <- c(1:ncol(qData))
      interC.list <- split(interC.list, sort(interC.list%%nCond))
      new.interC.list <- unlist(sample(interC.list,nCond), use.names = F)
      qData <- qData[,new.interC.list]
      pData <- pData[new.interC.list,]
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
      pData <- pData[new.order,]
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
      pData <- pData[new.order,]
      new.interC.list <- unlist(sample(er.ac,nCond), use.names = F)
      qData <- qData[,new.interC.list]
      pData <- pData[new.interC.list,]
      
    }
    
  }
  
  else { ## Full random option
    random_col_indices <- sample(ncol(qData),ncol(qData))
    qData <- qData[,random_col_indices]
    pData <- pData[random_col_indices,]
  }
  
  return(list(qData=qData, pData=pData))
}



##' This function is xxxxxx
##'
##' @title xxxxxx
##' @param ll
##' @param type xxxxx
##' @return An object of class MSnSet
##' @author Enora Fremy, Samuel Wieczorek
##' @examples
##' 
CreateMinimalistMSnset <- function(ll){
  qData <- ll$qData
  pData <- ll$pData
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


## Noms et parametres des fonctions d'imputation a tester
## utilises par do.call
##' This function is xxxxxx
##'
##' @title xxxxxx
##' @param xxx.
##' @param sTab xxxx 
##' @param xxxxx
##' @param type xxxxx
##' @return A list of two items : xxxxx
##' @author Enora Fremy, Samuel Wieczorek
##' @examples
GetListFuncs <- function(obj=NULL){
  
  ll <- NULL
  if (is.null(obj)) {
    ##liste des fonctions a tester
    ll <- c("wrapper.dapar.impute.mi",
            "wrapper.impute.mle",
            "wrapper.impute.slsa",
            "wrapper.impute.detQuant",
            "wrapper.impute.pa",
            "wrapper.impute.fixedValue",
            "wrapper.impute.KNN"
            
    ) }
  else {
    ll <- list(list(obj,nb.iter = 1, progress.bar = FALSE),
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
##' This function is xxxxxx
##'
##' @title xxxxxx
##' @param obj.original xxxx
##' @param obj.mixed xxxx 
##' @return A list of two items : xxxxx
##' @author Enora Fremy, Samuel Wieczorek
##' @examples
##' genDatasetArgs <- list(nbCond=3, nRep=3, mismatch.nRep=FALSE, prop.MV = 0.2, size = 100)
##' mixDatasetArgs <- list(do.interC=TRUE, do.intraC=TRUE, do.fullRandom=TRUE)
##' res <- test_imputation(nTest = 5,genDatasetArgs , mixDatasetArgs)
##' ll.original <- do.call("GenerateRandomDataset", genDatasetArgs)
##' ll.mixed <- do.call("mix_dataset_Enora",list(ll.original,mixDatasetArgs))
##' obj.original <- CreateMinimalistMSnset(ll.original)
##' obj.mixed <- CreateMinimalistMSnset(ll.mixed)
##' test_impute_functions(obj.original, obj.mixed)
test_impute_functions <- function(obj.original, obj.mixed){
  
  FUN <- GetListFuncs() # Liste des fonctions a tester
  df_list <- list()
  df_list$obj.original <- obj.original
  df_list$obj.mixed <- obj.mixed
  ll_params_original <- GetListFuncs(obj.original) # param a charger pour la fonction courante
  ll_params_mixed <- GetListFuncs(obj.mixed)
  
  # Pour chaque fonction d'imputation a tester
  for (i in 1:length(FUN)){
    #i=1
    cat("\nTest de la fonction : ",FUN[i])
    #recupere les fonctions a tester et leurs parametres
    
    
    tryCatch(
      {# execution de la fonction d'imputation a tester
        obj.original.imputed <- obj.original
        tmp.original.imputed <- do.call(FUN[i],ll_params_original[[i]])
        if (class(tmp.original.imputed)[1] != 'MSnSet'){
          Biobase::exprs(obj.original.imputed) <- tmp.original.imputed
        } else {
          obj.original.imputed <- tmp.original.imputed
        }
        
        obj.mixed.imputed <- obj.mixed
        tmp.mixed.imputed <- do.call(FUN[i],ll_params_mixed[[i]])
        if (class(tmp.mixed.imputed)[1] != 'MSnSet'){
          Biobase::exprs(obj.mixed.imputed) <- tmp.mixed.imputed
        } else {
          obj.mixed.imputed <- tmp.mixed.imputed
        }
        
        # tri des colonnes du dataset melange suivant l'ordre des 
        # colonnes du dataset original
        original.order <- colnames(Biobase::exprs(obj.original.imputed))
        tmp.qData <- Biobase::exprs(obj.mixed.imputed)[,original.order]
        tmp.pData <- Biobase::pData(obj.mixed.imputed)[original.order,]
        obj.mixed.imputed <- CreateMinimalistMSnset(list(qData=tmp.qData, pData=tmp.pData))

        # stocker les resultats d'imputation de orginal et mixed pour chaque methode
        
        method = paste0("obj.ori.",FUN[i])
        df_list[[method]] <- obj.original.imputed
        method = paste0("obj.mix.",FUN[i])
        df_list[[method]] <- obj.mixed.imputed
        
        # test de comparaison
        expect_equal(Biobase::exprs(obj.original.imputed), Biobase::exprs(obj.mixed.imputed), tolerance=1)
        
      },
      warning = function(w) {
        #message(w)
        },
      error = function(e) {message(e)},
      finally = {}
    )
    
  }
  return(df_list)
}




#------------------------------------------------------------
# Automatisation
#------------------------------------------------------------
##' This function is xxxxxx
##'
##' @title xxxxxx
##' @param xxx.
##' @param sTab xxxx 
##' @param xxxxx
##' @param type xxxxx
##' @return A list of two items : xxxxx
##' @author Enora Fremy, Samuel Wieczorek
##' @examples
##' genDatasetArgs <- list(nbCond=3, nRep=3, fixedDesign = FALSE,mismatch.nRep=TRUE, prop.MV = 0.2,  prop.MEC=0.3, prop.POV=0.7,size = 100)
##' mixDatasetArgs <- list(do.interC=TRUE, do.intraC=TRUE, do.fullRandom=TRUE)
##' res <- test_imputation(nTest = 20,genDatasetArgs , mixDatasetArgs)
test_imputation <- function(nTest = 5, genDatasetArgs, mixDatasetArgs) {
  res <- list()
  
  for (i in 1:nTest) {
    
    print(paste0("#######--------------------------- TEST ", i, " ---------------------------------#######"))
    
    
    
    # 1 - generation qData et pData
    ll.original <- do.call("GenerateRandomDataset", list(genDatasetArgs))
    print(unlist(genDatasetArgs))
    # 2 - etape de melange
    ll.mixed <- do.call("mix_dataset_Enora",list(ll.original,mixDatasetArgs))
    print(unlist(mixDatasetArgs))
    # 3 - mise sous forme de MSnset
    obj.original <- CreateMinimalistMSnset(ll.original)
    obj.mixed <- CreateMinimalistMSnset(ll.mixed)
    
    # 4 -test en serie des fonctions d'imputation
    impute <- test_impute_functions(obj.original, obj.mixed)
    tour <- paste0("tour",i)
    res[[tour]] <- impute
    
  }
  return(res)
}



showDatasets <- function(res, tour, size){
  
  data <- res[[paste0('tour', tour)]]
  for (i in 1:length(names(data))){
    print(paste0('---------',names(data)[i] ,'---------' ))
    print(exprs(data[[i]])[1:size,])
  }
  
}

genDatasetArgs <- list(nbCond=3, nRep=3, fixedDesign = FALSE,mismatch.nRep=FALSE, prop.MV = 0.2,  prop.MEC=0.3, prop.POV=0.7,size = 100)
mixDatasetArgs <- list(do.interC=FALSE, do.intraC=FALSE, do.fullRandom=TRUE)
res <- test_imputation(nTest = 5,genDatasetArgs , mixDatasetArgs)

showDatasets(res,2,5)
