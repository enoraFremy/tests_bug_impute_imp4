library(DAPAR)
library(DAPARdata)
library(imp4p)
library(testthat)

df_generation_Sam <- function(nbCond, nRep, mismatch.nRep = FALSE,prop.MV = 0.2, size = 100){
  qData <- pData <- NULL
  
  
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
  indices.MV <- sample(ncol(qData)*size, nb.MV)
  qData <- unlist(qData)
  qData[indices.MV] <- NA
  qData <- data.frame(matrix(qData, ncol=nb.samples))
  colnames(qData) <-sample.names 
  
  # introduction des lignes vides
  # TODO
  
  #introduction des MEC
  #qData <- IntroduceMEC(qData, nbMEC)
  # nb.MEC <- floor(nb.MV * prop.MEC)
  # indices.MEC <- sample(nbCond*size, nb.MEC)
  # for (i in 1:length(indices.MEC)){
  #   ind.cond <-floor(indices.MEC[i]/size) +1
  #   indices.samples <- which (base[ind.cond]==condition)
  #   qData[indices.MEC[i]%%size,indices.samples] <- NA
  # }
  
  
  #introduction des POV
  #qData <- IntroducePOV(qData, nbPOV)
  # nb.POV <- floor(nb.MV * (1-prop.MEC))
  # ensDeTest <- setdiff(nbCond*size, indices.MEC)
  # for (i in 1:length(ensDeTest)){
  #   ind.cond <-floor(indices.MEC[i]/size) +1
  #   indices.samples <- which (base[ind.cond]==condition)
  #   qData[ensDeTest[i],indices.samples] <- NA
  # }
  # 
  
  return (list(qData=qData, pData=pData))
}


mix_dataset_Enora <- function(ll, nCond, nRep, mismatch.nRep = FALSE, 
                              interC = FALSE, intraC = FALSE, fullRandom = FALSE){
  
  qData <- ll$qData
  pData <- ll$pData
  
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
  
  else {
    
    print("full random")
    random_col_indices <- sample(ncol(qData),ncol(qData))
    #print(paste0("random_col_indices: ",list(random_col_indices)))
    qData <- qData[,random_col_indices]
    pData <- pData[random_col_indices,]
  }
  
  
  #res$qData.mixed <- qData
  return(list(qData=qData, pData=pData))
}


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
  
  FUN <- GetListFuncs() # Liste des fonctions a tester
  df_list <- list()
  
  # Pour chaque fonction d'imputation a tester
  for (i in 1:length(FUN)){
    #i=1
    cat("\nTest de la fonction : ",FUN[i])
    #recupere les fonctions a tester et leurs parametres
    ll_params_original <- GetListFuncs(obj.original) # param a charger pour la fonction courante
    ll_params_mixed <- GetListFuncs(obj.mixed)
    
    df_list$res.original <- exprs(obj.original)
    df_list$res.mixed <- exprs(obj.mixed)
    
    tryCatch(
      {# execution de la fonction d'imputation a tester
        obj.original.imputed <- do.call(FUN[i],ll_params_original[[i]])
        obj.original.mixed <- do.call(FUN[i],ll_params_mixed[[i]])
        
        # tri des colonnes du dataset melange suivant l'ordre des 
        # colonnes du dataset original
        original.order <- colnames(exprs(obj.original.imputed))
        qData.original.mixed <- exprs(obj.original.mixed)
        qData.original.mixed <- (exprs(obj.original.mixed))[,original.order]
        
        # stocker les resultats d'imputation de orginal et mixed pour chaque methode
        
        method = paste0("res.ori.",FUN[i])
        df_list[[method]] <- exprs(obj.original.imputed)
        method = paste0("res.mix.",FUN[i])
        df_list[[method]] <- qData.original.mixed
        # exprs(obj.original.mixed) ne se modifie pas
        
        # head(exprs(obj.original))
        # head(exprs(obj.mixed))
        # head(exprs(obj.original.imputed))
        # #head(exprs(obj.original.mixed))
        # head(qData.original.mixed)
        
        # test de comparaison
        expect_equal(exprs(obj.original.imputed), qData.original.mixed,tolerance=1)

      },
      warning = function(w) {message(w)},
      error = function(e) {message(e)},
      finally = {}
    )
    
  }
  return(df_list)
}




#------------------------------------------------------------
# Automatisation
#------------------------------------------------------------
test_imputation <- function(nbCond, nRep, size = 100, mismatch.nRep = FALSE, interC, intraC, fullRandom, nTest = 5) {
  
  res <- list()

    
  for (i in 1:nTest) {
    
    cat("\n Caracteristiques du dataset :\n
        *** nCond: ",nbCond, ", nRep: ", nRep, ", interC: ", interC, ", intraC: ",intraC, ", fullRandom: ",fullRandom, " ***\n")
    
    # 1 - generation qData et pData
    res.original <- df_generation_Sam(nbCond, nRep, mismatch.nRep, prop.MV = 0.2, size = size)
    
    # 2 - etape de melange
    res.mixed <- mix_dataset_Enora(res.original, nbCond, nRep, mismatch.nRep = FALSE, 
                                   interC, intraC, fullRandom)
    #print(head(res.original$qData))
    #print(head(res.mixed$qData))
    
    # 3 - mise sous forme de MSnset
    obj.original <- CreateMinimalistMSnset(res.original)
    obj.mixed <- CreateMinimalistMSnset(res.mixed)
    
    # 4 -test en serie des fonctions d'imputation
    impute <- test_impute_functions(obj.original, obj.mixed)
    tour <- paste0("tour",i)
    res[[tour]] <- impute

  }
  return(res)
}

res <- test_imputation(3,2,100,F,0,1,0,5)
summary(res)
head(res$tour1$res.original)
head(res$tour1$res.mixed)
head(res$tour1$res.ori.wrapper.impute.mle)
head(res$tour1$res.mix.wrapper.impute.mle)


# nCond = sample(c(2:5),1)
# nRep = sample(c(2:4),1)
# interC = sample(c(0,1),1)
# intraC = sample(c(0,1),1)
# fullRandom = sample(c(0,1),1)
