#library(installr)
#uninstall.packages("DAPAR")
#uninstall.packages("imp4p")

#install.packages("~/Github/DAPAR_1.19.13.tar.gz", repo = NULL, type = "source")
#install.packages("~/R/imp4p_0.9.tar.gz", repo = NULL, type = "source")

library(DAPAR)
library(DAPARdata)
library(imp4p)
library(testthat)

## --------------------------------------------------------------- ##
# correction codes mle qd nb de rep different par condition. Attention à levels
# wrapper.impute.mle.toto <- function(obj){
#   cond <- as.factor(Biobase::pData(obj)$Condition)
#   
#   res <- impute.mle.toto(Biobase::exprs(obj), conditions=cond)
#   
#   Biobase::exprs(obj) <-res
#   return (obj)
# }
# 
# 
# impute.mle.toto <- function (tab, conditions) 
# {
#   tab_imp = as.matrix(tab)
#   nb_cond = length(levels(conditions))
#   nb_rep = rep(0, nb_cond)
#   k = 1
#   for (n in 1:nb_cond) {
#     nb_rep[n] = sum((conditions == unique(conditions)[n]))
#     #xincomplete = as.matrix(tab[, (k:(k + nb_rep[n] - 1))])
#     indices <- which((conditions == unique(conditions)[n]))
#     xincomplete = as.matrix(tab[, indices])
#     nbna = apply(xincomplete, 1, function(x) {
#       sum(is.na(x))
#     })
#     xincomplete1 = xincomplete[which(nbna != nb_rep[n]), 
#                                ]
#     s <- prelim.norm(xincomplete1)
#     thetahat <- em.norm(s, showits = FALSE)
#     rngseed(1234567)
#     xcomplete1 <- imp.norm(s, thetahat, xincomplete1)
#     tab_imp[which(nbna != nb_rep[n]), indices] = xcomplete1
#     # k = k + nb_rep[n]
#   }
#   return(tab_imp)
# }
## --------------------------------------------------------------- ##


#------------------------------------------------------------
# Fonctions
#------------------------------------------------------------

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
##' ll <- GenerateRandomDataset(nbCond=3, nRep=3, mismatch.nRep =TRUE, prop.MV=0.2,size=1000)
GenerateRandomDataset <- function(params){
  
  qData <- pData <- NULL
  for(i in 1:length(params)){
    assign(names(params[i]), params[[i]])    
  }
  
  if (prop.MV==0 && (1 != (prop.MEC + prop.POV))){
    warning("The sum of probability of POV and MEC missing values must be equal to 1.")
    return(NULL)
  }
  
  
  nbCond <- nbCond
  nRep <- nRep
  
  if (nRep==2){ # eviter les conditions a un seul replicat
    print("Only two replicats, mismatch.nRep set to FALSE")
    mismatch.nRep == FALSE
  }
  
  
  #creation dataset sans missing values
  base <- LETTERS[1:nbCond]
  FC.factor = 100
  if (!mismatch.nRep) {
    sample.names <- unlist(lapply(base, function(x)paste0(x,"_",1:nRep)))
    
  } else {
    sample.names <- unlist(lapply(base, function(x)paste0(x,"_",1:sample(2:nRep,1))))
    
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
  
  
  qData <- ll$qData
  pData <- ll$pData
  nCond <- length(unique(pData$Condition))
  interC <- intraC <- fullRandom <- 0
  # if (isTRUE(do.interC)){interC <- sample(c(TRUE, FALSE), 1)}
  # if (isTRUE(do.intraC)){intraC <- sample(c(TRUE, FALSE), 1)}
  # if (isTRUE(do.fullRandom)){fullRandom <- sample(c(TRUE, FALSE), 1)}
  if (isTRUE(do.interC)){interC <- TRUE}
  if (isTRUE(do.intraC)){intraC <- TRUE}
  if (isTRUE(do.fullRandom)){fullRandom <- TRUE}
  #### Mix columns qData ####
  
  if (fullRandom == 0) {
    
    if (interC == 1 && intraC == 0) { 
      
      print("conditions shuffled, replicates unchanged")
      interC.list <- c(1:ncol(qData))
      repl <- rep(1:nCond,table(pData$Condition)[unique(pData$Condition)])
      interC.list <- split(interC.list, repl) # general with mismatch.nRep FALSE or TRUE
      #interC.list <- split(interC.list, sort(interC.list%%nCond)) # for mismatch.nRep == FALSE
      new.interC.list <- unlist(sample(interC.list,nCond), use.names = F)
      qData <- qData[,new.interC.list]
      pData <- pData[new.interC.list,]
    }
    
    if (interC == 0 && intraC == 1) { 
      
      print("conditions unchanged, replicates shuffled")
      intraC.list <- c(1:ncol(qData))
      repl <- rep(1:nCond,table(pData$Condition)[unique(pData$Condition)])
      intraC.list <- split(intraC.list, repl)
      #intraC.list <- split(intraC.list, sort(intraC.list%%nCond))
      
      new.order <- vector()
      
      for (i in (1:nCond) ) {
        
        #print(i)
        #i=2
        new.order <- c(new.order,sample(intraC.list[[i]],length(intraC.list[[i]])))
      }
      qData <- qData[,new.order]
      pData <- pData[new.order,]
    }
    
    if (interC == 1 && intraC == 1) { 
      
      print("conditions and replicates shuffled")
      er.ra <- c(1:ncol(qData))
      repl <- rep(1:nCond,table(pData$Condition)[unique(pData$Condition)])
      er.ra <- split(er.ra, repl)
      #er.ra <- split(er.ra, sort(er.ra%%nCond))
      new.order <- vector()
      
      for (i in (1:nCond) ) {
        #print(i)
        #i=1
        new.order <- c(new.order,sample(er.ra[[i]]))
        
      }
      qData <- qData[,new.order]
      pData <- pData[new.order,]
      new.interC.list <- unlist(sample(er.ra,nCond), use.names = F)
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


##' This function generates datasets from DAPARdata with extra conditions and replicates
##'
##' @title xxxxxx
##' @param ll
##' @param nCond xxxx 
##' @param xxxxx
##' @param type xxxxx
##' @return A list of two items : qData and pData
##' @author Enora Fremy
##' @examples
##' ll <- realData(data = Exp1_R25_pept,plusCdt = 0,plusRep = 1)
realData <- function(data,plusCdt,plusRep){
  
  qData <- exprs(data)
  pData <- pData(data)
  nbCond.i <- length(unique(pData$Condition))
  nRep.i <- nrow(pData) / nbCond.i
  
  nbCond.f <- nbCond.i + plusCdt
  nRep.f <- nRep.i + plusRep
  
  if (plusCdt>0) {
    print(paste0(plusCdt," condition(s) added..."))
    # ajouter n condition
    
    # a qData
    plusCol <- matrix(nrow = nrow(qData), ncol = nRep.i*plusCdt)
    cdtDuplicated <- split(c(1:ncol(qData)), rep(1:nbCond.i,table(pData$Condition)))
    for (i in 1:plusCdt) {
      
      seq <- rep(1:plusCdt, each = nRep.i)
      coord_plusCol <- split(c(1:(ncol(plusCol))), seq )
      plusCol[,coord_plusCol[[i]]] <- jitter(qData[,unlist(sample(cdtDuplicated,1 ))],2)
    }
    qData <- cbind(qData,plusCol)
    
    # a pData
    plusRow <- as.data.frame(matrix(nrow = nRep.i*plusCdt, ncol = ncol(pData)))
    colnames(plusRow) <- colnames(pData)
    plusRow$Condition <- paste0("condition",rep(c(1:plusCdt), each = nRep.i))
    pData <- rbind(pData,plusRow)
    
  }
  
  
  base <- LETTERS[1:nbCond.f]
  sample.names <- unlist(lapply(base, function(x)paste0(x,"_",1:nRep.i)))
  colnames(qData) <- sample.names
  
  pData$Sample.name <- sample.names
  pData$Bio.Rep <- 1:nrow(pData)

  # ajouter n replicats (au hasard)
  if (plusRep>0) {
    print(paste0(plusRep," replicate(s) added..."))
    
    # 1- selectionner autant de col a dupliquer que de replicas voulus
    plusCol <- as.data.frame(matrix(NA, nrow=nrow(qData), ncol = plusRep))
    
    repSup <- sort(sample(1:ncol(qData),plusRep, replace = T) )
    print(paste0("repSup: ",repSup))
    
    for (i in 1:length(repSup)) {
      plusCol[,i] <- jitter(qData[,repSup[i]],2)
      names(plusCol)[i] <- repSup[i]
    }
    
    # 2 - placer les colonnes
    cond <- split(c(1:ncol(qData)), rep(1:nbCond.f,table(pData$Condition)))
    coord_cdt <- sapply(cond, tail, 1)
    compt = 0
    
    
    for (i in 1:length(cond)) {
      #i=2
      combien <- sum( (repSup %in% cond[[i]]) , na.rm = TRUE) 
      
      if (combien>0) {
        qData <- cbind(qData, plusCol[,1:combien, drop = F])
        
        if (i != length(cond)) {
          plusCol <- plusCol[,-c(1:combien), drop = F]
          qData <- qData[, c( 1:(tail(cond[[i]],1)+compt), (ncol(qData)- combien + 1) : ncol(qData)  , (tail(cond[[i]],1)+1+compt) : (ncol(qData)- combien))]
        }
      }
      compt = compt + combien 
    }
   
  }
  
  
  indices <- which(!grepl("^[A-Z]..",names(qData)))
  for (i in indices) {
    #i=4
    toCorrect <- unlist(strsplit(names(qData)[i-1], "_"))
    names(qData)[i] <- paste0(toCorrect[1],"_",as.numeric(toCorrect[2])+1)
  }
  
                                   
  plusRow <- as.data.frame(matrix(nrow = ncol(qData)-nrow(pData), ncol = ncol(pData)))
  colnames(plusRow) <- colnames(pData)
  plusRow$Sample.name <- setdiff(colnames(qData), pData$Sample.name )
  pData <- rbind(pData,plusRow)
  pData <- pData[order(pData$Sample.name),]
  
  indices <- which(is.na(pData$Condition))
  for (i in indices) {
    pData$Condition[i] <- pData$Condition[i-1]
  }
  
  pData$Bio.Rep <- 1:nrow(pData)
  rownames(pData) <- colnames(qData)
  
  ll <- list(qData=qData, pData=pData)
  ll <- CreateMinimalistMSnset(ll)
  return ( ll )
  
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
    ll <- c("wrapper.dapar.impute.mi"#,
            #"wrapper.impute.mle.toto"#,
            #"wrapper.impute.slsa"#,
            #"wrapper.impute.detQuant",
            #"wrapper.impute.pa"#,
            #"wrapper.impute.fixedValue",
            #"wrapper.impute.KNN"
            
    ) }
  else {
    ll <- list(list(obj,nb.iter = 1, progress.bar = FALSE)#,
               #list(obj)#,
               #list(obj)#,
               #list(obj,qval=0.025, factor=1),
               #list(obj,q.min = 0.025)#,
               #list(obj,fixVal=0),
               #list(obj,K=10)
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
##' genDatasetArgs <- list(nbCond=3, nRep=3, mismatch.nRep=FALSE, prop.MV = 0.2, size = 1000)
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
        message(w)
        },
      error = function(e) {
        message(e)
        },
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
##' genDatasetArgs <- list(nbCond=3, nRep=3,mismatch.nRep=TRUE, prop.MV = 0.2,  prop.MEC=0.3, prop.POV=0.7,size = 1000)
##' mixDatasetArgs <- list(do.interC=TRUE, do.intraC=TRUE, do.fullRandom=TRUE)
##' res <- test_imputation(nTest = 20,genDatasetArgs , mixDatasetArgs)
test_imputation <- function(nTest = 5, genDatasetArgs, mixDatasetArgs, realData, data) {
  res <- list()
  
    for (i in 1:nTest) {
      
      print(paste0("#######--------------------------- TEST ", i, " ---------------------------------#######"))
      
      if (realData) {
        ll.original <- list(qData=exprs(data),pData=pData(data))
      }
      
      
      else {
        # 1 - generation qData et pData 
        ll.original <- do.call("GenerateRandomDataset", list(genDatasetArgs))
        print(unlist(genDatasetArgs))
      }  
        
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


showDatasets <- function(res, tour,  size){
  
  data <- res[[paste0('tour', tour)]]
  for (i in 1:length(names(data))){
    print(paste0('---------',names(data)[i] ,'---------' ))
    print(exprs(data[[i]])[size,])
  }
  
}


#------------------------------------------------------------
# Tests
#------------------------------------------------------------
data("Exp1_R25_pept")

ll <- realData(Exp1_R25_pept,plusCdt = 2, plusRep = 3)
pData(ll)
# line 13 of condition A1 <- 100 to see if imputation accurated
exprs(ll)[14,1] <- 50

#genDatasetArgs <- list(nbCond=3, nRep=3, mismatch.nRep=TRUE, prop.MV = 0.4,  prop.MEC=0.3, prop.POV=0.7,size = 5000)
mixDatasetArgs <- list(do.interC=TRUE, do.intraC=TRUE, do.fullRandom=FALSE)
start_time <- Sys.time()
res.1 <- test_imputation(nTest=1,genDatasetArgs , mixDatasetArgs, realData = T, data = ll)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
mixDatasetArgs <- list(do.interC=TRUE, do.intraC=TRUE, do.fullRandom=TRUE)
res.2 <- test_imputation(nTest=1,genDatasetArgs , mixDatasetArgs, realData = T, data = ll)
end_time <- Sys.time()
end_time - start_time

# testS
#summary(res.1)
showDatasets(res.1,1,5:16)
showDatasets(res.2,2,20:23)


#save(res,file = "res_mi_test2.RData")

# Imputation complete ?
# method = c("wrapper.dapar.impute.mi","wrapper.impute.mle","wrapper.impute.slsa","wrapper.impute.detQuant","wrapper.impute.pa",
#            "wrapper.impute.fixedValue","wrapper.impute.KNN")
method="wrapper.dapar.impute.mi"
tour=1
res <- res.1
#res <- res.2
data <- res[[paste0("tour",tour)]][[paste0("obj.ori.",method)]]
(exprs(data))[!complete.cases(exprs(data)),]
data <- res[[paste0("tour",tour)]][[paste0("obj.mix.",method)]]
(exprs(data))[!complete.cases(exprs(data)),]
data <- res$tour1$obj.original
(exprs(data))[!complete.cases(exprs(data)),]