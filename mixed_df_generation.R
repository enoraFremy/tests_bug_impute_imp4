library(DAPAR)
library(DAPARdata)
library(imp4p)
library(testthat)

# "wrapper.impute.mle" fonctionne pas, NA pas remplaces

## Noms et parametres des fonctions d'imputation a tester
## utilises par do.call
GetListFuncs <- function(obj=NULL){
  
  ll <- NULL
  if (is.null(obj)) {
    ##liste des fonctions a tester
    ll <- c("wrapper.dapar.impute.mi",
            "wrapper.impute.slsa",
            "wrapper.impute.detQuant",
            "wrapper.impute.pa",
            "wrapper.impute.fixedValue",
            "wrapper.impute.KNN"#,
            #"wrapper.impute.mle"
    ) }
  else {
    ll <- list(list(obj,nb.iter = 1),
               
               list(obj),
               list(obj,qval=0.025, factor=1),
               list(obj,q.min = 0.025),
               list(obj,fixVal=0),
               list(obj,K=10)#,
               #list(obj)
    )
  }
  return(ll)
}

##Teste les fonctions d'imputation sur 2 datasets
test_impute_functions <- function(obj.original, obj.mixed){
  FUN <- GetListFuncs() # Liste des fonctions a tester
  # Pour chaque fonction d'imputation ÃÂ  tester
  for (i in 1:length(FUN)){
    #i=1
    print(paste0("test de la fonction : ",FUN[i]))
    #recupere les fonctions a tester et leurs parametres
    ll_params_original <- GetListFuncs(obj.original) # param a charger pour la fonction courante
    ll_params_mixed <- GetListFuncs(obj.mixed)
    
    # execution de la fonction d'imputation a tester
    obj.original.imputed <- do.call(FUN[i],ll_params_original[[i]])
    obj.original.mixed <- do.call(FUN[i],ll_params_mixed[[i]])
    
    # tri des colonnes du dataset melange suivant l'ordre des 
    # colonnes du dataset original
    original.order <- colnames(exprs(obj.original.imputed))
    qData.original.mixed <- exprs(obj.original.mixed)
    qData.original.mixed <- (exprs(obj.original.mixed))[,original.order]
    # exprs(obj.original.mixed) ne se modifie pas
    
    head(exprs(obj.original))
    head(exprs(obj.mixed))
    head(exprs(obj.original.imputed))
    #head(exprs(obj.original.mixed))
    head(qData.original.mixed)
    
    # test de comparaison
    expect_equal(exprs(obj.original.imputed), qData.original.mixed,tolerance=1)
    
  }
}

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
    
    for (i in (1:ncol(qData.plus))) {
      random_col_qData <- sample(ncol(qData),1)
      #print(paste0("Column qData random: ", random_col_qData))
      qData.plus[,i] <- qData[,random_col_qData]
    }
    
    colnames(qData.plus) <- pData$Sample.name
    qData <- as.data.frame(qData.plus)
    
    res <- list(pData=pData, qData.original=qData)
    
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
  
  res$qData.mixed <- qData
  return(res)
  
}


# #------------------------------------------------------------
# data("Exp1_R25_pept")
# qData <- (Biobase::exprs(Exp1_R25_pept))[1:1000,]
# pData <- Biobase::pData(Exp1_R25_pept)
# 
# res <- df_generation(qData, pData, nCond = 3, nRep = 3, mismatch.nRep = FALSE, interC = 0, intraC = 0, fullRandom = 1)
# #View(res$pData)
# 
# head(res$qData.original)
# head(res$qData.mixed)
# res$pData
# 
# obj.original <- CreateMinimalistMSnset(res$qData.original, res$pData)
# obj.mixed <- CreateMinimalistMSnset(res$qData.mixed, res$pData[colnames(res$qData.mixed),])
# 
# test_impute_functions(obj.original, obj.mixed)

#------------------------------------------------------------
# Automatisation
#------------------------------------------------------------
data("Exp1_R25_pept")
qData <- (Biobase::exprs(Exp1_R25_pept))[1:1000,]
pData <- Biobase::pData(Exp1_R25_pept)

test_imputation <- function(qData, pData, nCond, nRep, mismatch.nRep = FALSE, interC, intraC, fullRandom) {
  
  cat("\n *** nCond: ",nCond, ", nRep: ", nRep, ", interC: ", interC, ", intraC: ",intraC, ", fullRandom: ",fullRandom, " ***\n")
  
  res <- df_generation(qData, pData, nCond, nRep, mismatch.nRep = FALSE, interC, intraC, fullRandom)
  
  head(res$qData.original)
  head(res$qData.mixed)
  res$pData
  
  obj.original <- CreateMinimalistMSnset(res$qData.original, res$pData)
  obj.mixed <- CreateMinimalistMSnset(res$qData.mixed, res$pData[colnames(res$qData.mixed),])
  
  test_impute_functions(obj.original, obj.mixed)
  
  
}

# nCond = 3
nRep = 2
# interC = 1
# intraC = 0
# fullRandom = 1

for (i in 1:5) {
  nCond = sample(c(2:5),1)
  #nRep = sample(c(2:4),1)
  interC = sample(c(0,1),1)
  intraC = sample(c(0,1),1)
  fullRandom = sample(c(0,1),1)
  
  test_imputation(qData, pData, nCond, nRep, mismatch.nRep = FALSE, interC, intraC, fullRandom)
  
}
#------------------------------------------------------------
# impute.mi.test <- function(qData, pData, nb.iter = 3, 
#                            nknn = 15, selec = 600, siz = 500, weight = 1, ind.comp = 1, 
#                            progress.bar = TRUE, x.step.mod = 300, 
#                            x.step.pi = 300, nb.rei = 100, method = 4, gridsize = 300, 
#                            q = 0.95, q.min = 0, q.norm = 3, eps = 0, methodi = "slsa",
#                            lapala = TRUE,
#                            distribution="unif"){
#   
#   ## order exp and pData table before using imp4p functions
#   conds <- factor(pData$Condition, levels=unique(pData$Condition))
#   sample.names.old <- pData$Sample.name
#   sTab <- pData
#   new.order <- unlist(lapply(split(sTab, conds), function(x) {x['Sample.name']}))
#   qData <- qData[,new.order]
#   sTab <- pData[new.order,]
#   
#   
#   
#   conditions <- as.factor(sTab$Condition)
#   repbio <- as.factor(sTab$Bio.Rep)
#   reptech <-as.factor(sTab$Tech.Rep)
#   
#   tab <- qData
#   
#   if (progress.bar == TRUE) {
#     cat(paste("\n 1/ Initial imputation under the MCAR assumption with impute.rand ... \n  "))
#   }
#   dat.slsa = imp4p::impute.rand(tab = tab, conditions = conditions)
#   
#   if (progress.bar == TRUE) {
#     cat(paste("\n 2/ Estimation of the mixture model in each sample... \n  "))
#   }
#   res = estim.mix(tab = tab, tab.imp = dat.slsa, conditions = conditions, 
#                   x.step.mod = x.step.mod, 
#                   x.step.pi = x.step.pi, nb.rei = nb.rei)
#   
#   
#   if (progress.bar == TRUE) {
#     cat(paste("\n 3/ Estimation of the probabilities each missing value is MCAR... \n  "))
#   }
#   born = estim.bound(tab = tab, conditions = conditions, q = q)
#   proba = prob.mcar.tab(born$tab.upper, res)
#   
#   
#   if (progress.bar == TRUE) {
#     cat(paste("\n 4/ Multiple imputation strategy with mi.mix ... \n  "))
#   }
#   data.mi = mi.mix(tab = tab, tab.imp = dat.slsa, prob.MCAR = proba, 
#                    conditions = conditions, repbio = repbio, reptech = reptech, 
#                    nb.iter = nb.iter, nknn = nknn, weight = weight, selec = selec, 
#                    siz = siz, ind.comp = ind.comp, methodi = methodi, q = q, 
#                    progress.bar = progress.bar)
#   
#   if (lapala == TRUE){
#     if (progress.bar == TRUE) {
#       cat(paste("\n\n 5/ Imputation of rows with only missing values in a condition with impute.pa ... \n  "))
#     }
#     data.final = impute.pa2(tab = data.mi, conditions = conditions, 
#                             q.min = q.min, q.norm = q.norm, eps = eps, distribution = distribution)
#   } else {
#     data.final <- data.mi
#   }
#   
#   
#   # restore previous order
#   colnames(data.final) <- new.order
#   data.final <- data.final[,sample.names.old]
#   
#   return(data.final)
# }
# qData.original <- impute.mi.test(res$qData.original, res$pData)
# qData.mixed <- impute.mi.test(res$qData.mixed, res$pData)
# head(qData.original)
# head(qData.mixed)
# testSpecialDatasets(qData.original, qData.mixed)
# testSpecialDatasets(qData.original, qData.mixed)
