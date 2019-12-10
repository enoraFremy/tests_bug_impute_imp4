

GetListFuncs_HypothesisTest <- function(obj=NULL){
  
  ll <- NULL
  if (is.null(obj)) {
    ##liste des fonctions a tester
    ll <- c("limmaCompleteTest",
            "wrapper.t_test_Complete"
    ) }
  else {
    ll <- list(list(exprs(obj), pData(obj)),
               list(obj)
    )
  }
  return(ll)
}



##Teste les fonctions d'imputation sur 2 datasets
test_HypothesisTest_functions <- function(obj.original, obj.mixed){
  
  FUN <- GetListFuncs_HypothesisTest() # Liste des fonctions a tester
  df_list <- list()
  df_list$obj.original <- obj.original
  df_list$obj.mixed <- obj.mixed
  ll_params_original <- GetListFuncs_HypothesisTest(obj.original) # param a charger pour la fonction courante
  ll_params_mixed <- GetListFuncs_HypothesisTest(obj.mixed)
  
  # Pour chaque fonction d'imputation a tester
  for (i in 1:length(FUN)){
    #i=1
    cat("\nTest de la fonction : ",FUN[i])
    #recupere les fonctions a tester et leurs parametres
    
    
    tryCatch(
      {# execution de la fonction d'imputation a tester
        obj.original.hypo <- obj.original
        tmp.original.hypo <- do.call(FUN[i],ll_params_original[[i]])
        if (class(qData.obj.original.imputed)[1] != 'MSnSet'){
          Biobase::exprs(obj.original.imputed) <- tmp.original.imputed
        } else {
          obj.original.imputed <- tmp.original.imputed
        }
        
        obj.mixed.imputed <- obj.mixed
        tmp.mixed.imputed <- do.call(FUN[i],ll_params_mixed[[i]])
        if (class(qData.obj.mixed.imputed)[1] != 'MSnSet'){
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



test_HypothesisTest <- function(maxNbCond=3, MaxnRep=3, fixedDesign = FALSE, size = 1000, mismatch.nRep = FALSE, interC, intraC, fullRandom, nTest = 5) {
  
  res <- list()
  
  
  for (i in 1:nTest) {
    
    print(paste0("#######--------------------------- TEST ", i, " ---------------------------------#######"))
    if (isTRUE(fixedDesign)){
      nbCond <- maxNbCond
      nRep <- MaxnRep
    } else {
      nbCond = sample(c(2:maxNbCond),1)
      nRep = sample(c(1:MaxnRep),1)
    }
    
    interC = sample(c(0,1),1)
    intraC = sample(c(0,1),1)
    fullRandom = sample(c(0,1),1)
    
    cat("\n Caracteristiques du dataset :\n
        *** nCond: ",nbCond, ", nRep: ", nRep, ", interC: ", interC, ", intraC: ",intraC, ", fullRandom: ",fullRandom, " ***\n")
    
    # 1 - generation qData et pData
    res.original <- GenerateRandomDataset(nbCond=nbCond, nRep=nRep,mismatch.nRep, prop.MV = 0, size = size)
    
    # 2 - etape de melange
    res.mixed <- mix_dataset_Enora(res.original, nbCond, nRep, mismatch.nRep = FALSE, interC, intraC, fullRandom)
    
    # 3 - mise sous forme de MSnset
    obj.original <- CreateMinimalistMSnset(res.original)
    obj.mixed <- CreateMinimalistMSnset(res.mixed)
    
    # 4 -test en serie des fonctions d'imputation
    impute <- test_HypothesisTest_functions(obj.original, obj.mixed)
    tour <- paste0("tour",i)
    res[[tour]] <- impute
    
  }
  return(res)
}
