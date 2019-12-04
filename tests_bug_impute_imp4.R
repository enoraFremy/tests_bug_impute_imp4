library(DAPAR)
library(DAPARdata)
library(imp4p)
library(testthat)

#### Function to test the egality between the original (ordered) dataset and
#### the mixed one which has been imputed
testSpecialDatasets <- function(qData.original, qData.mixed){
  
  original.order <- colnames(qData.original)
  qData.mixed.reordered <- qData.mixed[,original.order]
  
  testthat::expect_equal(qData.mixed.reordered, qData.original,tolerance=1)
}


impute.mi.test <- function(qData, pData, nb.iter = 3, 
                               nknn = 15, selec = 600, siz = 500, weight = 1, ind.comp = 1, 
                               progress.bar = TRUE, x.step.mod = 300, 
                               x.step.pi = 300, nb.rei = 100, method = 4, gridsize = 300, 
                               q = 0.95, q.min = 0, q.norm = 3, eps = 0, methodi = "slsa",
                               lapala = TRUE,
                               distribution="unif"){
  
  ## order exp and pData table before using imp4p functions
  conds <- factor(pData$Condition, levels=unique(pData$Condition))
  sample.names.old <- pData$Sample.name
  sTab <- pData
  new.order <- unlist(lapply(split(sTab, conds), function(x) {x['Sample.name']}))
  qData <- qData[,new.order]
  sTab <- pData[new.order,]
  
  
  
  conditions <- as.factor(sTab$Condition)
  repbio <- as.factor(sTab$Bio.Rep)
  reptech <-as.factor(sTab$Tech.Rep)
  
  tab <- qData
  
  if (progress.bar == TRUE) {
    cat(paste("\n 1/ Initial imputation under the MCAR assumption with impute.rand ... \n  "))
  }
  dat.slsa = imp4p::impute.rand(tab = tab, conditions = conditions)
  
  if (progress.bar == TRUE) {
    cat(paste("\n 2/ Estimation of the mixture model in each sample... \n  "))
  }
  res = estim.mix(tab = tab, tab.imp = dat.slsa, conditions = conditions, 
                  x.step.mod = x.step.mod, 
                  x.step.pi = x.step.pi, nb.rei = nb.rei)
  
  
  if (progress.bar == TRUE) {
    cat(paste("\n 3/ Estimation of the probabilities each missing value is MCAR... \n  "))
  }
  born = estim.bound(tab = tab, conditions = conditions, q = q)
  proba = prob.mcar.tab(born$tab.upper, res)
  
  
  if (progress.bar == TRUE) {
    cat(paste("\n 4/ Multiple imputation strategy with mi.mix ... \n  "))
  }
  data.mi = mi.mix(tab = tab, tab.imp = dat.slsa, prob.MCAR = proba, 
                   conditions = conditions, repbio = repbio, reptech = reptech, 
                   nb.iter = nb.iter, nknn = nknn, weight = weight, selec = selec, 
                   siz = siz, ind.comp = ind.comp, methodi = methodi, q = q, 
                   progress.bar = progress.bar)
  
  if (lapala == TRUE){
    if (progress.bar == TRUE) {
      cat(paste("\n\n 5/ Imputation of rows with only missing values in a condition with impute.pa ... \n  "))
    }
    data.final = impute.pa2(tab = data.mi, conditions = conditions, 
                            q.min = q.min, q.norm = q.norm, eps = eps, distribution = distribution)
  } else {
    data.final <- data.mi
  }
  
  
  # restore previous order
  colnames(data.final) <- new.order
  data.final <- data.final[,sample.names.old]
  
  return(data.final)
}

#### imputation of missing values in qData, method slsa or impute.mi
imputation <- function(qData, method) { # "slsa" ou "impute.mi"
  
  if (method == "slsa") { # wrapper.impute.slsa
    
    conds <- as.factor(c("A","A","A","B","B","B","C","C","C"))
    res.slsa <- impute.slsa(qData, conditions=conds, nknn=15, selec="all", weight=1, ind.comp=1)
    
    return(res.slsa)
  }
  
  else { # wrapper.dapar.impute.mi
    # gere pas plus de deux cdts ?
    tab <- abs(qData[,-c(7:9)])
    head(tab)
    conditions <- factor(c("A","A","A","B","B","B"))
    repbio <- factor(c("1","2","3","4","5","6"))
    reptech <- as.factor(NULL)
    
    # marche pas meme avec deux conditions -> pb avec les rnorm
    qData <- Biobase::exprs(Exp1_R25_prot)
    tab <- qData
    # la, estim.mix fonctionne
    
    # essai de rajouter une troisieme condition a qData qui fonctionne
    tab <- cbind(qData, qData[,c(1,3,5)])
    head(tab)
    conditions <- factor(c("A","A","A","B","B","B","C","C","C"))
    repbio <- factor(c("1","2","3","4","5","6","7","8","9"))
    reptech <- as.factor(NULL)
    # Fonctionne aussi -> pb initial vient de rnorm ou colnames(qData)
    
    # essai donnees rnorm avec colnames differentes
    colnames(qData) <- c("Intensity_C_R1","Intensity_C_R2","Intensity_C_R3","Intensity_D_R1",
                         "Intensity_D_R2","Intensity_D_R3", "Intensity_E_R3","Intensity_E_R3",
                         "Intensity_E_R3")
    tab <- qData
    # Encore erreur optim() -> pb avec rnorm...comment ?
    
    dat.slsa = impute.rand(tab = tab, conditions = conditions)
    
    res = estim.mix(tab = tab, tab.imp = dat.slsa, conditions = conditions, 
                    x.step.mod = 300, 
                    x.step.pi = 300, nb.rei = 100)
    # 
    # Error in optim(init, fr, gr = NULL, x = (abs[, j] - xmin), pi_est = pi_init,  : 
    # L-BFGS-B n?cessite des valeurs finies de 'fn'
    born = estim.bound(tab = tab, conditions = conditions, q = 0.95)
    proba = prob.mcar.tab(born$tab.upper, res)
    data.mi = mi.mix(tab = tab, tab.imp = dat.slsa, prob.MCAR = proba, 
                     conditions = conditions, repbio = repbio, reptech = reptech, 
                     nb.iter = nb.iter, nknn = nknn, weight = weight, selec = selec, 
                     siz = siz, ind.comp = ind.comp, methodi = methodi, q = q, 
                     progress.bar = progress.bar)
    data.final = impute.pa2(tab = data.mi, conditions = conditions, 
                            q.min = q.min, q.norm = q.norm, eps = eps, distribution = distribution)
    
    qData.mixed <- data.final
    
  }
  return(qData.mixed)
  
}


#### data 1 #### [ "A1" "A2" "A3" "B1" "B2" "B3" "C1" "C2" "C3" ]

data("Exp1_R25_prot")


## Creation du dataset original
nbCond <- 3
qData <- Biobase::exprs(Exp1_R25_prot)
pData <- Biobase::pData(Exp1_R25_prot)
qData <- cbind(qData, qData[,c(1,3,5)]) # Third condition
colnames(qData) <- c("A1","A2","A3","B1","B2","B3","C1","C2","C3")
pData <- rbind(pData, pData[1:nbCond,])
pData$Condition <- c(rep("A",nbCond), rep("B", nbCond), rep("C", nbCond))
pData$Sample.name <- paste0(pData$Condition, 1:nbCond)
rownames(pData) <- pData$Sample.name

qData.NA <- is.na(qData) # keep NA for imputation

# qData[,1:3] <- apply(qData[,1:3], 2, function(x) rnorm(x, mean = 10, sd = 1))
# qData[,4:6] <- apply(qData[,4:6], 2, function(x) rnorm(x, mean = 50, sd = 1))
# qData[,7:9] <- apply(qData[,7:9], 2, function(x) rnorm(x, mean = 100, sd = 1))
qData[,1:3] <- apply(qData[,1:3], 2, function(x) sample(1:20, nrow(qData), replace = T))
qData[,4:6] <- apply(qData[,4:6], 2, function(x) sample(49:51, nrow(qData), replace = T))
qData[,7:9] <- apply(qData[,7:9], 2, function(x) sample(91:110, nrow(qData), replace = T))
qData[qData.NA] <- NA
rm(qData.NA)

qData.original <- qData
pData.original <- pData


## creation du dataset melange
#new.order <- xxx
qData.mixed <- qData
pData.mixed <- pData


## Imputation des deux datasets

qData.original <- impute.mi.test(qData.original, pData.original)
qData.mixed <- impute.mi.test(qData.mixed, pData.mixed)


## Test d'egalite des deux resultats
testSpecialDatasets(qData.original, qData.mixed)


--------------




df.test <- list(qData.original = qData.original, qData.mixed = qData.mixed)

#### data 2 #### [ "B1" "B2" "B3" "C1" "C2" "C3" "A1" "A2" "A3" ]
qData.original <- qData[,c(4:6,7:9,1:3)]
qData.mixed <- imputation(qData.original, "slsa")
testSpecialDatasets(qData.original, qData.mixed)

df.test$qData.original.desorder <- qData.original
df.test$qData.mixed.desorder <- qData.mixed


#### data 3 #### [ "A1" "B1" "C1" "A2" "B2" "C2" "A3" "B3" "C3" ]
head(qData)
qData.original <- qData[,c(1,4,7,2,5,8,3,6,9)]
qData.mixed <- imputation(qData.original, "slsa")
testSpecialDatasets(qData.original, qData.mixed)

df.test$qData.original.desorder2 <- qData.original
df.test$qData.mixed.desorder2 <- qData.mixed

boxplot(df.test$qData.original)
boxplot(df.test$qData.mixed)
boxplot(df.test$qData.original.desorder)
boxplot(df.test$qData.mixed.desorder)
boxplot(df.test$qData.original.desorder2)
boxplot(df.test$qData.mixed.desorder2)
