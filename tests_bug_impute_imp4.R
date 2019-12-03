library(DAPAR)
library(DAPARdata)
library(imp4p)

#### data loading ####


data("Exp1_R25_prot")
data("Exp1_R25_pept")
dat <- mvFilter(Exp1_R25_prot, type="allCond", th = 1)
#dat <- mvFilter(Exp1_R25_pept, type="allCond", th = 1)
qData <- Biobase::exprs(dat)

qData <- cbind(qData, qData[,c(1,3,5)]) # Third condition
colnames(qData) <- c("A1","A2","A3","B1","B2","B3","C1","C2","C3")

qData.NA <- is.na(qData) # keep NA for imputation

qData[,1:3] <- apply(qData[,1:3], 2, function(x) rnorm(x, mean = 1, sd = 1))
qData[,4:6] <- apply(qData[,4:6], 2, function(x) rnorm(x, mean = 500, sd = 1))
qData[,7:9] <- apply(qData[,7:9], 2, function(x) rnorm(x, mean = 1000, sd = 1))

qData[qData.NA] <- NA
rm(qData.NA)

df.init <- list(qData) # save NA tab in df.init



#### imputation ####


res.df <- list(slsa="",mi="")

# wrapper.impute.slsa
conds <- as.factor(c("A","A","A","B","B","B","C","C","C"))
res.slsa <- impute.slsa(qData, conditions=conds, nknn=15, selec="all", weight=1, ind.comp=1)

res.df[[1]] <- res.slsa
colnames(res.df$slsa) <- paste("slsa.protein",colnames(res.df$slsa), sep = "_")

# wrapper.dapar.impute.mi
conditions <- as.factor(c("A","A","A","B","B","B","C","C","C"))
repbio <- as.factor(c("1","2","3","4","5","6","7","8","9"))
reptech <- as.factor(NULL)

tab <- qData

dat.slsa = impute.rand(tab = tab, conditions = conditions)
res = estim.mix(tab = tab, tab.imp = dat.slsa, conditions = conditions, 
                x.step.mod = 300, 
                x.step.pi = 300, nb.rei = 100)
# Error in optim(init, fr, gr = NULL, x = (abs[, j] - xmin), pi_est = pi_init,  : 
# L-BFGS-B nécessite des valeurs finies de 'fn'
born = estim.bound(tab = tab, conditions = conditions, q = 0.95)
proba = prob.mcar.tab(born$tab.upper, res)
data.mi = mi.mix(tab = tab, tab.imp = dat.slsa, prob.MCAR = proba, 
                 conditions = conditions, repbio = repbio, reptech = reptech, 
                 nb.iter = nb.iter, nknn = nknn, weight = weight, selec = selec, 
                 siz = siz, ind.comp = ind.comp, methodi = methodi, q = q, 
                 progress.bar = progress.bar)
data.final = impute.pa2(tab = data.mi, conditions = conditions, 
                        q.min = q.min, q.norm = q.norm, eps = eps, distribution = distribution)

res.df[[2]] <- data.final
colnames(res.df$mi) <- paste("impute.mi.protein",colnames(res.df$mi), sep = "_")

# res.mi <- impute.mi(qData, conditions=conds, nknn=15, selec="all", weight=1, ind.comp=1)
# Error in optim(init, fr, gr = NULL, x = (abs[, j] - xmin), pi_est = pi_init,  : 
# L-BFGS-B nécessite des valeurs finies de 'fn'