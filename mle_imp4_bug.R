#------------------------------------------------------------
# Test imputation imp4p mle
#------------------------------------------------------------

# va avec mixed_df_generation.R

wrapper.impute.mle <- function(obj){
  cond <- as.factor(Biobase::pData(obj)$Condition)
  
  res <- impute.mle(Biobase::exprs(obj), conditions=cond)
  
  Biobase::exprs(obj) <-res
  return (obj)
}

obj.original.imputed <- wrapper.impute.mle(obj.original)
obj.original.mixed <- wrapper.impute.mle(obj.mixed)
original.order <- colnames(exprs(obj.original.imputed))
qData.original.mixed <- exprs(obj.original.mixed)
qData.original.mixed <- (exprs(obj.original.mixed))[,original.order]
# exprs(obj.original.mixed) ne se modifie pas

head(exprs(obj.original))
head(exprs(obj.mixed))
head(exprs(obj.original.imputed))
#head(exprs(obj.original.mixed))
head(qData.original.mixed)
expect_equal(exprs(obj.original.imputed), qData.original.mixed,tolerance=1)
# NA pas imputes

#------------------------------------------------------------
# Test imputation imp4p mle, dernier code Quentin
#------------------------------------------------------------

impute.mle=function (tab, conditions) {
  
  tab_imp=as.matrix(tab);
  
  conditions=factor(as.character(conditions),levels=as.character(unique(conditions)));
  
  nb_cond=length(levels(conditions));
  nb_rep=rep(0,nb_cond);
  k=1;
  
  for (n in 1:nb_cond){
    nb_rep[n]=sum((conditions==levels(conditions)[n]));
    xincomplete=as.matrix(tab[,(k:(k+nb_rep[n]-1))]);
    nbna=fast_apply_nb_na(xincomplete,1);
    if (sum(nbna)>0){
      xincomplete1=xincomplete[which(nbna!=nb_rep[n]),];
      nbna2=fast_apply_nb_na(xincomplete1,1);
      if (sum(nbna2)>0){
        s <- prelim.norm(xincomplete1);
        thetahat <- em.norm(s, showits = FALSE);
        rngseed(1234567);
        xcomplete1 <- imp.norm(s, thetahat, xincomplete1);
        tab_imp[which(nbna!=nb_rep[n]),(k:(k+nb_rep[n]-1))]=xcomplete1;
      }
    }
    k=k+nb_rep[n];
  }
  
  return(tab_imp)
}