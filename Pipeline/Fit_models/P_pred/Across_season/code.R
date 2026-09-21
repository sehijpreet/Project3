rm(list=ls())
library(BGLR)
setwd('C:/Users/sidhu/OneDrive - University of Florida/Desktop')

ROOT <- "../Desktop/PhD data/P3/Pipeline"


Pheno <- read.csv(file.path(ROOT,"3.Phenomic_profiling/output/Final_outs1/Pheno_std_matched.csv"))
Res   <- read.csv(file.path(ROOT,"3.Phenomic_profiling/output/Final_outs1/Res1_matched.csv"))
load(file.path(ROOT,"4.G/output/ZGZt.rda"))
load(file.path(ROOT,"6.P/output/P_std.rda"))
dim(ZGZt)
dim(P_std)
ZGZt[1:5, 1:5]
P_std[1:5, 1:5]

Pheno$id <- paste(Pheno$Season,Pheno$Genotype,sep="_")
Res$id   <- paste(Res$Season,Res$Genotype,sep="_")

targets <- sort(intersect(unique(Pheno$Season),unique(Res$Season)))
nIter <- 12000; burnIn <- 2000

for(target in targets){
  
  cat("\n================",target,"================\n")
  
  outdir <- file.path(ROOT,"Fit_models/P_pred/Across_season",target,"output")
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  
  id <- sort(Reduce(intersect,list(Pheno$id,Res$id,rownames(ZGZt),rownames(P_std))))
  ph <- Pheno[match(id,Pheno$id),]
  rs <- Res[match(id,Res$id),]
  
  Y <- as.matrix(rs[,grep("_cum$",names(rs))]); rownames(Y)<-id
  G <- ZGZt[id,id]; G <- G/mean(diag(G))
  P <- P_std[id,id]; P <- P/mean(diag(P))
  
  test <- which(ph$Season==target)
  train <- which(ph$Season!=target)
  Set <- ifelse(seq_along(id)%in%test,"Test","Train")
  
  cat("Test:",target,"\nTraining:",paste(unique(ph$Season[train]),collapse=", "),"\n")
  
  # -------- Predict UAV traits using ALL OTHER seasons --------
  traits <- setdiff(names(ph),c("Season","Genotype","id"))
  traits <- traits[sapply(ph[,traits,drop=FALSE],is.numeric)]
  
  pred_X <- matrix(NA,nrow(Y),length(traits),dimnames=list(id,traits))
  
  for(j in seq_along(traits)){
    y <- ph[[traits[j]]]; y[test] <- NA
    if(sum(!is.na(y[train]))<10) next
    
    fit <- BGLR(y,ETA=list(G=list(K=G,model="RKHS")),
                nIter=nIter,burnIn=burnIn,verbose=FALSE)
    
    pred_X[,j] <- fit$yHat
  }
  
  # -------- P-hat --------
  keep <- apply(pred_X,2,function(x)
    all(is.finite(x)) && is.finite(sd(x[train])) && sd(x[train])>0)
  
  Z <- pred_X[,keep,drop=FALSE]
  Z <- scale(Z,colMeans(Z[train,,drop=FALSE]),
             apply(Z[train,,drop=FALSE],2,sd))
  
  Phat <- tcrossprod(Z)/ncol(Z)
  Phat <- Phat/mean(diag(Phat))
  rownames(Phat)<-colnames(Phat)<-id
  
  # -------- Residual = observed UAV - predicted UAV --------
  Xobs <- as.matrix(ph[,colnames(pred_X),drop=FALSE])
  storage.mode(Xobs) <- "numeric"
  Xres <- Xobs-pred_X
  
  sdR <- apply(Xres[train,,drop=FALSE],2,sd,na.rm=TRUE)
  keepR <- keep & is.finite(sdR) & sdR>0
  
  R <- Xres[,keepR,drop=FALSE]
  muR <- colMeans(R[train,,drop=FALSE],na.rm=TRUE)
  sdR <- apply(R[train,,drop=FALSE],2,sd,na.rm=TRUE)
  
  R <- scale(R,muR,sdR)
  R[!is.finite(R)] <- 0
  
  Pres <- tcrossprod(R)/ncol(R)
  Pres <- Pres/mean(diag(Pres))
  rownames(Pres)<-colnames(Pres)<-id
  
  # -------- Models --------
  mods <- list(
    G=list(G=list(K=G,model="RKHS")),
    P=list(P=list(K=P,model="RKHS")),
    `G+P`=list(G=list(K=G,model="RKHS"),P=list(K=P,model="RKHS")),
    Phat=list(PH=list(K=Phat,model="RKHS")),
    `G+Phat`=list(G=list(K=G,model="RKHS"),PH=list(K=Phat,model="RKHS")),
    Residual=list(R=list(K=Pres,model="RKHS")),
    `G+Residual`=list(G=list(K=G,model="RKHS"),R=list(K=Pres,model="RKHS"))
  )
  
  Pred <- lapply(mods,function(x) Y*NA_real_)
  
  for(m in names(mods)){
    cat("Running",target,m,"\n")
    
    for(j in seq_len(ncol(Y))){
      y <- Y[,j]; y[test] <- NA
      
      fit <- BGLR(y,ETA=mods[[m]],
                  nIter=nIter,burnIn=burnIn,verbose=FALSE)
      
      Pred[[m]][,j] <- fit$yHat
    }
  }
  
  # -------- Accuracy --------
  safe_cor <- function(y,p){
    ok <- is.finite(y)&is.finite(p)
    if(sum(ok)<3) return(NA_real_)
    cor(y[ok],p[ok])
  }
  
  acc <- do.call(rbind,lapply(names(Pred),function(m)
    do.call(rbind,lapply(c("Train","Test"),function(s){
      ii <- which(Set==s)
      
      data.frame(
        Season=target,Trait=colnames(Y),Set=s,Model=m,
        Correlation=sapply(seq_len(ncol(Y)),function(j)
          safe_cor(Y[ii,j],Pred[[m]][ii,j])),
        RMSE=sapply(seq_len(ncol(Y)),function(j)
          sqrt(mean((Y[ii,j]-Pred[[m]][ii,j])^2,na.rm=TRUE)))
      )
    }))
  ))
  
  Predictions <- do.call(rbind,lapply(names(Pred),function(m)
    data.frame(
      id=rep(id,ncol(Y)),
      Season=rep(ph$Season,ncol(Y)),
      Set=rep(Set,ncol(Y)),
      Model=m,
      Trait=rep(colnames(Y),each=length(id)),
      Observed=as.vector(Y),
      Predicted=as.vector(Pred[[m]])
    )
  ))
  
  summary <- aggregate(
    Correlation~Model,
    subset(acc,Set=="Test" & !is.na(Correlation)),
    mean
  )
  
  write.csv(data.frame(id=id,Season=ph$Season,Set=Set,pred_X),
            file.path(outdir,"Predicted_all_UAV_traits.csv"),row.names=FALSE)
  
  write.csv(data.frame(id=id,Season=ph$Season,Set=Set,Xres),
            file.path(outdir,"Residual_UAV_traits.csv"),row.names=FALSE)
  
  write.csv(acc,file.path(outdir,"Across_season_accuracy.csv"),row.names=FALSE)
  write.csv(Predictions,file.path(outdir,"Across_season_predictions.csv"),row.names=FALSE)
  write.csv(summary,file.path(outdir,"Mean_test_accuracy.csv"),row.names=FALSE)
  
  save(pred_X,Phat,Pres,G,P,Pred,Y,id,train,test,Set,
       file=file.path(outdir,"Across_season_models.RData"))
  
  print(summary[order(-summary$Correlation),])
}
