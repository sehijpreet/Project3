rm(list=ls())
library(BGLR)
getwd()
season <- "2022-23"
base <- "../Desktop/PhD data/P3/Pipeline/Fit_models/P_pred/Within_S"
setwd(file.path(base,season,"output"))

dir.create("Res",showWarnings=FALSE)
load("Predicted_phenomics_and_Phat_matrices.RData")

set.seed(2526); nIter<-12000; burnIn<-2000

Pheno_std <- read.csv("../../../../../3.Phenomic_profiling/output/Final_outs1/Pheno_std_matched.csv")
Pheno <- Pheno_std[Pheno_std$Season==season,]
Pheno$id <- paste(Pheno$Season,Pheno$Genotype,sep="_")
Pheno <- Pheno[match(id,Pheno$id),]

traits <- intersect(colnames(pred_X),names(Pheno))
Xobs <- as.matrix(Pheno[,traits,drop=FALSE])
Xpred <- pred_X[,traits,drop=FALSE]
storage.mode(Xobs)<-"numeric"

Xres <- Xobs-Xpred
sdv <- apply(Xres[train,,drop=FALSE],2,sd,na.rm=TRUE)
keep <- is.finite(sdv)&sdv>0

Z <- Xres[,keep,drop=FALSE]
mu <- colMeans(Z[train,,drop=FALSE],na.rm=TRUE)
sdv <- apply(Z[train,,drop=FALSE],2,sd,na.rm=TRUE)
Z <- scale(Z,mu,sdv)
Z[!is.finite(Z)] <- 0

Pres <- tcrossprod(Z)/ncol(Z)
rownames(Pres)<-colnames(Pres)<-id
Pres <- Pres/mean(diag(Pres))

mods <- list(
  Residual=list(R=list(K=Pres,model="RKHS")),
  `G+Residual`=list(G=list(K=G,model="RKHS"),R=list(K=Pres,model="RKHS"))
)

Pred <- lapply(mods,function(x) Y*NA_real_)

for(m in names(mods))
  for(j in seq_len(ncol(Y))){
    y<-Y[,j]; y[test]<-NA
    fit<-BGLR(y,ETA=mods[[m]],nIter=nIter,burnIn=burnIn,verbose=FALSE)
    Pred[[m]][,j]<-fit$yHat
  }

acc <- do.call(rbind,lapply(names(Pred),function(m)
  do.call(rbind,lapply(c("Train","Test"),function(s){
    ii<-which(Set==s)
    data.frame(Season=season,Trait=colnames(Y),Set=s,Model=m,
               Correlation=sapply(seq_len(ncol(Y)),function(j) cor(Y[ii,j],Pred[[m]][ii,j],use="complete.obs")),
               RMSE=sapply(seq_len(ncol(Y)),function(j) sqrt(mean((Y[ii,j]-Pred[[m]][ii,j])^2,na.rm=TRUE))))
  }))
))

Predictions <- do.call(rbind,lapply(names(Pred),function(m)
  data.frame(id=rep(id,ncol(Y)),Set=rep(Set,ncol(Y)),Model=m,
             Trait=rep(colnames(Y),each=length(id)),
             Observed=as.vector(Y),Predicted=as.vector(Pred[[m]]))))

write.csv(data.frame(id=id,Set=Set,Xres),file.path("Res","Residual_UAV_traits.csv"),row.names=FALSE)
write.csv(acc,file.path("Res","Residual_models_accuracy.csv"),row.names=FALSE)
write.csv(Predictions,file.path("Res","Residual_models_predictions.csv"),row.names=FALSE)

save(Pres,Pred,file=file.path("Res","Residual_models.RData"))
print(aggregate(Correlation~Model,subset(acc,Set=="Test"),mean,na.rm=TRUE))

