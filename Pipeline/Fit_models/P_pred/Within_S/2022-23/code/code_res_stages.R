rm(list=ls())
library(BGLR)
getwd()
season <- "2022-23"
base <- "../Desktop/PhD data/P3/Pipeline/Fit_models/P_pred/Within_S"
setwd(file.path(base,season,"output"))

load("Predicted_phenomics_and_Phat_matrices.RData")
set.seed(2526); nIter<-12000; burnIn<-2000

Pheno_std <- read.csv("../../../../../3.Phenomic_profiling/output/Final_outs1/Pheno_std_matched.csv")
Pheno <- Pheno_std[Pheno_std$Season==season,]
Pheno$id <- paste(Pheno$Season,Pheno$Genotype,sep="_")
Pheno <- Pheno[match(id,Pheno$id),]

traits <- intersect(colnames(pred_X),names(Pheno))
Xres <- as.matrix(Pheno[,traits,drop=FALSE])-pred_X[,traits,drop=FALSE]
storage.mode(Xres)<-"numeric"

Pres_S <- pred_R_S <- pred_GR_S <- vector("list",12)

for(s in 1:12){
  
  cols <- grep(paste0("^S",s,"_"),colnames(Xres))
  Z <- Xres[,cols,drop=FALSE]
  
  mu <- colMeans(Z[train,,drop=FALSE],na.rm=TRUE)
  sdv <- apply(Z[train,,drop=FALSE],2,sd,na.rm=TRUE)
  keep <- is.finite(mu)&is.finite(sdv)&sdv>0
  
  Z <- scale(Z[,keep,drop=FALSE],mu[keep],sdv[keep])
  Z[!is.finite(Z)] <- 0
  
  Pres_S[[s]] <- tcrossprod(Z)/ncol(Z)
  Pres_S[[s]] <- Pres_S[[s]]/mean(diag(Pres_S[[s]]))
  rownames(Pres_S[[s]])<-colnames(Pres_S[[s]])<-id
  
  pred_R_S[[s]] <- pred_GR_S[[s]] <- Y*NA_real_
  
  for(j in seq_len(ncol(Y))){
    y<-Y[,j]; y[test]<-NA
    
    f1<-BGLR(y,ETA=list(R=list(K=Pres_S[[s]],model="RKHS")),
             nIter=nIter,burnIn=burnIn,verbose=FALSE)
    
    f2<-BGLR(y,ETA=list(G=list(K=G,model="RKHS"),
                        R=list(K=Pres_S[[s]],model="RKHS")),
             nIter=nIter,burnIn=burnIn,verbose=FALSE)
    
    pred_R_S[[s]][,j] <- f1$yHat
    pred_GR_S[[s]][,j] <- f2$yHat
  }
}

accuracy <- do.call(rbind,lapply(1:12,function(s)
  do.call(rbind,lapply(seq_len(ncol(Y)),function(j)
    do.call(rbind,lapply(c("Train","Test"),function(z){
      ii<-which(Set==z)
      data.frame(Stage=paste0("S",s),Trait=colnames(Y)[j],Set=z,
                 Model=c("Residual","G+Residual"),
                 Correlation=c(
                   cor(Y[ii,j],pred_R_S[[s]][ii,j],use="complete.obs"),
                   cor(Y[ii,j],pred_GR_S[[s]][ii,j],use="complete.obs")),
                 RMSE=c(
                   sqrt(mean((Y[ii,j]-pred_R_S[[s]][ii,j])^2,na.rm=TRUE)),
                   sqrt(mean((Y[ii,j]-pred_GR_S[[s]][ii,j])^2,na.rm=TRUE))))
    }))
  ))
))

Predictions <- do.call(rbind,lapply(1:12,function(s)
  data.frame(id=rep(id,ncol(Y)),Set=rep(Set,ncol(Y)),Stage=paste0("S",s),
             Trait=rep(colnames(Y),each=length(id)),
             Observed=as.vector(Y),
             Predicted_Residual=as.vector(pred_R_S[[s]]),
             Predicted_GResidual=as.vector(pred_GR_S[[s]]))))

summary_stage <- aggregate(Correlation~Stage+Model,
                           subset(accuracy,Set=="Test"),mean,na.rm=TRUE)

write.csv(accuracy,file.path("Res","Stage_residual_accuracy.csv"),row.names=FALSE)
write.csv(Predictions,file.path("Res","Stage_residual_predictions.csv"),row.names=FALSE)
write.csv(summary_stage,file.path("Res","Stage_residual_summary.csv"),row.names=FALSE)

save(Pres_S,pred_R_S,pred_GR_S,
     file=file.path("Res","Stage_residual_models.RData"))

print(summary_stage[order(-summary_stage$Correlation),])

