rm(list=ls())
setwd("../Desktop/PhD data/P3/Pipeline/Fit_models/P_pred/Within_S/2025-26/output")
getwd()

dir.create("Res",showWarnings=FALSE); setwd("Res")
library(BGLR)

set.seed(2526); nIter<-12000; burnIn<-2000

Pheno_std<-read.csv("../../../../../../3.Phenomic_profiling/output/Final_outs1/Pheno_std_matched.csv")
Res<-read.csv("../../../../../../3.Phenomic_profiling/output/Final_outs1/Res1_matched.csv")
load("../../../../../../4.G/output/ZGZt.rda")
load("../Predicted_phenomics_and_Phat_matrices.RData")

Pheno<-Pheno_std[Pheno_std$Season=="2025-26",]
Res25<-Res[Res$Season=="2025-26",]
Pheno$id<-paste(Pheno$Season,Pheno$Genotype,sep="_")
Res25$id<-paste(Res25$Season,Res25$Genotype,sep="_")

id<-sort(Reduce(intersect,list(Pheno$id,Res25$id,rownames(ZGZt))))
Pheno<-Pheno[match(id,Pheno$id),]; Res25<-Res25[match(id,Res25$id),]

Y<-as.matrix(Res25[,grep("_cum$",names(Res25))]); rownames(Y)<-id
G<-ZGZt[id,id]; G<-G/mean(diag(G))
stopifnot(nrow(pred_X)==length(id)); rownames(pred_X)<-id

set.seed(2526)
test<-sample(seq_len(nrow(Y)),round(.30*nrow(Y)))
train<-setdiff(seq_len(nrow(Y)),test)
Set<-ifelse(seq_len(nrow(Y))%in%test,"Test","Train")

traits<-intersect(colnames(pred_X),names(Pheno))
Xres<-as.matrix(Pheno[,traits])-as.matrix(pred_X[,traits])
storage.mode(Xres)<-"numeric"; rownames(Xres)<-id

Pres_S<-pred_Res_S<-pred_GRes_S<-vector("list",12)

for(s in 1:12){
  cols<-grep(paste0("^S",s,"_"),colnames(Xres))
  Z<-Xres[,cols,drop=FALSE]
  mu<-colMeans(Z[train,,drop=FALSE],na.rm=TRUE)
  sdv<-apply(Z[train,,drop=FALSE],2,sd,na.rm=TRUE)
  keep<-is.finite(mu)&is.finite(sdv)&sdv>0
  Z<-scale(Z[,keep,drop=FALSE],mu[keep],sdv[keep])
  Z[!is.finite(Z)]<-0
  
  Pres_S[[s]]<-tcrossprod(Z)/ncol(Z)
  Pres_S[[s]]<-Pres_S[[s]]/mean(diag(Pres_S[[s]]))
  rownames(Pres_S[[s]])<-colnames(Pres_S[[s]])<-id
  
  pred_Res_S[[s]]<-pred_GRes_S[[s]]<-Y*NA_real_
  
  for(j in seq_len(ncol(Y))){
    y<-Y[,j]; y[test]<-NA
    
    fit1<-BGLR(y,ETA=list(list(K=Pres_S[[s]],model="RKHS")),
               nIter=nIter,burnIn=burnIn,verbose=FALSE)
    
    fit2<-BGLR(y,ETA=list(list(K=G,model="RKHS"),
                          list(K=Pres_S[[s]],model="RKHS")),
               nIter=nIter,burnIn=burnIn,verbose=FALSE)
    
    pred_Res_S[[s]][,j]<-fit1$yHat
    pred_GRes_S[[s]][,j]<-fit2$yHat
  }
}

accuracy_stage<-do.call(rbind,lapply(1:12,function(s)
  do.call(rbind,lapply(seq_len(ncol(Y)),function(j)
    do.call(rbind,lapply(c("Train","Test"),function(set_name){
      ii<-which(Set==set_name)
      data.frame(
        Stage=paste0("S",s),TimePoint=colnames(Y)[j],Set=set_name,
        Model=c("Residual","G+Residual"),
        Correlation=c(
          cor(Y[ii,j],pred_Res_S[[s]][ii,j],use="complete.obs"),
          cor(Y[ii,j],pred_GRes_S[[s]][ii,j],use="complete.obs")),
        RMSE=c(
          sqrt(mean((Y[ii,j]-pred_Res_S[[s]][ii,j])^2,na.rm=TRUE)),
          sqrt(mean((Y[ii,j]-pred_GRes_S[[s]][ii,j])^2,na.rm=TRUE)))
      )
    }))
  ))
))

Predictions_stage<-do.call(rbind,lapply(1:12,function(s)
  data.frame(
    id=rep(id,ncol(Y)),Set=rep(Set,ncol(Y)),Stage=paste0("S",s),
    TimePoint=rep(colnames(Y),each=length(id)),
    Observed=as.vector(Y),
    Predicted_Residual=as.vector(pred_Res_S[[s]]),
    Predicted_GResidual=as.vector(pred_GRes_S[[s]])
  )
))

stage_summary<-aggregate(Correlation~Stage+Model,
                         subset(accuracy_stage,Set=="Test"),
                         mean,na.rm=TRUE)

write.csv(accuracy_stage,"Stage_specific_residual_accuracy.csv",row.names=FALSE)
write.csv(Predictions_stage,"Stage_specific_residual_predictions.csv",row.names=FALSE)
write.csv(stage_summary,"Stage_specific_residual_mean_test_accuracy.csv",row.names=FALSE)

save(Pres_S,pred_Res_S,pred_GRes_S,
     file="Stage_specific_residual_models.RData")

print(stage_summary[order(-stage_summary$Correlation),])