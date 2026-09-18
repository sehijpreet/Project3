rm(list=ls())
library(BGLR)
getwd()


season <- "2020-21"
base <- "../Desktop/PhD data/P3/Pipeline/Fit_models/P_pred/Within_S"
getwd()
setwd('../../2020-21/output')
set.seed(2526); nIter<-12000; burnIn<-2000

Pheno_std <- read.csv("../../../../../3.Phenomic_profiling/output/Final_outs1/Pheno_std_matched.csv")
Res <- read.csv("../../../../../3.Phenomic_profiling/output/Final_outs1/Res1_matched.csv")
load("../../../../../4.G/output/ZGZt.rda")
load("../../../../../6.P/output/P_std.rda")

Pheno <- Pheno_std[Pheno_std$Season==season,]
ResS <- Res[Res$Season==season,]
Pheno$id <- paste(Pheno$Season,Pheno$Genotype,sep="_")
ResS$id <- paste(ResS$Season,ResS$Genotype,sep="_")

id <- sort(Reduce(intersect,list(Pheno$id,ResS$id,rownames(ZGZt),rownames(P_std))))
Pheno <- Pheno[match(id,Pheno$id),]
ResS <- ResS[match(id,ResS$id),]

Y <- as.matrix(ResS[,grep("_cum$",names(ResS))]); rownames(Y)<-id
G <- ZGZt[id,id]; G <- G/mean(diag(G))
P <- P_std[id,id]; P <- P/mean(diag(P))

set.seed(2526)
test <- sample(seq_len(nrow(Y)),round(.30*nrow(Y)))
train <- setdiff(seq_len(nrow(Y)),test)
Set <- ifelse(seq_len(nrow(Y))%in%test,"Test","Train")

# ---------- predict all phenomic/UAV traits from G ----------
traits <- setdiff(names(Pheno),c("Season","Genotype","id"))
traits <- traits[sapply(Pheno[,traits,drop=FALSE],is.numeric)]

pred_X <- matrix(NA,nrow(Y),length(traits),dimnames=list(id,traits))

for(j in seq_along(traits)){
  y <- Pheno[[traits[j]]]; y[test] <- NA
  if(sum(!is.na(y[train]))<10) next
  fit <- BGLR(y,ETA=list(G=list(K=G,model="RKHS")),
              nIter=nIter,burnIn=burnIn,verbose=FALSE)
  pred_X[,j] <- fit$yHat
}

# ---------- P-hat from predicted phenomics ----------
keep <- apply(pred_X,2,function(x) all(is.finite(x)) && sd(x[train])>0)
Z <- pred_X[,keep,drop=FALSE]
mu <- colMeans(Z[train,,drop=FALSE])
sdv <- apply(Z[train,,drop=FALSE],2,sd)
Z <- scale(Z,mu,sdv)

Phat <- tcrossprod(Z)/ncol(Z)
rownames(Phat)<-colnames(Phat)<-id
Phat <- Phat/mean(diag(Phat))

# ---------- fit basic yield models ----------
mods <- list(
  G=list(G=list(K=G,model="RKHS")),
  P=list(P=list(K=P,model="RKHS")),
  `G+P`=list(G=list(K=G,model="RKHS"),P=list(K=P,model="RKHS")),
  Phat=list(Phat=list(K=Phat,model="RKHS")),
  `G+Phat`=list(G=list(K=G,model="RKHS"),Phat=list(K=Phat,model="RKHS"))
)

Pred <- lapply(mods,function(x) Y*NA_real_)

for(m in names(mods)){
  cat("\n",m,"\n")
  for(j in seq_len(ncol(Y))){
    y<-Y[,j]; y[test]<-NA
    fit<-BGLR(y,ETA=mods[[m]],nIter=nIter,burnIn=burnIn,verbose=FALSE)
    Pred[[m]][,j]<-fit$yHat
  }
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

write.csv(data.frame(id=id,Set=Set,pred_X),"Predicted_all_UAV_traits.csv",row.names=FALSE)
write.csv(acc,"Basic_models_accuracy.csv",row.names=FALSE)
write.csv(Predictions,"Basic_models_predictions.csv",row.names=FALSE)

save(pred_X,Phat,G,P,Y,id,train,test,Set,traits,
     file="Predicted_phenomics_and_Phat_matrices.RData")

print(aggregate(Correlation~Model,subset(acc,Set=="Test"),mean,na.rm=TRUE))

