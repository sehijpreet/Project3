rm(list=ls())
setwd("../Desktop/PhD data/P3/Pipeline/Fit_models/P_pred/Within_S/2025-26/output")
getwd()
dir.create('Res', showWarnings = FALSE)
setwd('./Res')
library(BGLR)

set.seed(2526)
nIter <- 12000; burnIn <- 2000

Pheno_std <- read.csv("../../../../../../3.Phenomic_profiling/output/Final_outs1/Pheno_std_matched.csv")
Res <- read.csv("../../../../../../3.Phenomic_profiling/output/Final_outs1/Res1_matched.csv")
load("../../../../../../4.G/output/ZGZt.rda")
load("../Predicted_phenomics_and_Phat_matrices.RData")   # pred_X

# ---------- Match 2025-26 individuals ----------
Pheno <- Pheno_std[Pheno_std$Season=="2025-26",]
Res25 <- Res[Res$Season=="2025-26",]
Pheno$id <- paste(Pheno$Season,Pheno$Genotype,sep="_")
Res25$id <- paste(Res25$Season,Res25$Genotype,sep="_")

id <- sort(Reduce(intersect,list(Pheno$id,Res25$id,rownames(ZGZt))))
Pheno <- Pheno[match(id,Pheno$id),]
Res25 <- Res25[match(id,Res25$id),]

Y <- as.matrix(Res25[,grep("_cum$",names(Res25))]); rownames(Y) <- id
G <- ZGZt[id,id]; G <- G/mean(diag(G))

stopifnot(nrow(pred_X)==length(id))
rownames(pred_X) <- id

# ---------- Same 70/30 split ----------
set.seed(2526)
test <- sample(seq_len(nrow(Y)),round(.30*nrow(Y)))
train <- setdiff(seq_len(nrow(Y)),test)
Set <- ifelse(seq_len(nrow(Y)) %in% test,"Test","Train")

# ---------- Observed - predicted UAV ----------
traits <- intersect(colnames(pred_X),names(Pheno))
Xobs <- as.matrix(Pheno[,traits,drop=FALSE])
Xpred <- as.matrix(pred_X[,traits,drop=FALSE])
storage.mode(Xobs) <- storage.mode(Xpred) <- "numeric"

Xres <- Xobs-Xpred

# remove useless traits and standardize
s <- apply(Xres,2,sd,na.rm=TRUE)
Xres <- Xres[,is.finite(s) & s>0,drop=FALSE]
Xres <- scale(Xres)

# mean-impute remaining missing residuals
for(j in seq_len(ncol(Xres))){
  ii <- is.na(Xres[,j])
  if(any(ii)) Xres[ii,j] <- mean(Xres[,j],na.rm=TRUE)
}

write.csv(data.frame(id=id,Set=Set,Xres),
          "Residual_UAV_traits_obs_minus_pred.csv",row.names=FALSE)

# ---------- Residual phenomic kernel ----------
P_resid <- tcrossprod(Xres)/ncol(Xres)
rownames(P_resid) <- colnames(P_resid) <- id
P_resid <- P_resid/mean(diag(P_resid))

# ---------- Fit Residual and G + Residual ----------
pred_R <- pred_GR <- matrix(NA,nrow(Y),ncol(Y),
                            dimnames=list(id,colnames(Y)))

for(j in seq_len(ncol(Y))){
  
  cat("Trait",j,"/",ncol(Y),colnames(Y)[j],"\n")
  y <- Y[,j]; y[test] <- NA
  
  fit_R <- BGLR(y,
                ETA=list(R=list(K=P_resid,model="RKHS")),
                nIter=nIter,burnIn=burnIn,verbose=FALSE)
  
  fit_GR <- BGLR(y,
                 ETA=list(G=list(K=G,model="RKHS"),
                          R=list(K=P_resid,model="RKHS")),
                 nIter=nIter,burnIn=burnIn,verbose=FALSE)
  
  pred_R[,j] <- fit_R$yHat
  pred_GR[,j] <- fit_GR$yHat
}

# ---------- Accuracy ----------
acc <- function(pred,model){
  do.call(rbind,lapply(c("Train","Test"),function(s){
    ii <- which(Set==s)
    data.frame(
      Trait=colnames(Y),Set=s,Model=model,
      Correlation=sapply(seq_len(ncol(Y)),
                         function(j) cor(Y[ii,j],pred[ii,j],use="complete.obs")),
      RMSE=sapply(seq_len(ncol(Y)),
                  function(j) sqrt(mean((Y[ii,j]-pred[ii,j])^2,na.rm=TRUE)))
    )
  }))
}

accuracy <- rbind(acc(pred_R,"Residual"),
                  acc(pred_GR,"G+Residual"))

# ---------- Save predictions ----------
Predictions <- data.frame(
  id=id,Set=Set,
  setNames(as.data.frame(Y),paste0(colnames(Y),"_Observed")),
  setNames(as.data.frame(pred_R),paste0(colnames(Y),"_Residual")),
  setNames(as.data.frame(pred_GR),paste0(colnames(Y),"_G_Residual"))
)

write.csv(accuracy,"Accuracy_Residual_GResidual.csv",row.names=FALSE)
write.csv(Predictions,"Predictions_Residual_GResidual.csv",row.names=FALSE)

print(accuracy)
