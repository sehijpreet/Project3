rm(list=ls())

getwd()
setwd('../Desktop/PhD data/P3/Pipeline/3.Phenomic_profiling/output/Final_outs/')
setwd('../../../6.P/output/')

#X <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Geno_matched.csv')
Pheno_orig <- read.csv('../../3.Phenomic_profiling/output/Final_outs1/Pheno_orig_matched.csv')
Pheno_std <- read.csv('../../3.Phenomic_profiling/output/Final_outs1/Pheno_std_matched.csv')
Res <- read.csv('../../3.Phenomic_profiling/output/Final_outs1/Res1_matched.csv')

load('../../4.G/output/G.rda')
load('../../4.G/output/ZGZt.rda')

G[1:5, 1:5]
ZGZt[1:5, 1:5]

sum(is.na(Pheno_std[, 3:ncol(Pheno_std)]))
sum(is.na(Pheno_orig[, 3:ncol(Pheno_orig)]))
tapply(is.na(Pheno_std$S3_CanopyArea_z),  Pheno_std$Season, mean) * 100

rownames(Pheno_std)

uid <- paste(Pheno_std$Season, Pheno_std$Genotype, sep = "_")
rownames(Pheno_std)  <- uid
rownames(Pheno_orig) <- uid
rownames(Res)        <- uid


Pheno_std_imp  <- Pheno_std
Pheno_orig_imp <- Pheno_orig

for(i in 3:ncol(Pheno_std_imp)){
  m <- mean(Pheno_std_imp[, i], na.rm = TRUE)
  Pheno_std_imp[is.na(Pheno_std_imp[, i]), i] <- m
}

for(i in 3:ncol(Pheno_orig_imp)){
  m <- mean(Pheno_orig_imp[, i], na.rm = TRUE)
  Pheno_orig_imp[is.na(Pheno_orig_imp[, i]), i] <- m
}

# Verify no NAs remain
sum(is.na(Pheno_std_imp[,  3:ncol(Pheno_std_imp)]))
sum(is.na(Pheno_orig_imp[, 3:ncol(Pheno_orig_imp)]))

#  Extract trait matrices (drop Season + Genotype cols) ---
X_std  <- as.matrix(Pheno_std_imp[,  3:ncol(Pheno_std_imp)])
X_orig <- as.matrix(Pheno_orig_imp[, 3:ncol(Pheno_orig_imp)])

# Scale and build kernels ---

# All stages combined
X_std_scaled  <- scale(X_std)
X_orig_scaled <- scale(X_orig)

P_std  <- tcrossprod(X_std_scaled)  / ncol(X_std_scaled)
P_orig <- tcrossprod(X_orig_scaled) / ncol(X_orig_scaled)
dim(P_std)
dim(P_orig)

rownames(P_std)  <- uid
colnames(P_std)  <- uid
rownames(P_orig) <- uid
colnames(P_orig) <- uid

mean(diag(P_std))
mean(diag(P_orig))

stages_list <- paste0("S", 1:12)

P_std_stages  <- vector("list", length(stages_list))
P_orig_stages <- vector("list", length(stages_list))
names(P_std_stages)  <- stages_list
names(P_orig_stages) <- stages_list

for(stg in stages_list){
  
  # Columns belonging to this stage
  stg_cols <- grep(paste0("^", stg, "_"), colnames(X_std))
  
  if(length(stg_cols) == 0) next
  
  Xs <- scale(X_std[,  stg_cols])
  Xo <- scale(X_orig[, stg_cols])
  
  Ks <- tcrossprod(Xs) / ncol(Xs)
  Ko <- tcrossprod(Xo) / ncol(Xo)
  
  rownames(Ks) <- uid;  colnames(Ks) <- uid
  rownames(Ko) <- uid;  colnames(Ko) <- uid
  
  P_std_stages[[stg]]  <- Ks
  P_orig_stages[[stg]] <- Ko
}

# Check diagonal means per stage
round(sapply(P_std_stages,  function(k) mean(diag(k))), 3)
round(sapply(P_orig_stages, function(k) mean(diag(k))), 3)

P_std_1 <- P_std
getwd()
P_std_2  <- load('P_std.rda')

all.equal(P_std_1, P_std_2)
identical(P_std_1, P_std_2)

# --- 5.7: Save ---
save(P_std,        file = "P_std.rda")
save(P_orig,       file = "P_orig.rda")
save(P_std_stages, file = "P_std_stages.rda")
save(P_orig_stages, file = "P_orig_stages.rda")


#############

# --- EVD for full phenomic kernels ---
EVD_P_std  <- eigen(P_std)
EVD_P_orig <- eigen(P_orig)

rownames(EVD_P_std$vectors)  <- uid
rownames(EVD_P_orig$vectors) <- uid

# --- EVD for per stage kernels ---
EVD_P_std_stages  <- vector("list", length(stages_list))
EVD_P_orig_stages <- vector("list", length(stages_list))
names(EVD_P_std_stages)  <- stages_list
names(EVD_P_orig_stages) <- stages_list

for(stg in stages_list){
  
  if(is.null(P_std_stages[[stg]])) next
  
  evd_s <- eigen(P_std_stages[[stg]])
  evd_o <- eigen(P_orig_stages[[stg]])
  
  rownames(evd_s$vectors) <- uid
  rownames(evd_o$vectors) <- uid
  
  EVD_P_std_stages[[stg]]  <- evd_s
  EVD_P_orig_stages[[stg]] <- evd_o
}

# --- Verify ---
names(EVD_P_std_stages)
rownames(EVD_P_std$vectors)[1:3]
rownames(EVD_P_std_stages[["S1"]]$vectors)[1:3]

# --- Save ---
save(EVD_P_std,        file = "EVD_P_std.rda")
save(EVD_P_orig,       file = "EVD_P_orig.rda")
save(EVD_P_std_stages, file = "EVD_P_std_stages.rda")
save(EVD_P_orig_stages,file = "EVD_P_orig_stages.rda")

