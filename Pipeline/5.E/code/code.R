rm(list=ls())

getwd()
setwd('../../5.E/output/')

#X <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Geno_matched.csv')
Pheno_orig <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Pheno_orig_matched.csv')
Pheno_std <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Pheno_std_matched.csv')
Res <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Res1_matched.csv')


uid <- paste(Pheno_std$Season, Pheno_std$Genotype, sep = "_")

# --- Environment factor ---
EID <- as.factor(Res$Season)
levels(EID)
length(levels(EID))   

Z_E <- as.matrix(model.matrix(~ EID - 1))
rownames(Z_E) <- uid
colnames(Z_E) <- levels(EID)

dim(Z_E)   # n_obs x 6


E <- diag(length(levels(EID)))
rownames(E) <- levels(EID)
colnames(E) <- levels(EID)

# --- ZEZ' ---
ZEZt <- tcrossprod(tcrossprod(Z_E, E), Z_E)
rownames(ZEZt) <- uid
colnames(ZEZt) <- uid

dim(ZEZt)
mean(diag(ZEZt))
ZEZt[1:5, 1:5]

# --- EVD ---
EVD_E <- eigen(ZEZt)
rownames(EVD_E$vectors) <- uid

# --- Save ---
save(Z_E,   file = "Z_E.rda")
save(E,     file = "E.rda")
save(ZEZt,  file = "ZEZt.rda")
save(EVD_E, file = "EVD_E.rda")

# --- Check ---
rownames(ZEZt)[1:3]
colnames(ZEZt)[1:3]
rownames(EVD_E$vectors)[1:3]

