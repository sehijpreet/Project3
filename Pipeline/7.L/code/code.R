rm(list=ls())
setwd('../../7.L/output/')

Res       <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Res1_matched.csv')
Pheno_std <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Pheno_std_matched.csv')

uid <- paste(Pheno_std$Season, Pheno_std$Genotype, sep = "_")


# Grouping factor = unique genotypes
GID <- factor(Res$Genotype)
levels(GID)[1:5]
length(levels(GID))   # number of unique genotypes

Z_L <- as.matrix(model.matrix(~ GID - 1))
rownames(Z_L) <- uid
colnames(Z_L) <- levels(GID)

dim(Z_L)   # n_obs x n_unique_genotypes

# d: number of observations per genotype (how many seasons each appears)
d <- colSums(Z_L)
table(d)   
sum(d)

# V: Z scaled by sqrt(d)
V <- Z_L
for(i in 1:ncol(Z_L)){
  V[, i] <- V[, i] / sqrt(d[i])
}

# EVD object (as per BGLR workflow)
EVD_L <- list(vectors = V, values = d)

# L = ZZ' (obs x obs)
# Same genotype across seasons gets 1, different genotypes get 0
L <- tcrossprod(Z_L)
rownames(L) <- uid
colnames(L) <- uid

# Verify
dim(L)
mean(diag(L))
L[1:5, 1:5]   # same genotype different season = 1

# Save
save(Z_L,   file = "Z_L.rda")
save(L,     file = "L.rda")
save(EVD_L, file = "EVD_L.rda")


# Check
rownames(L)[1:3]
colnames(L)[1:3]
length(levels(GID))
