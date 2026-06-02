rm(list=ls())

getwd()
setwd('../../4.G/output/')

X <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Geno_matched.csv')
Pheno_orig <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Pheno_orig_matched.csv')
Pheno_std <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Pheno_std_matched.csv')
Res <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Res1_matched.csv')

for(i in 1:ncol(X)){
  meanXi <- mean(X[, i], na.rm = TRUE)
  X[, i] <- ifelse(is.na(X[, i]), meanXi, X[, i])
}
rownames(X) <- X$X
 X <- X[, -1]
 str(X)
 
p     <- colMeans(X, na.rm = TRUE) / 2
p     <- ifelse(p <= 0.5, p, 1 - p)
index <- which(p >= 0.05)

cat("Markers before MAF:", ncol(X), "\n")
cat("Markers after  MAF:", length(index), "\n")

X <- X[, index]
X <- scale(X)

any(duplicated(rownames(X)))

X <- as.matrix(X)
  
G <- tcrossprod(X) / ncol(X)
rownames(G) <- colnames(G) <- rownames(X)

cat("G dimensions:      ", dim(G), "\n")
cat("Mean diagonal of G:", round(mean(diag(G)), 3), "\n")
cat("Max  diagonal of G:", round(max(diag(G)),  3), "\n")

uid <- paste(Pheno_std$Season, Pheno_std$Genotype, sep = "_")
names(uid) <- seq_along(uid)

rownames(Pheno_std)  <- uid
rownames(Pheno_orig) <- uid
rownames(Res)       <- uid
all.equal(Pheno_std$Genotype, Pheno_orig$Genotype)

# Check
head(uid)
sum(duplicated(uid))


obs_genos <- Pheno_std$Genotype
uni_genos <- rownames(G)

IDs <- factor(obs_genos, levels = uni_genos)
Z   <- as.matrix(model.matrix(~ IDs - 1))
colnames(Z) <- uni_genos
rownames(Z) <- paste(Pheno_std$Season, Pheno_std$Genotype, sep = "_")

cat("Z dimensions:", dim(Z), "\n")


ZGZt <- tcrossprod(tcrossprod(Z, G), Z)

rownames(ZGZt) <- paste(Pheno_std$Season, Pheno_std$Genotype, sep = "_")
colnames(ZGZt) <- paste(Pheno_std$Season, Pheno_std$Genotype, sep = "_")

# --- 4.6: EVD on G ---
EVD <- eigen(ZGZt)
rownames(EVD$vectors) <- rownames(ZGZt)

# --- 4.7: PCA plots ---
par(mfrow = c(2, 3))
plot(EVD$vectors[, 1], EVD$vectors[, 2], xlab = "PC1", ylab = "PC2")
plot(EVD$vectors[, 1], EVD$vectors[, 3], xlab = "PC1", ylab = "PC3")
plot(EVD$vectors[, 2], EVD$vectors[, 3], xlab = "PC2", ylab = "PC3")
plot(EVD$values,       xlab = "Component", ylab = "Eigenvalue")
plot(cumsum(EVD$values)[1:100] / sum(EVD$values),
     xlab = "Components", ylab = "Cumulative Variance")
abline(h = 0.8, lty = 3)
par(mfrow = c(1, 1))

mean(diag(G))
mean(diag(ZGZt))
# --- 4.8: Save ---

# G
save(G, file = "G.rda")
saveRDS(G,  "G_matched.rds")
write.csv(G,  "G_matched.csv", row.names = TRUE)

# ZGZt
save(ZGZt, file = "ZGZt.rda")
saveRDS(ZGZt, "ZGZt_matched.rds")
# write.csv(ZGZt, "ZGZt_matched.csv", row.names = TRUE)  # skip - very large

# Z
#save(Z, file = "Z.rda")
#saveRDS(Z,"Z_matched.rds")
#write.csv(Z,    "Z_matched.csv",    row.names = TRUE)

# EVD
save(EVD, file = "EVD.rda")

rownames(G)[1:3]
rownames(ZGZt)[1:3]
colnames(ZGZt)[1:3]
rownames(Z)[1:3]
rownames(EVD$vectors)[1:3]
