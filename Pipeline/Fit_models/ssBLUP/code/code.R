
rm(list=ls())
setwd('../output/')
library(BGLR)
library(parallel)

# --- Load ---
Res <- read.csv('../../../3.Phenomic_profiling/output/Final_outs/Res1_matched.csv')
load('../../../6.P/output/P_std.rda')
load('../../../4.G/output/ZGZt.rda')
tmp <- load('../../../5.E/output/EVD_E.rda'); eE <- get(tmp); rm(list=tmp)
tmp <- load('../../../7.L/output/EVD_L.rda'); eL <- get(tmp); rm(list=tmp)

# --- Setup ---
uid           <- rownames(P_std)
rownames(Res) <- uid
tr            <- which(Res$Season != "2025-26")
pr            <- which(Res$Season == "2025-26")
traits        <- colnames(Res)[3:13]
alphas        <- seq(0, 1, by=0.1)

dir.create("2025-26_pred", showWarnings=FALSE)

# --- Build H ---
buildH <- function(a, P, G, pr){
  H      <- a*P + (1-a)*G
  H[pr,] <- G[pr,]
  H[,pr] <- G[,pr]
  H
}

# --- Fit all traits for one alpha ---
run_a <- function(a){

  cat("Alpha:", a, "\n")

  # Build H + EVD once per alpha (not per trait)
  H  <- buildH(a, P_std, ZGZt, pr)
  eH <- eigen(H)
  rownames(eH$vectors) <- uid
  rm(H)

  dir.create(paste0("alpha_", a),              showWarnings=FALSE)
  dir.create(paste0("alpha_", a, "/pred"),     showWarnings=FALSE)

  ETA <- list(E = list(V=eE$vectors, d=eE$values, model='RKHS'),
    L = list(V=eL$vectors, d=eL$values, model='RKHS'),
    H = list(V=eH$vectors, d=eH$values, model='RKHS'))

  res_a <- vector("list", length(traits))
  names(res_a) <- traits

  for(trait in traits){
    cat("alpha:", a, "| trait:", trait, "\n")

    y   <- Res[, trait]
    yNA <- y
    yNA[pr] <- NA
    td  <- paste0("alpha_", a, "/tmp_", trait)
    dir.create(td, showWarnings=FALSE)

    fm  <- BGLR(y=yNA, ETA=ETA, nIter=12000, burnIn=2000,  saveAt=paste0(td,"/"), verbose=FALSE)

    acc <- cor(y[tr], fm$yHat[tr], use="complete.obs")

    out <- data.frame(Season=Res$Season, Genotype=Res$Genotype, y_obs=y, y_hat=fm$yHat)
    rownames(out) <- uid

    write.csv(out[pr,], paste0("alpha_", a, "/pred/", trait, ".csv"), row.names=TRUE)

    res_a[[trait]] <- list(trait=trait, alpha=a, cor=round(acc,3), varE=round(fm$varE,4), varH=round(fm$ETA[[3]]$varU,4))

    unlink(paste0(td,"/*.dat")); unlink(td, recursive=TRUE)
    rm(fm)
  }

  # Summary for this alpha
  sum_a <- data.frame(Trait = traits, Alpha = a,
    Cor   = sapply(res_a, function(x) x$cor),
    VarH  = sapply(res_a, function(x) x$varH),
    VarE  = sapply(res_a, function(x) x$varE))
  write.csv(sum_a, paste0("alpha_", a, "/summary_alpha_", a, ".csv"), row.names=FALSE)

  return(res_a)
}

# RUN - 11 alphas x 8 cores

all_res <- mclapply(alphas, run_a, mc.cores=11)
names(all_res) <- paste0("a_", alphas)


# Correlation table: traits x alphas

cor_tab <- do.call(cbind, lapply(alphas, function(a){
  sapply(traits, function(t) all_res[[paste0("a_",a)]][[t]]$cor)
}))
rownames(cor_tab) <- traits
colnames(cor_tab) <- paste0("a_", alphas)
cor_tab <- as.data.frame(cor_tab)
cor_tab$Best_alpha <- alphas[apply(cor_tab, 1, which.max)]
cor_tab$Best_cor   <- apply(cor_tab[, paste0("a_", alphas)], 1, max)

print(cor_tab)
write.csv(cor_tab, "ssBLUP_cor_all_alphas.csv", row.names=TRUE)

