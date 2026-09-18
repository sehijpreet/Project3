  rm(list = ls())
  getwd()
  setwd('../Desktop/PhD data/P3/Pipeline/Fit_models/P_pred/Within_S/2025-26/output')
  setwd('../../../Fit_models/MegaLMM/output/')

#X <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Geno_matched.csv')
  Pheno_orig <- read.csv('../../../../../3.Phenomic_profiling/output/Final_outs1/Pheno_orig_matched.csv')
  Pheno_std <- read.csv('../../../../../3.Phenomic_profiling/output/Final_outs1/Pheno_std_matched.csv')
  Res <- read.csv('../../../../../3.Phenomic_profiling/output/Final_outs1/Res1_matched.csv')
 X          <- read.csv("../../../3.Phenomic_profiling/output/Final_outs1/Geno_matched.csv")

 load("../../../../../4.G/output/ZGZt.rda")
 load("../../../../../6.P/output/P_std.rda")
 G <- ZGZt

## Extract the 2025-26 season
  Pheno2526 <- Pheno_std[Pheno_std$Season == "2025-26", ]
  Res2526   <- Res[Res$Season == "2025-26", ]

  set.seed(2526)
  nIter  <- 12000
  burnIn <- 2000

## ---- 1_prepare_yield_data ---------------------------------------------------
  Res2526$id <- paste(Res2526$Season, Res2526$Genotype, sep = "_")
  Pheno2526$id <- paste(Pheno2526$Season, Pheno2526$Genotype, sep = "_")

## Keep only IDs present in yield, UAV, G, and P
  id <- sort(Reduce(intersect, list(
  Res2526$id, Pheno2526$id, rownames(G), rownames(P_std)
  )))
 
  Res2526   <- Res2526[match(id, Res2526$id), ]
  Pheno2526 <- Pheno2526[match(id, Pheno2526$id), ]

  Y <- as.matrix(Res2526[, grep("_cum$", names(Res2526))])
  rownames(Y) <- id
  
## ---- 2_prepare_G_P_UAV ------------------------------------------------------
  G <- G[id, id, drop = FALSE]
  P <- P_std[id, id, drop = FALSE]

  G <- G / mean(diag(G))
  P <- P / mean(diag(P))

  stopifnot(identical(id, Pheno2526$id), nrow(Y) == nrow(G), nrow(Y) == nrow(P))

## ---- 3_create_train_test_split ---------------------------------------------
  test  <- sample(seq_len(nrow(Y)), round(0.30 * nrow(Y)))
  train <- setdiff(seq_len(nrow(Y)), test)
  Set <- rep("Train", nrow(Y)); Set[test] <- "Test"
  library(BGLR)

## ---- 4_models_G_P_GP --------------------------------------------------------

  pred_G <- pred_P <- pred_GP <- Y * NA_real_

  for (j in seq_len(ncol(Y))) {
   y <- Y[, j]; y[test] <- NA

   fit_G  <- BGLR(y, ETA = list(list(K = G, model = "RKHS")),
                 nIter = nIter, burnIn = burnIn, verbose = FALSE)
   fit_P  <- BGLR(y, ETA = list(list(K = P, model = "RKHS")),
                 nIter = nIter, burnIn = burnIn, verbose = FALSE)
   fit_GP <- BGLR(y, ETA = list(list(K = G, model = "RKHS"),
                               list(K = P, model = "RKHS")),
                 nIter = nIter, burnIn = burnIn, verbose = FALSE)

   pred_G[, j]  <- fit_G$yHat
   pred_P[, j]  <- fit_P$yHat
   pred_GP[, j] <- fit_GP$yHat
  }

## ---- 5_stagewise_PCA --------------------------------------------------------
  PC_stage <- NULL

  for (s in 1:12) {
   Xs <- as.matrix(Pheno2526[, grep(paste0("^S", s, "_"), names(Pheno2526))])
   Xs[!is.finite(Xs)] <- NA
   mu <- colMeans(Xs[train, , drop = FALSE], na.rm = TRUE)
   sd_x <- apply(Xs[train, , drop = FALSE], 2, sd, na.rm = TRUE)
   keep <- is.finite(mu) & is.finite(sd_x) & sd_x > 0
   Xs <- Xs[, keep, drop = FALSE]; mu <- mu[keep]
   for (m in seq_len(ncol(Xs))) Xs[is.na(Xs[, m]), m] <- mu[m]

   pca <- prcomp(Xs[train, , drop = FALSE], center = TRUE, scale. = TRUE)
   scores <- predict(pca, newdata = Xs)[, 1:5, drop = FALSE]
   colnames(scores) <- paste0("S", s, "_PC", 1:5)
   PC_stage <- cbind(PC_stage, scores)
 }

  rownames(PC_stage) <- id

## ---- 6_predict_stagewise_PCs ------------------------------------------------

  pred_PC_stage <- PC_stage * NA_real_

  for (j in seq_len(ncol(PC_stage))) {
    y <- PC_stage[, j]; y[test] <- NA
    fit <- BGLR(y, ETA = list(list(K = G, model = "RKHS")),
              nIter = nIter, burnIn = burnIn, verbose = FALSE)
    pred_PC_stage[, j] <- fit$yHat
  }

  Z <- scale(pred_PC_stage,
           center = colMeans(pred_PC_stage[train, , drop = FALSE]),
           scale = apply(pred_PC_stage[train, , drop = FALSE], 2, sd))
  Phat_stage <- tcrossprod(Z) / ncol(Z)
  rownames(Phat_stage) <- colnames(Phat_stage) <- id
  Phat_stage <- Phat_stage / mean(diag(Phat_stage))

## ---- 7_models_stagewise_Phat ------------------------------------------------
 
  pred_Phat_stage <- pred_GPhat_stage <- Y * NA_real_

  for (j in seq_len(ncol(Y))) {
  y <- Y[, j]; y[test] <- NA
  fit_1 <- BGLR(y, ETA = list(list(K = Phat_stage, model = "RKHS")),
                nIter = nIter, burnIn = burnIn, verbose = FALSE)
  fit_2 <- BGLR(y, ETA = list(list(K = G, model = "RKHS"),
                              list(K = Phat_stage, model = "RKHS")),
                nIter = nIter, burnIn = burnIn, verbose = FALSE)
  pred_Phat_stage[, j]  <- fit_1$yHat
  pred_GPhat_stage[, j] <- fit_2$yHat
}

## ---- 8_global_PCA -----------------------------------------------------------
X <- as.matrix(Pheno2526[, grep("^S[0-9]+_", names(Pheno2526))])
X[!is.finite(X)] <- NA
mu <- colMeans(X[train, , drop = FALSE], na.rm = TRUE)
sd_x <- apply(X[train, , drop = FALSE], 2, sd, na.rm = TRUE)
keep <- is.finite(mu) & is.finite(sd_x) & sd_x > 0
X <- X[, keep, drop = FALSE]; mu <- mu[keep]; sd_x <- sd_x[keep]
for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- mu[j]

pca_all <- prcomp(X[train, , drop = FALSE], center = TRUE, scale. = TRUE)
nPC <- min(60, ncol(pca_all$rotation))
PC_all <- predict(pca_all, newdata = X)[, 1:nPC, drop = FALSE]
rownames(PC_all) <- id; colnames(PC_all) <- paste0("PC", 1:nPC)

## ---- 9_predict_global_PCs ---------------------------------------------------
pred_PC_all <- PC_all * NA_real_

for (j in seq_len(ncol(PC_all))) {
  y <- PC_all[, j]; y[test] <- NA
  fit <- BGLR(y, ETA = list(list(K = G, model = "RKHS")),
              nIter = nIter, burnIn = burnIn, verbose = FALSE)
  pred_PC_all[, j] <- fit$yHat
}

Z <- scale(pred_PC_all,
           center = colMeans(pred_PC_all[train, , drop = FALSE]),
           scale = apply(pred_PC_all[train, , drop = FALSE], 2, sd))
Phat_all <- tcrossprod(Z) / ncol(Z)
rownames(Phat_all) <- colnames(Phat_all) <- id
Phat_all <- Phat_all / mean(diag(Phat_all))

## ---- 10_models_global_Phat --------------------------------------------------
pred_Phat_all <- pred_GPhat_all <- Y * NA_real_

for (j in seq_len(ncol(Y))) {
  y <- Y[, j]; y[test] <- NA
  fit_1 <- BGLR(y, ETA = list(list(K = Phat_all, model = "RKHS")),
                nIter = nIter, burnIn = burnIn, verbose = FALSE)
  fit_2 <- BGLR(y, ETA = list(list(K = G, model = "RKHS"),
                              list(K = Phat_all, model = "RKHS")),
                nIter = nIter, burnIn = burnIn, verbose = FALSE)
  pred_Phat_all[, j]  <- fit_1$yHat
  pred_GPhat_all[, j] <- fit_2$yHat
}

## ---- 11_predict_all_UAV_traits ---------------------------------------------
## X is already imputed; standardize it using training-set parameters.
X_std <- scale(X, center = mu, scale = sd_x)
pred_X <- X_std * NA_real_

for (j in seq_len(ncol(X_std))) {
  y <- X_std[, j]; y[test] <- NA
  fit <- BGLR(y, ETA = list(list(K = G, model = "RKHS")),
              nIter = nIter, burnIn = burnIn, verbose = FALSE)
  pred_X[, j] <- fit$yHat
}

Z <- scale(pred_X,
           center = colMeans(pred_X[train, , drop = FALSE]),
           scale = apply(pred_X[train, , drop = FALSE], 2, sd))
keep <- apply(Z, 2, function(z) all(is.finite(z)))
Z <- Z[, keep, drop = FALSE]
Phat_traits <- tcrossprod(Z) / ncol(Z)
rownames(Phat_traits) <- colnames(Phat_traits) <- id
Phat_traits <- Phat_traits / mean(diag(Phat_traits))

## ---- 12_models_all_trait_Phat ----------------------------------------------
pred_Phat_traits <- pred_GPhat_traits <- Y * NA_real_

for (j in seq_len(ncol(Y))) {
  y <- Y[, j]; y[test] <- NA
  fit_1 <- BGLR(y, ETA = list(list(K = Phat_traits, model = "RKHS")),
                nIter = nIter, burnIn = burnIn, verbose = FALSE)
  fit_2 <- BGLR(y, ETA = list(list(K = G, model = "RKHS"),
                              list(K = Phat_traits, model = "RKHS")),
                nIter = nIter, burnIn = burnIn, verbose = FALSE)
  pred_Phat_traits[, j]  <- fit_1$yHat
  pred_GPhat_traits[, j] <- fit_2$yHat
}

## ---- 13_combine_accuracy ----------------------------------------------------
model_names <- c("G", "P", "G+P", "Phat_stage", "G+Phat_stage",
                 "Phat_all", "G+Phat_all", "Phat_traits", "G+Phat_traits")
pred_list <- list(pred_G, pred_P, pred_GP, pred_Phat_stage, pred_GPhat_stage,
                  pred_Phat_all, pred_GPhat_all, pred_Phat_traits,
                  pred_GPhat_traits)

accuracy <- data.frame()

for (j in seq_len(ncol(Y))) {
  for (s in c("Train", "Test")) {
    ii <- which(Set == s)
    accuracy <- rbind(accuracy, data.frame(
      TimePoint = colnames(Y)[j], Set = s, Model = model_names,
      Correlation = sapply(pred_list, function(x)
        cor(Y[ii, j], x[ii, j], use = "complete.obs")),
      RMSE = sapply(pred_list, function(x)
        sqrt(mean((Y[ii, j] - x[ii, j])^2, na.rm = TRUE)))
    ))
  }
}

## ---- 14_save_predictions_and_results ---------------------------------------
Predictions <- data.frame(
  id = id, Set = Set, Y,
  setNames(as.data.frame(pred_G), paste0(colnames(Y), "_G")),
  setNames(as.data.frame(pred_P), paste0(colnames(Y), "_P")),
  setNames(as.data.frame(pred_GP), paste0(colnames(Y), "_GP")),
  setNames(as.data.frame(pred_Phat_stage), paste0(colnames(Y), "_Phat_stage")),
  setNames(as.data.frame(pred_GPhat_stage), paste0(colnames(Y), "_GPhat_stage")),
  setNames(as.data.frame(pred_Phat_all), paste0(colnames(Y), "_Phat_all")),
  setNames(as.data.frame(pred_GPhat_all), paste0(colnames(Y), "_GPhat_all")),
  setNames(as.data.frame(pred_Phat_traits), paste0(colnames(Y), "_Phat_traits")),
  setNames(as.data.frame(pred_GPhat_traits), paste0(colnames(Y), "_GPhat_traits"))
)

write.csv(Predictions, "All_model_predictions.csv", row.names = FALSE)
write.csv(accuracy, "All_model_accuracy.csv", row.names = FALSE)
save(PC_stage, PC_all, pred_PC_stage, pred_PC_all, pred_X,
     Phat_stage, Phat_all, Phat_traits,
     file = "Predicted_phenomics_and_Phat_matrices.RData")

print(subset(accuracy, Set == "Test"))