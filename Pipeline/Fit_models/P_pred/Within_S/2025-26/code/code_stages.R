
 rm(list = ls())

 getwd()
 setwd("../Desktop/PhD data/P3/Pipeline/Fit_models/P_pred/Within_S/2025-26/output")

 library(BGLR)

 set.seed(2526)
 nIter  <- 12000
 burnIn <- 2000


 Pheno_std <- read.csv("../../../../../3.Phenomic_profiling/output/Final_outs1/Pheno_std_matched.csv")
 Res <- read.csv("../../../../../3.Phenomic_profiling/output/Final_outs1/Res1_matched.csv")

 load("../../../../../4.G/output/ZGZt.rda")
 load("../../../../../6.P/output/P_std.rda")
 load("Predicted_phenomics_and_Phat_matrices.RData")  # loads pred_X
 dim(pred_X)

  G <- ZGZt
  Pheno2526 <- Pheno_std[Pheno_std$Season == "2025-26", ]
  Res2526   <- Res[Res$Season == "2025-26", ]

  Res2526$id <- paste(Res2526$Season, Res2526$Genotype, sep = "_")
  Pheno2526$id <- paste(Pheno2526$Season, Pheno2526$Genotype, sep = "_")

  id <- sort(Reduce(intersect, list(
  Res2526$id, Pheno2526$id, rownames(G), rownames(P_std)
  )))
 
  Res2526 <- Res2526[match(id, Res2526$id), ]
  Y <- as.matrix(Res2526[, grep("_cum$", names(Res2526))])
  rownames(Y) <- id

  G <- G[id, id, drop = FALSE]
  G <- G / mean(diag(G))
  
  stopifnot(nrow(pred_X) == length(id))
  rownames(pred_X) <- id
  


  test  <- sample(seq_len(nrow(Y)), round(0.30 * nrow(Y)))
  train <- setdiff(seq_len(nrow(Y)), test)
  Set <- rep("Train", nrow(Y)); Set[test] <- "Test"



  Predicted_UAV <- data.frame(id = id, Set = Set, pred_X)
  write.csv(Predicted_UAV, "Predicted_all_UAV_traits.csv", row.names = FALSE)


   
  Phat_S <- vector("list", 12)
  names(Phat_S) <- paste0("S", 1:12)

  for (s in 1:12) {
    cols <- grep(paste0("^S", s, "_"), colnames(pred_X))
    Zs <- pred_X[, cols, drop = FALSE]
 
    mu <- colMeans(Zs[train, , drop = FALSE])
    sd_s <- apply(Zs[train, , drop = FALSE], 2, sd)
    keep <- is.finite(mu) & is.finite(sd_s) & sd_s > 0
    Zs <- scale(Zs[, keep, drop = FALSE], mu[keep], sd_s[keep])

    Phat_S[[s]] <- tcrossprod(Zs) / ncol(Zs)
    rownames(Phat_S[[s]]) <- colnames(Phat_S[[s]]) <- id
    Phat_S[[s]] <- Phat_S[[s]] / mean(diag(Phat_S[[s]]))
  }

  lapply(Phat_S, dim)


  pred_Phat_S <- pred_GPhat_S <- vector("list", 12)

  for (s in 1:12) {
    pred_Phat_S[[s]] <- Y * NA_real_
    pred_GPhat_S[[s]] <- Y * NA_real_

    for (j in seq_len(ncol(Y))) {
      y <- Y[, j]; y[test] <- NA

      fit_1 <- BGLR(
        y, ETA = list(list(K = Phat_S[[s]], model = "RKHS")),
        nIter = nIter, burnIn = burnIn, verbose = FALSE
        )

      fit_2 <- BGLR(
        y, ETA = list(list(K = G, model = "RKHS"),
                      list(K = Phat_S[[s]], model = "RKHS")),
        nIter = nIter, burnIn = burnIn, verbose = FALSE
      )

      pred_Phat_S[[s]][, j] <- fit_1$yHat
      pred_GPhat_S[[s]][, j] <- fit_2$yHat
    }

  save(Phat_S, pred_Phat_S, pred_GPhat_S,
       file = "Stage_specific_models_checkpoint.RData")
  }

 
  accuracy_stage <- data.frame()

  for (s in 1:12) {
    for (j in seq_len(ncol(Y))) {
      for (set_name in c("Train", "Test")) {
        ii <- which(Set == set_name)

      accuracy_stage <- rbind(accuracy_stage, data.frame(
        Stage = paste0("S", s),
        TimePoint = colnames(Y)[j],
        Set = set_name,
        Model = c("Phat", "G+Phat"),
        Correlation = c(
          cor(Y[ii, j], pred_Phat_S[[s]][ii, j], use = "complete.obs"),
          cor(Y[ii, j], pred_GPhat_S[[s]][ii, j], use = "complete.obs")
        ),
        RMSE = c(
          sqrt(mean((Y[ii, j] - pred_Phat_S[[s]][ii, j])^2, na.rm = TRUE)),
          sqrt(mean((Y[ii, j] - pred_GPhat_S[[s]][ii, j])^2, na.rm = TRUE))
         )
       ))
    }
   }
  }



  Predictions_stage <- do.call(rbind, lapply(1:12, function(s) {
    data.frame(
    id = rep(id, ncol(Y)),
    Set = rep(Set, ncol(Y)),
    Stage = paste0("S", s),
    TimePoint = rep(colnames(Y), each = length(id)),
    Observed = as.vector(Y),
    Predicted_Phat = as.vector(pred_Phat_S[[s]]),
    Predicted_GPhat = as.vector(pred_GPhat_S[[s]])
  )
  }))

  stage_summary <- aggregate(
    Correlation ~ Stage + Model,
    data = subset(accuracy_stage, Set == "Test"),
    FUN = mean
  )

  write.csv(accuracy_stage, "Stage_specific_accuracy.csv", row.names = FALSE)
  write.csv(Predictions_stage, "Stage_specific_predictions.csv", row.names = FALSE)
  write.csv(stage_summary, "Stage_specific_mean_test_accuracy.csv", row.names = FALSE)

  save(Phat_S, pred_Phat_S, pred_GPhat_S,
     file = "Stage_specific_Phat_models.RData")

  print(stage_summary[order(-stage_summary$Correlation), ])