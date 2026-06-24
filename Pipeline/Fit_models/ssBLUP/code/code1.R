rm(list=ls())

library(BGLR)
library(parallel)

# --- Load ---
Res <- read.csv('../../../3.Phenomic_profiling/output/Final_outs/Res1_matched.csv')
load('../../../6.P/output/P_std.rda')
load('../../../4.G/output/ZGZt.rda')
tmp <- load('../../../5.E/output/EVD_E.rda'); eE <- get(tmp); rm(list=tmp)
tmp <- load('../../../7.L/output/EVD_L.rda'); eL <- get(tmp); rm(list=tmp)

# --- Setup ---
uid        <- rownames(P_std)
tr         <- which(Res$Season != '2025-26')
pr         <- which(Res$Season == '2025-26')
seasons_tr <- unique(Res$Season[tr])
traits     <- colnames(Res)[3:22]
alphas     <- seq(0, 1, by=0.1)

# --- Step 1: Leave-one-season-out CV within training to find best alpha ---
run_alpha_cv <- function(a) {
  cat('Alpha:', a, '\n')
  res_a <- lapply(traits, function(trait) {
    cat('  trait:', trait, '\n')
    y_tr    <- Res[tr, trait]
    fold_stats <- lapply(seasons_tr, function(s) {
      fold         <- which(Res$Season[tr] == s)
      # P = 0 for fold (mimics no P for 2025-26)
      P_cv         <- P_std[tr, tr]
      P_cv[fold, ] <- 0
      P_cv[, fold] <- 0
      H_cv  <- a * P_cv + (1-a) * ZGZt[tr, tr]
      eH_cv <- eigen(H_cv, symmetric=TRUE)
      ETA <- list(
        E = list(V=eE$vectors[tr,], d=eE$values, model='RKHS'),
        L = list(V=eL$vectors[tr,], d=eL$values, model='RKHS'),
        H = list(V=eH_cv$vectors,   d=eH_cv$values, model='RKHS')
      )
      y_cv       <- y_tr
      y_cv[fold] <- NA
      fm         <- BGLR(y=y_cv, ETA=ETA, nIter=12000, burnIn=2000, verbose=FALSE)
      obs  <- y_tr[fold]
      pred <- fm$yHat[fold]
      list(
        cor  = cor(obs, pred, use='complete.obs'),
        rmse = sqrt(mean((obs - pred)^2))
      )
    })
    list(
      trait    = trait,
      alpha    = a,
      mean_cor  = round(mean(sapply(fold_stats, `[[`, 'cor'),  na.rm=TRUE), 4),
      mean_rmse = round(mean(sapply(fold_stats, `[[`, 'rmse'), na.rm=TRUE), 4)
    )
  })
  names(res_a) <- traits
  res_a
}

# Run all alphas in parallel
all_cv <- mclapply(alphas, run_alpha_cv, mc.cores=11)
names(all_cv) <- paste0('a_', alphas)

# --- Build summary tables ---
cor_tab <- do.call(cbind, lapply(alphas, function(a)
  sapply(traits, function(t) all_cv[[paste0('a_',a)]][[t]]$mean_cor)
))
rownames(cor_tab) <- traits
colnames(cor_tab) <- paste0('a_', alphas)
cor_tab           <- as.data.frame(cor_tab)
cor_tab$best_alpha <- alphas[apply(cor_tab, 1, which.max)]
cor_tab$best_cor   <- apply(cor_tab[, paste0('a_', alphas)], 1, max)

rmse_tab <- do.call(cbind, lapply(alphas, function(a)
  sapply(traits, function(t) all_cv[[paste0('a_',a)]][[t]]$mean_rmse)
))
rownames(rmse_tab) <- traits
colnames(rmse_tab) <- paste0('a_', alphas)
rmse_tab           <- as.data.frame(rmse_tab)
rmse_tab$best_alpha <- alphas[apply(rmse_tab, 1, which.min)]
rmse_tab$best_rmse  <- apply(rmse_tab[, paste0('a_', alphas)], 1, min)

cat('\n--- CV Correlation (higher = better) ---\n'); print(cor_tab)
cat('\n--- CV RMSE (lower = better) ---\n');         print(rmse_tab)

write.csv(cor_tab,  'cv_cor_all_alphas.csv',  row.names=TRUE)
write.csv(rmse_tab, 'cv_rmse_all_alphas.csv', row.names=TRUE)

# --- Step 2: Predict 2025-26 using best alpha per trait (from cor) ---
dir.create('2025-26_pred', showWarnings=FALSE)
preds_best <- Res[pr, c('Season','Genotype')]

for (trait in traits) {
  best_a <- cor_tab[trait, 'best_alpha']
  cat('\nPredicting:', trait, '| best alpha:', best_a, '\n')
  
  # Full H: P = 0 for test rows/cols
  P_full           <- matrix(0, nrow=length(uid), ncol=length(uid))
  rownames(P_full) <- colnames(P_full) <- uid
  P_full[tr, tr]   <- P_std[tr, tr]
  
  H_full <- best_a * P_full + (1-best_a) * ZGZt
  eH     <- eigen(H_full, symmetric=TRUE)
  
  ETA <- list(
    E = list(V=eE$vectors, d=eE$values, model='RKHS'),
    L = list(V=eL$vectors, d=eL$values, model='RKHS'),
    H = list(V=eH$vectors, d=eH$values, model='RKHS')
  )
  
  y      <- Res[, trait]
  y[pr]  <- NA
  fm     <- BGLR(y=y, ETA=ETA, nIter=12000, burnIn=2000, verbose=FALSE)
  preds_best[[trait]] <- fm$yHat[pr]
}

write.csv(preds_best, '2025-26_pred/preds_best_alpha.csv', row.names=FALSE)
cat('\nDone.\n')