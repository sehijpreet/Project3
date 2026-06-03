rm(list=ls())

getwd()
setwd('../../../OUTS/E+L+G+P_std/output/')

Res <- read.csv('../../../3.Phenomic_profiling/output/Final_outs/Res1_matched.csv')

trait_names <- colnames(Res)[3:ncol(Res)]

# --- Load all predictions in a loop ---
pred_path <- '../../../Fit_models/E+L+G+P_std/output/2025-26_pred/'

Pred_2526 <- Res[Res$Season == "2025-26", c("Season", "Genotype")]
rownames(Pred_2526) <- paste(Pred_2526$Season, Pred_2526$Genotype, sep = "_")

for(trait in trait_names){
  tmp <- read.csv(paste0(pred_path, trait, "_2025-26.csv"))
  
  key_pred <- paste(tmp$Season, tmp$Genotype, sep = "_")
  key_out  <- rownames(Pred_2526)
  
  idx <- match(key_out, key_pred)
  Pred_2526[, paste0(trait, "_hat")] <- tmp$y_hat[idx]
}


Obs_2526 <- Res[Res$Season == "2025-26", ]
rownames(Obs_2526) <- paste(Obs_2526$Season, Obs_2526$Genotype, sep = "_")

# Combine obs + pred side by side
Out <- merge(Obs_2526, Pred_2526[, !colnames(Pred_2526) %in% c("Season", "Genotype")], by = "row.names")

rownames(Out) <- Out$Row.names
Out$Row.names <- NULL

obs_cols  <- trait_names
pred_cols <- paste0(trait_names, "_hat")
col_order <- c("Season", "Genotype", as.vector(rbind(obs_cols, pred_cols)))

Out <- Out[, col_order]

Out <- Out[order(Out$Genotype), ]

# --- Verify ---
dim(Out)
head(Out[, 1:8])

write.csv(Out, "ELGP_std_2025-26_obs_vs_pred.csv", row.names = FALSE)



cor_df <- data.frame( Trait   = trait_names, 
                      Cor     = sapply(trait_names, function(t){
                        round(cor(Out[, t], Out[, paste0(t, "_hat")], use = "complete.obs"), 3)
                      })
)

print(cor_df)
write.csv(cor_df, "ELGP_std_2025-26_correlations.csv", row.names = FALSE)
