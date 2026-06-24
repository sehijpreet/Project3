rm(list=ls())
getwd()
setwd('../../../OUTS/ssBLUP/output/')

Res <- read.csv('../../../3.Phenomic_profiling/output/Final_outs/Res1_matched.csv')
trait_names <- colnames(Res)[3:ncol(Res)]
pred_path <- '../../../Fit_models/ssBLUP/output/2025-26_pred/'

Pred_2526 <- Res[Res$Season == "2025-26", c("Season", "Genotype")]
rownames(Pred_2526) <- paste(Pred_2526$Season, Pred_2526$Genotype, sep = "_")

tmp <- read.csv(paste0(pred_path, "preds_best_alpha.csv"))

key_pred <- paste(tmp$Season, tmp$Genotype, sep = "_")
key_out  <- rownames(Pred_2526)
idx <- match(key_out, key_pred)

pred_cols <- c(paste0("X", 1:10), paste0("X", 1:10, "_cum"))
Pred_2526[, paste0(pred_cols, "_hat")] <- tmp[idx, pred_cols]

Obs_2526 <- Res[Res$Season == "2025-26", ]
rownames(Obs_2526) <- paste(Obs_2526$Season, Obs_2526$Genotype, sep = "_")

# Combine obs + pred side by side (rownames already aligned, no merge needed)
Out <- cbind(Obs_2526, Pred_2526[rownames(Obs_2526), !colnames(Pred_2526) %in% c("Season", "Genotype")])

obs_cols  <- trait_names
hat_cols  <- paste0(pred_cols, "_hat")
col_order <- c("Season", "Genotype", as.vector(rbind(obs_cols, hat_cols)))
Out <- Out[, col_order]
Out <- Out[order(Out$Genotype), ]

# --- Verify ---
dim(Out)
head(Out[, 1:8])

write.csv(Out, "ssBLUP_2025-26_obs_vs_pred.csv", row.names = FALSE)



cor_df <- data.frame( Trait   = trait_names, 
                      Cor     = sapply(trait_names, function(t){
                        round(cor(Out[, t], Out[, paste0(t, "_hat")], use = "complete.obs"), 3)
                      })
)

print(cor_df)
write.csv(cor_df, "ssBLUP_2025-26_correlations.csv", row.names = FALSE)
