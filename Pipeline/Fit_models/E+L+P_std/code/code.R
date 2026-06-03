####E+L+P_std
rm(list=ls())
getwd()
setwd('../output/')
library(BGLR)

source('../input/parameters.R')

# --- Load data ---
Y   <- read.csv(phenotype.file, header = TRUE, stringsAsFactors = FALSE)
uid <- paste(Y[, colSEA], Y[, colVAR], sep = "_")
head(Y)

# --- Trait columns ---
trait_cols  <- 3:ncol(Y)
trait_names <- colnames(Y)[trait_cols]

# --- Prediction set: 2025-26 ---
pred_idx  <- which(Y[, colSEA] == "2025-26")
train_idx <- which(Y[, colSEA] != "2025-26")


# --- Build ETA from AB list ---
nk  <- length(AB)
ETA <- vector("list", nk)


for(i in 1:nk){
  if(type[i] == 'RKHS'){
    tmp <- load(AB[[i]])   # returns object name as string
    EVD <- get(tmp)        # get object regardless of name
    ETA[[i]] <- list(V = EVD$vectors, d = EVD$values, model = 'RKHS')
    rm(list = tmp)         # removes EVD_E, EVD_L, EVD etc.
  }
}

# Clean up EVD if it lingered
if(exists("EVD")) rm(EVD)


# --- Output folder for 2025-26 predictions ---
dir.create("2025-26_pred", showWarnings = FALSE)

summary_list <- vector("list", length(trait_names))
names(summary_list) <- trait_names


for(trait in trait_names){
  
  # Per trait folder
  dir.create(trait, showWarnings = FALSE)
  
  y <- Y[, trait]
  if(ESC){ y <- scale(y, center = TRUE, scale = TRUE) }
  
  # Mask 2025-26
  yNA           <- y
  yNA[pred_idx] <- NA
  
  # Fit model
  fm <- BGLR( y = yNA,ETA = ETA, nIter= nIter, burnIn = burnIn, saveAt = paste0(trait, "/"),verbose = FALSE  )
  fm$y <- y
  
  # Training accuracy
  acc <- cor(y[train_idx], fm$yHat[train_idx], use = "complete.obs")
  
  # All predictions
  all_preds <- data.frame(Season= Y[, colSEA], Genotype = Y[, colVAR],y_obs= y,y_hat= fm$yHat )
  rownames(all_preds) <- uid
  
  write.csv(all_preds, paste0(trait, "/predictions_", trait, ".csv"), row.names = TRUE)
  
  # 2025-26 predictions only
  pred_out <- all_preds[pred_idx, ]
  write.csv(pred_out,  paste0("2025-26_pred/", trait, "_2025-26.csv"), row.names = TRUE)
  
  # Summary
  summary_list[[trait]] <- data.frame(Trait  = trait,
    Cor    = round(acc, 3),
    VarE   = round(fm$varE, 4),
    VarK1  = round(fm$ETA[[1]]$varU, 4),
    VarK2  = round(fm$ETA[[2]]$varU, 4),
    VarK3  = round(fm$ETA[[3]]$varU, 4)
  )
  
  # Clean BGLR temp files
  unlink(paste0(trait, "/*.dat"))
  
  print(str(fm))
  rm(fm)
}


