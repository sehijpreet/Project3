rm(list=ls())
library(MegaLMM)
getwd()

setwd("../../Fit_models/MegaLMM/output2")

# Load data
Pheno_std <- read.csv("../../../3.Phenomic_profiling/output/Final_outs1/Pheno_std_matched.csv")
Res <- read.csv("../../../3.Phenomic_profiling/output/Final_outs1/Res1_matched.csv")
load("../../4.G/output/ZGZt.rda")
G <- ZGZt

# Merge and select only cumulative yield
Y <- merge(Pheno_std,Res,by=c("Season","Genotype"))
cum_names <- grep("^X[0-9]+_cum$",names(Y),value=TRUE)
Y_cum <- Y[,c("Season","Genotype",cum_names)]

# Preserve observed data and hide 2025-26
Y_cum_obs <- Y_cum
Y_cum_tr <- Y_cum
Y_cum_tr[Y_cum_tr$Season=="2025-26",-(1:2)] <- NA

# Align genomic matrix
id <- paste(Y_cum$Season,Y_cum$Genotype,sep="_")
stopifnot(all(id %in% rownames(G)))
G <- G[id,id]
stopifnot(identical(id,rownames(G)))

# MegaLMM inputs
Y_tr <- as.matrix(Y_cum_tr[,-(1:2)])
sample_data <- data.frame(ID=factor(id,levels=id))

# Model settings
run_parameters <- MegaLMM_control(h2_divisions=10,
  burn=0,
  thin=2,
  K=8,
  scale_Y=TRUE
)

# Construct model
MegaLMM_state <- setup_model_MegaLMM(
  Y_tr,~1+(1|ID),
  data=sample_data,
  relmat=list(ID=G),
  run_parameters=run_parameters,
  run_ID="MegaLMM_yield_only_2025_26"
)

# Priors
Lambda_prior <- list(
  sampler=sample_Lambda_prec_ARD,
  Lambda_df=3,
  delta_1=list(shape=2,rate=1),
  delta_2=list(shape=3,rate=1),
  delta_iterations_factor=100
)

priors <- MegaLMM_priors(
  tot_Y_var=list(V=.5,nu=5),
  tot_F_var=list(V=.9,nu=20),
  h2_priors_resids_fun=function(h2s,n) 1,
  h2_priors_factors_fun=function(h2s,n) 1,
  Lambda_prior=Lambda_prior
)

MegaLMM_state <- set_priors_MegaLMM(MegaLMM_state,priors)

# Missing-data structure
maps <- make_Missing_data_map(
  MegaLMM_state,
  max_NA_groups=3,
  verbose=FALSE
)

MegaLMM_state <- set_Missing_data_map(
  MegaLMM_state,
  maps$Missing_data_map
)

# Initialize
MegaLMM_state <- initialize_variables_MegaLMM(MegaLMM_state)
estimate_memory_initialization_MegaLMM(MegaLMM_state)
MegaLMM_state <- initialize_MegaLMM(MegaLMM_state,verbose=TRUE)

# Results to save
MegaLMM_state$Posterior$posteriorSample_params <-
  c("Lambda","F_h2","resid_h2","tot_Eta_prec","B1")

MegaLMM_state$Posterior$posteriorMean_params <- "Eta_mean"

MegaLMM_state$Posterior$posteriorFunctions <- list(
  U="U_F %*% Lambda + U_R + X1 %*% B1"
)

MegaLMM_state <- clear_Posterior(MegaLMM_state)

# Burn-in
for(i in 1:5){
  MegaLMM_state <- reorder_factors(
    MegaLMM_state,
    drop_cor_threshold=.6
  )
  MegaLMM_state <- clear_Posterior(MegaLMM_state)
  MegaLMM_state <- sample_MegaLMM(MegaLMM_state,100)
  print(paste("Burn-in round",i,"completed"))
}

MegaLMM_state <- clear_Posterior(MegaLMM_state)

# Final posterior sampling
for(i in 1:4){
  MegaLMM_state <- sample_MegaLMM(MegaLMM_state,250)
  MegaLMM_state <- save_posterior_chunk(MegaLMM_state)
  print(paste("Sampling round",i,"completed"))
}

# Extract predictions
U_samples <- load_posterior_param(MegaLMM_state,"U")
U_hat <- get_posterior_mean(U_samples)

# Evaluate 2025-26
test <- Y_cum$Season=="2025-26"
pred <- U_hat[test,,drop=FALSE]
obs <- as.matrix(Y_cum_obs[test,-(1:2)])

accuracy <- diag(cor(
  obs,pred,
  use="pairwise.complete.obs"
))

# Save predictions
pred_out <- cbind(
  Y_cum[test,c("Season","Genotype")],
  as.data.frame(pred)
)

write.csv(
  pred_out,
  "MegaLMM_yield_only_predictions_2025_26.csv",
  row.names=FALSE
)

# Save correlations
accuracy_out <- data.frame(
  Trait=colnames(Y_tr),
  Yield_only_correlation=accuracy
)

write.csv(
  accuracy_out,
  "MegaLMM_yield_only_accuracy_2025_26.csv",
  row.names=FALSE
)

# Save model
saveRDS(
  MegaLMM_state,
  "MegaLMM_yield_only_state_2025_26.rds"
)
