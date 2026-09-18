rm(list=ls())

getwd()
setwd('../Desktop/PhD data/P3/Pipeline/3.Phenomic_profiling/output/Final_outs/')
setwd('../../../Fit_models/MegaLMM/output/')

#X <- read.csv('../../3.Phenomic_profiling/output/Final_outs/Geno_matched.csv')
Pheno_orig <- read.csv('../../../3.Phenomic_profiling/output/Final_outs1/Pheno_orig_matched.csv')
Pheno_std <- read.csv('../../../3.Phenomic_profiling/output/Final_outs1/Pheno_std_matched.csv')
Res <- read.csv('../../../3.Phenomic_profiling/output/Final_outs1/Res1_matched.csv')

X <- read.csv('../../../3.Phenomic_profiling/output/Final_outs1/Geno_matched.csv')
X[1:5, 1:5]

Pheno_std[1:8, 1:8]
Res[1:5, 1:5]

Y <- merge( Pheno_std, Res,  by = c("Season", "Genotype"),  all = FALSE)
dim(Y)
head(Y)

load('../../../4.G/output/ZGZt.rda')
dim(ZGZt)
G<- ZGZt

Y2 <- Y[, !grepl("^X[0-9]+$", colnames(Y))]
table(Y2$Season)

Y2_tr <- Y2
Y2_tr[Y2_tr$Season=='2020-21',-(1:2)] <- NA

G[1:5, 1:5]

id <- paste(Y2$Season, Y2$Genotype, sep="_")
G <- G[id, id]
Y_tr <- as.matrix(Y2_tr[, -(1:2)])
identical(id, rownames(G))
dim(Y_tr)
dim(G)

if(!require(devtools)) { install.packages("devtools"); library(devtools) }
if(!require(MegaLMM)) { 
  devtools::install_github('deruncie/MegaLMM')
  library(MegaLMM) 
}


sample_data <- data.frame(ID=factor(id, levels=id))
run_parameters <- MegaLMM_control(h2_divisions=10, burn=0, thin=2, K=15)


MegaLMM_state <- setup_model_MegaLMM(Y_tr, ~1+(1|ID), data=sample_data,
                                     relmat=list(ID=G), run_parameters=run_parameters, run_ID="MegaLMM_2020_21")


##setting priors

Lambda_prior <- list(sampler=sample_Lambda_prec_ARD, Lambda_df=3, delta_1=list(shape=2,rate=1), delta_2=list(shape=3,rate=1), delta_iterations_factor=100)

priors <- MegaLMM_priors(tot_Y_var=list(V=.5,nu=5), tot_F_var=list(V=.9,nu=20), h2_priors_resids_fun=function(h2s,n) 1, h2_priors_factors_fun=function(h2s,n) 1, Lambda_prior=Lambda_prior)

MegaLMM_state <- set_priors_MegaLMM(MegaLMM_state, priors)


####

maps <- make_Missing_data_map(MegaLMM_state, max_NA_groups=3, verbose=FALSE)
MegaLMM_state <- set_Missing_data_map(MegaLMM_state, maps$Missing_data_map)
MegaLMM_state <- initialize_variables_MegaLMM(MegaLMM_state)

estimate_memory_initialization_MegaLMM(MegaLMM_state)


####
MegaLMM_state <- initialize_MegaLMM(MegaLMM_state, verbose=TRUE)
MegaLMM_state$Posterior$posteriorSample_params <- c("Lambda","F_h2","resid_h2","tot_Eta_prec","B1")
MegaLMM_state$Posterior$posteriorMean_params <- "Eta_mean"
MegaLMM_state$Posterior$posteriorFunctions <- list(U="U_F %*% Lambda + U_R + X1 %*% B1")


MegaLMM_state <- clear_Posterior(MegaLMM_state)
estimate_memory_posterior(MegaLMM_state, 100)


for(i in 1:5){
  MegaLMM_state <- reorder_factors(MegaLMM_state, drop_cor_threshold=.6)
  MegaLMM_state <- clear_Posterior(MegaLMM_state)
  MegaLMM_state <- sample_MegaLMM(MegaLMM_state,100)
  print(paste("Burn-in round",i,"completed"))
}

MegaLMM_state <- clear_Posterior(MegaLMM_state)


for(i in 1:4){
  MegaLMM_state <- sample_MegaLMM(MegaLMM_state,250)
  MegaLMM_state <- save_posterior_chunk(MegaLMM_state)
  print(paste("Sampling round",i,"completed"))
}

MegaLMM_state$Posterior <- reload_Posterior(MegaLMM_state)
U_hat <- get_posterior_mean(MegaLMM_state,U)


U_samples <- load_posterior_param(MegaLMM_state,"U")
U_hat <- get_posterior_mean(U_samples)
dim(U_hat)

test <- Y2$Season=="2020-21"
cum <- grep("^X[0-9]+_cum$", colnames(Y_tr))

pred <- U_hat[test,cum,drop=FALSE]
obs <- as.matrix(Y2[test,-(1:2)])[ ,cum,drop=FALSE]

accuracy <- diag(cor(obs,pred,use="pairwise.complete.obs"))
data.frame(Trait=colnames(Y_tr)[cum], Accuracy=accuracy)

saveRDS(MegaLMM_state,"MegaLMM_state_2020-21.rds")

pred_out <- cbind(Y2[test,c("Season","Genotype")], as.data.frame(pred))
write.csv(pred_out,"MegaLMM_predictions_2020-21.csv",row.names=FALSE)

pred_all <- cbind(Y2[test,c("Season","Genotype")],as.data.frame(U_hat[test,]))
write.csv(pred_all,"MegaLMM_all_predictions_2020-21.csv",row.names=FALSE)

accuracy_out <- data.frame(Trait=colnames(Y_tr)[cum],Correlation=accuracy)
write.csv(accuracy_out,"MegaLMM_accuracy_2020-21.csv",row.names=FALSE)


###Check convergence
Lambda_samples <- load_posterior_param(MegaLMM_state,"Lambda")
traceplot_array(Lambda_samples,facet_dim=2,name="Lambda_trace")
traceplot_array(U_samples,facet_dim=3,name="U_trace")

print(MegaLMM_state)
summary(MegaLMM_state)
