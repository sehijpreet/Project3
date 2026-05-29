rm(list=ls())
library(lme4)

getwd()
setwd('/blue/jhernandezjarqui/skaur4/P3/Pipeline/2.BLUEs/output/')

Pheno <- read.csv('../../1. Standardizing_phenomics/output/Pheno_all.csv')
Res <- read.csv('../../../Preparing_data/output/Response.csv')

summary(Pheno)
summary(Res)

Pheno[1:5, 1:8]
Pheno$Season <- as.factor(Pheno$Season)
Pheno$Date <- as.factor(Pheno$Date)
Pheno$SeasonTime <- as.factor(Pheno$SeasonTime)
Pheno$Rep <- as.factor(Pheno$Rep)
Pheno$Genotype <- as.factor(Pheno$Genotype)
Pheno$Stage<- as.factor(Pheno$Stage)
str(Pheno$Stage)

Pheno1 <- Pheno
#sum(is.na(Pheno1$Stage))

head(Res)

Res <- Res[,-c(17:24)]
head(Res)

Res$Season <- as.factor(Res$Season)
Res$Rep <- as.factor(Res$Rep)
Res$Row <- as.factor(Res$Row)
Res$Column <- as.factor(Res$Column)
Res$Bed <- as.factor(Res$Bed)
Res$Genotype <- as.factor(Res$Genotype)

Res1 <- Res

########### BLUEs for PHeno file  #####
colnames(Pheno1)

library(lme4)
library(dplyr)

head(Pheno1)
traits <- colnames(Pheno1)[7:103]

sum(is.na(Pheno1$Rep))
Pheno2 <- Pheno1

library(parallel)

run_trait <- function(tr){

  blue_list <- list()

  for(season_i in unique(Pheno2$Season)){
    for(stage_i in levels(Pheno2$Stage)){

      tmp <- subset(Pheno2, Season==season_i & Stage==stage_i)
      tmp <- tmp[!is.na(tmp[[tr]]),]

      if(nrow(tmp)>0){

        fit <- lmer(as.formula(paste0(tr," ~ Genotype + (1|Rep)")), data=tmp)

        coeffs <- fixef(fit)
        mu <- coeffs[1]
        geno_eff <- coeffs[grepl("^Genotype",names(coeffs))]
        yH <- mu + geno_eff

        out <- data.frame(
          Genotype=c(levels(tmp$Genotype)[1],gsub("Genotype","",names(yH))),
          BLUE=c(mu,yH),
          Season=season_i,
          Stage=stage_i,
          Trait=tr
        )

        blue_list[[length(blue_list)+1]] <- out
      }
    }
  }
  do.call(rbind,blue_list)
}

Pheno_BLUEs_long <- do.call(rbind,  mclapply(traits, run_trait, mc.cores=8))
write.csv(Pheno_BLUEs_long, file= 'Pheno_BLUEs_long.csv', row.names= FALSE)

# convert to wide dataframe

Pheno_BLUEs_wide <- Pheno_BLUEs_long %>%
  tidyr::pivot_wider(
    id_cols = c(Season, Stage, Genotype),
    names_from = Trait,
    values_from = BLUE)

head(Pheno_BLUEs_wide)

dim(Pheno_BLUEs_wide)

# save

write.csv(Pheno_BLUEs_long, "Phenomic_BLUEs_long.csv",row.names = FALSE)

write.csv(Pheno_BLUEs_wide,  "Phenomic_BLUEs_wide.csv",  row.names = FALSE)
