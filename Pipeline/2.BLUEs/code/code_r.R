rm(list=ls())
library(lme4)

getwd()
setwd('../output/')

Pheno <- read.csv('../../1. Standardizing_phenomics/output/Pheno_all.csv')
Res <- read.csv('../../../Preparing_data/output/Response.csv')

summary(Pheno)
summary(Res)

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

##### BLUEs for Response #######


library(lme4)
library(dplyr)
library(tidyr)

head(Res)

# response traits
traits <- colnames(Res)[7:ncol(Res)]

Res2 <- Res[  !is.na(Res$Rep) &    !is.na(Res$Row) &    !is.na(Res$Column),]


library(parallel)

run_trait <- function(tr){

  blue_list <- list()

  for(season_i in unique(Res2$Season)) {

      tmp <- Res2[Res2$Season == season_i,]

      tmp <- tmp[!is.na(tmp[[tr]]),]

      if(nrow(tmp) > 0) {

        fit <- lmer( as.formula( paste0(tr," ~ Genotype + (1|Rep) + (1|Bed) + (1|Row) + (1|Column)") ),   data = tmp)

        coeffs <- fixef(fit)

        mu <- coeffs[1]

        geno_eff <- coeffs[
          grepl("^Genotype", names(coeffs))
        ]

        yH <- mu + geno_eff

        nams <- gsub("Genotype","", names(yH))

        ref_geno <- levels(as.factor(tmp$Genotype))[1]

        ref_row <- data.frame( Genotype = ref_geno,   BLUE = mu )

        other_rows <- data.frame(Genotype = nams, BLUE = yH )

        out <- rbind(ref_row, other_rows)
        out$Season <- season_i
        out$Trait  <- tr
        blue_list[[length(blue_list)+1]] <- out
      
    }
  }

  do.call(rbind, blue_list)
}

res <- mclapply( traits,  run_trait,  mc.cores = 8)

Res_BLUEs_long <- do.call(rbind, res)
write.csv(  Res_BLUEs_long,  "Response_BLUEs_long.csv",  row.names = FALSE)

# convert to wide dataframe
Res_BLUEs_wide <- Res_BLUEs_long %>%
  pivot_wider(id_cols = c(Season, Genotype),
    names_from = Trait,
    values_from = BLUE
  )

head(Res_BLUEs_wide)

dim(Res_BLUEs_wide)

# save

write.csv(  Res_BLUEs_wide,  "Response_BLUEs_wide.csv",  row.names = FALSE)

