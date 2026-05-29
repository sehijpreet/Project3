rm(list=ls())
library(lme4)

getwd()
setwd('../../2.BLUEs/output/')

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


blue_list <- list()

for(tr in traits) {
  
  for(season_i in unique(Pheno2$Season)) {
    
    for(stage_i in levels(Pheno2$Stage)) {
      
      tmp <- Pheno2[Pheno2$Season == season_i &  Pheno2$Stage == stage_i,]
      
      tmp <- tmp[!is.na(tmp[[tr]]),]
      
      if(nrow(tmp) > 0) {
        
      # mixed model
      fit <- lmer(as.formula( paste0(tr," ~ Genotype + (1|Rep)" )),data = tmp)
        
      # fixed effects
      coeffs <- fixef(fit)
        
      mu <- coeffs[1]
        
      geno_eff <- coeffs[grepl("^Genotype", names(coeffs))]
        
      yH <- mu + geno_eff
        
      # genotype names
      nams <- gsub( "Genotype", "",  names(yH))
        
      names(yH) <- nams
        
    # reference genotype
      ref_geno <- levels(as.factor(tmp$Genotype))[1]
        
        ref_row <- data.frame(Genotype = ref_geno,BLUE = mu )
        
        other_rows <- data.frame( Genotype = nams,  BLUE = yH )
        
        out <- rbind( ref_row,  other_rows)
        out$Season <- season_i
        out$Stage  <- stage_i
        out$Trait  <- tr
        
        blue_list[[length(blue_list)+1]] <- out
        
        print(c(tr, season_i, stage_i))
      }
    }
  }
}

# combine all BLUEs
Pheno_BLUEs_long <- do.call(rbind,  blue_list)

head(Pheno_BLUEs_long)

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


##### BLUEs for Response #######


library(lme4)
library(dplyr)
library(tidyr)

head(Res)

# response traits
traits <- colnames(Res)[7:ncol(Res)]

Res2 <- Res[
  !is.na(Res$Rep) &
    !is.na(Res$Row) &
    !is.na(Res$Column),
]

blue_list <- list()

for(tr in traits) {
  
  for(season_i in unique(Res2$Season)) {
    
    for(stage_i in unique(Res2$Stage)) {
      
      tmp <- Res2[Res2$Season == season_i &  Res2$Stage == stage_i, ]
      
      tmp <- tmp[!is.na(tmp[[tr]]),]
      
      if(nrow(tmp) > 0) {
        
        # mixed model
        fit <- lmer(as.formula( paste0( tr, " ~ Genotype + ","(1|Rep) + ", "(1|Bed) + ", "(1|Row) + ","(1|Column)")), data = tmp)        
        
        # fixed effects
        coeffs <- fixef(fit)
        
        mu <- coeffs[1]
        
        geno_eff <- coeffs[ grepl("^Genotype", names(coeffs))]
        
        yH <- mu + geno_eff
        
        # genotype names
        nams <- gsub("Genotype","", names(yH))
        
        names(yH) <- nams
        
        # reference genotype
        ref_geno <- levels(as.factor(tmp$Genotype))[1]
        
        ref_row <- data.frame(Genotype = ref_geno,BLUE = mu)
        
        other_rows <- data.frame(Genotype = nams, LUE = yH)
        
        out <- rbind( ref_row, other_rows  )
        
        out$Season <- season_i
        out$Stage  <- stage_i
        out$Trait  <- tr
        
        blue_list[[length(blue_list)+1]] <- out
        
        print(c(tr, season_i, stage_i))
      }
    }
  }
}

# combine all BLUEs
Res_BLUEs_long <- do.call(rbind,blue_list)

head(Res_BLUEs_long)

# convert to wide dataframe
Res_BLUEs_wide <- Res_BLUEs_long %>%
  pivot_wider(id_cols = c(Season,  Stage,  Genotype  ),
    names_from = Trait,
    values_from = BLUE
  )

head(Res_BLUEs_wide)

dim(Res_BLUEs_wide)

# save
write.csv(  Res_BLUEs_long,  "Response_BLUEs_long.csv",  row.names = FALSE)

write.csv(  Res_BLUEs_wide,  "Response_BLUEs_wide.csv",  row.names = FALSE)

