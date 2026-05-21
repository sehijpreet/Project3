rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)

getwd()

setwd('../../Pipeline/1. Standardizing_phenomics/output/')

Res <- read.csv('../../../Preparing_data/output/Response.csv')
Pheno<- read.csv('../../../Preparing_data/output/Phenomics.csv')
Geno <- read.csv('../../../Preparing_data/output/Genomics.csv')

head(Pheno)
head(Geno)
head(Res)

table(Pheno$SeasonTime)
table(Pheno$Stage)
199646/12

# structural traits to standardize
struct_traits <- c("CanopyArea",  "CanopyVolume",  "Height")

# BEFORE standardization ----------------------------

Pheno$Stage <- factor(Pheno$Stage,levels = paste0("S", 1:12))

plot_before <- Pheno %>%
  select(Season, Stage, all_of(struct_traits)) %>%
  pivot_longer(
    cols = all_of(struct_traits),
    names_to = "Trait",
    values_to = "Value"
  )

ggplot(plot_before,  aes(x = Stage, y = Value,color = Season, group = Season)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point") +
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw()


############  STANDARDIZE WITHIN SEASON X STAGE ####################

St_Pheno <- Pheno %>%
  group_by(Season) %>%
  mutate(
    across(
      all_of(struct_traits),
      ~ as.numeric(scale(.)),
      .names = "{.col}_z"
    )
  ) %>%
  ungroup()





scaled_traits <- c( "CanopyArea_z",  "CanopyVolume_z",  "Height_z")

plot_after <- St_Pheno %>%
  select(Season, Stage, all_of(scaled_traits)) %>%
  pivot_longer(
    cols = all_of(scaled_traits),
    names_to = "Trait",
    values_to = "Value"
  )

ggplot(plot_after, aes(x = Stage,   y = Value,  color = Season, group = Season)) +
  stat_summary(fun = mean, geom = "line",  linewidth = 1) +
  stat_summary(fun = mean, geom = "point") +
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw()


################################# Check indices ##############

head(Pheno)



# structural traits to standardize
struct_traits <- c("GDVI",  "GEMI",  "G", "RRI2", "NDVI", "R", "RE", "NIR", "DVI")

# BEFORE standardization ----------------------------

Pheno$Stage <- factor(Pheno$Stage,levels = paste0("S", 1:12))

plot_before <- Pheno %>%
  select(Season, Stage, all_of(struct_traits)) %>%
  pivot_longer(
    cols = all_of(struct_traits),
    names_to = "Trait",
    values_to = "Value"
  )

ggplot(plot_before,  aes(x = Stage, y = Value,color = Season, group = Season)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point") +
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw()



summary(St_Pheno)






