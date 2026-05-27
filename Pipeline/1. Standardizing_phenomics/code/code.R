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
table(Pheno$Season, Pheno$Stage)

id_cols <- c("Season",  "Date",  "SeasonTime",  "Stage",  "Rep",  "Genotype")
sum(duplicated(Pheno[,id_cols]))


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

tiff("st_plot_before.tif",width = 12, height = 6, units = "in", res = 300)

ggplot(plot_before,  aes(x = Stage, y = Value,color = Season, group = Season)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point") +
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw()

dev.off()

Pheno[Pheno$Season=='2021-22' & Pheno$Stage=='S11',]
############  STANDARDIZE WITHIN SEASON  ####################

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

################ PLOT AFTER#############

scaled_traits <- c( "CanopyArea_z",  "CanopyVolume_z",  "Height_z")

plot_after <- St_Pheno %>%
  select(Season, Stage, all_of(scaled_traits)) %>%
  pivot_longer(
    cols = all_of(scaled_traits),
    names_to = "Trait",
    values_to = "Value"
  )

tiff("st_plot_after.tif",width = 12, height = 6, units = "in", res = 300)

ggplot(plot_after, aes(x = Stage,   y = Value,  color = Season, group = Season)) +
  stat_summary(fun = mean, geom = "line",  linewidth = 1) +
  stat_summary(fun = mean, geom = "point") +
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw()

dev.off()
################################# Check indices ##############

head(Pheno)



# indices to standardize
colnames(Pheno)
ind_traits <- names(Pheno[11:55])

# BEFORE standardization ----------------------------

Pheno$Stage <- factor(Pheno$Stage,levels = paste0("S", 1:12))

plot_before <- Pheno %>%
  select(Season, Stage, all_of(ind_traits)) %>%
  pivot_longer(
    cols = all_of(ind_traits),
    names_to = "Trait",
    values_to = "Value"
  )

tiff("ind_plot_before.tif",width = 12, height = 6, units = "in", res = 300)
par(mfrow=c(7,7))

ggplot(plot_before,  aes(x = Stage, y = Value,color = Season, group = Season)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point") +
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw()

dev.off()
####### STANDARDISE ###########

Ind_Pheno <- Pheno %>%
  group_by(Season) %>%
  mutate(
    across(
      all_of(ind_traits),
      ~ as.numeric(scale(.)),
      .names = "{.col}_f"
    )
  ) %>%
  ungroup()

ind_traits <- names(Ind_Pheno[56:100])


plot_after <- Ind_Pheno %>%
  select(Season, Stage, all_of(ind_traits)) %>%
  pivot_longer(
    cols = all_of(ind_traits),
    names_to = "Trait",
    values_to = "Value"
  )

tiff("ind_plot_after.tif",width = 12, height = 6, units = "in", res = 300)
par(mfrow=c(7,7))

ggplot(plot_after, aes(x = Stage,   y = Value,  color = Season, group = Season)) +
  stat_summary(fun = mean, geom = "line",  linewidth = 1) +
  stat_summary(fun = mean, geom = "point") +
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw()

dev.off()

############ Save after Standardize ###########
id_cols <- c("Season",  "Date",  "SeasonTime",  "Stage",  "Rep",  "Genotype")

Pheno_st1 <- merge(St_Pheno[, c(id_cols, scaled_traits)], Ind_Pheno[, c(id_cols, ind_traits)],  by = id_cols,  all = TRUE) 
Pheno_all1 <- merge(Pheno, Pheno_st1,  by = id_cols,  all.x = TRUE)
colnames(Pheno_all1)

sum(duplicated(Pheno_st1[,id_cols]))
sum(duplicated(Pheno_all1[,id_cols]))


write.csv(Pheno_st1, file= 'Pheno_st.csv', row.names= FALSE)
write.csv(Pheno, file = 'Pheno.csv', row.names=FALSE)                  
write.csv(Pheno_all1, file = 'Pheno_all.csv', row.names=FALSE)                  



#####VAlidate if correctly merged


Pheno_all[Pheno_all$Season == "2022-23" &  Pheno_all$Genotype == "Radiance" &  Pheno_all$Stage == "S5",]

Pheno[Pheno$Season == "2022-23" &  Pheno$Genotype == "Radiance" &  Pheno$Stage == "S5",]

Pheno_st1[Pheno_st1$Season == "2022-23" &  Pheno_st1$Genotype == "Radiance" &  Pheno_st1$Stage == "S5",]


id_col <- c('Season', 'Stage', 'Rep', 'Genotype')
sum(duplicated(Pheno_all1[, id_col]))




