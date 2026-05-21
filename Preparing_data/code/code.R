
rm(list=ls())

getwd()

setwd('../../Preparing_data/output/')
pheno <- read.csv('../data/BVI_mean.csv')

load('../data/X4.rda')

geno <- X4

res0 <- read.csv('../data/Phenos.csv')

table(res0$Season)

res <- res0[res0$Season %in% c('2020-21', '2021-22', '2022-23', '2023-24', '2024-25', '2025-26'), ]
table(res$Season)

length(unique(res$Genotype))
length(unique(pheno$Genotype))
length(unique(rownames(geno)))

a <- Reduce(intersect, list(res$Genotype, pheno$Genotype, rownames(geno)))

pheno1 <- pheno[pheno$Genotype %in% a,]
geno1<- geno[rownames(geno) %in% a,]
res1<- res[res$Genotype %in% a, ]
dim(geno1)
g_pheno1 <- unique(pheno1$Genotype)
g_res1 <- unique(res1$Genotype)

library(tidyr)
library(dplyr)
head(res1)
res2 <- res1 %>%
  pivot_wider(
    id_cols = c(Season, Rep, Bed, Column, Row, Genotype),
    names_from = TraitDate,
    values_from = Value
  )

table(pheno1$Season)
table(res2$Season)
table(pheno1$Season, pheno1$SeasonTime)

#####time summary
time_summary <- pheno1 %>%
  group_by(Season, SeasonTime) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(Season, SeasonTime)
time_summary


library(ggplot2)

ggplot(time_summary, aes(x = SeasonTime, y = Season,size = n)) +geom_point() +theme_bw()

ref_dap <- c(45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145, 155)

stage_names <- paste0("S", 1:length(ref_dap))

ref_table <- data.frame(Stage = stage_names, RefDAP = ref_dap)
ref_table

flight_table <- unique(pheno1[, c("Season", "SeasonTime")])
flight_table <- flight_table[order(flight_table$Season, flight_table$SeasonTime),]

nearest_list <- list()

for(season_i in unique(flight_table$Season))
{
  season_times <- sort(unique(flight_table$SeasonTime[flight_table$Season == season_i]))
  for(i in seq_along(ref_dap))
  {
    target <- ref_dap[i]
    diffs <- abs(season_times - target)
    nearest <- season_times[which.min(diffs)]
    # if tie, choose earlier flight
    nearest <- min(season_times[diffs == min(diffs)])
    nearest_list[[length(nearest_list) + 1]] <- data.frame(Season = season_i, Stage = stage_names[i],
                                                           RefDAP = target, SelectedDAP = nearest)
  }
}
nearest_table <- do.call(rbind, nearest_list)
nearest_table

nearest_table[order(nearest_table$Stage, nearest_table$Season), ]



pheno2 <- merge(pheno1,
                nearest_table[, c("Season", "Stage", "SelectedDAP")],
                by.x = c("Season", "SeasonTime"),
                by.y = c("Season", "SelectedDAP"),
                all = FALSE
)
head(pheno2)


pheno2$Stage <- factor(pheno2$Stage,levels = paste0("S", 1:12))

pheno2 <- pheno2[order(pheno2$Season, pheno2$Stage),]
head(pheno2)

stage_col <- which(names(pheno2) == "Stage")

pheno2 <- pheno2[, c(1:2, stage_col,setdiff(seq_along(pheno2), c(1:4, stage_col)))]

table(pheno2$Season, pheno2$Stage)

#####time summary
time_summary <- pheno2 %>%
  group_by(Season, Stage) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(Season, Stage)
time_summary
library(ggplot2)
ggplot(time_summary, aes(x = Stage, y = Season,size = n)) +geom_point() +theme_bw()

head(pheno2)

pheno2 <- pheno2[,-7]
head(pheno2)


struct_traits <- c("CanopyArea","CanopyVolume", "Height")

plot_data <- pheno2 %>%
  select(Season, Stage, all_of(struct_traits)) %>%
  pivot_longer(
    cols = all_of(struct_traits),
    names_to = "Trait",
    values_to = "Value"
  )

ggplot(plot_data, aes(x = Stage, y = Value, color = Season, group = Season)) +
  stat_summary(fun = mean,geom = "line",linewidth = 1) +
  stat_summary(fun = mean, geom = "point") +
  facet_wrap(~Trait,scales = "free_y") +
  theme_bw()

pheno2[pheno2$Stage=='S12',]

pheno2[pheno2$Stage=='S12' & pheno2$Season=='2021-22',]

ggplot(plot_data %>%
         filter(Trait %in% c("CanopyArea", "CanopyVolume")),aes(x = Season, y = Value,fill = Season)) +
  geom_boxplot() +
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw()

ggplot(plot_data %>%
         filter(Trait %in% c("CanopyArea","CanopyVolume", "Height")),
       aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot() +
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw()

index_traits <- c("NDVI", "ARI2", "NDRE","RRI2", "TNDVI","WDRVI")

plot_indices <- pheno2 %>%
  select(Season, Stage, all_of(index_traits)) %>%
  pivot_longer(
    cols = all_of(index_traits),
    names_to = "Trait",
    values_to = "Value"
  )

ggplot(plot_indices, aes(x = Stage,  y = Value, color = Season,   group = Season)) +
  stat_summary(fun = mean,   geom = "line",       linewidth = 1) +
  stat_summary(fun = mean,   geom = "point") +
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw()

ggplot(plot_indices, aes(x = Season,   y = Value,     fill = Season)) +
  geom_boxplot() +
  facet_wrap(~Trait, scales = "free_y") +
  theme_bw()

View(ref_table)
View(nearest_table)

write.csv(pheno2, file= 'Phenomics.csv', row.names= FALSE)
write.csv(geno1, file= 'Genomics.csv', row.names= FALSE)
write.csv(nearest_table, file= 'Phenomic_Stages.csv', row.names= FALSE)
write.csv(res2, file= 'Response.csv', row.names= FALSE)
