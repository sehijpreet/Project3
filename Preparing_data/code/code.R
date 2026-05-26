
rm(list=ls())

getwd()

setwd('../../../Preparing_data/output/')
pheno <- read.csv('../data/BVI_mean.csv')
pheno0 <- read.csv('../../Exploratory_DA/output/BVI.csv')

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
                                                           RefDAP = target, SeasonTime = nearest)
  }
}
nearest_table <- do.call(rbind, nearest_list)
nearest_table

nearest_table[order(nearest_table$Stage, nearest_table$Season), ]

str(pheno1$SeasonTime)
str(nearest_table$SeasonTime)

nearest_table[nearest_table$Season=='2024-25',]
table(pheno1$SeasonTime[pheno1$Season=='2024-25'])

nearest_table$SeasonTime[nearest_table$Season == "2024-25" & nearest_table$Stage == "S2"] <- 50
nearest_table$SeasonTime[nearest_table$Season == "2024-25" & nearest_table$Stage == "S3"] <- 55
nearest_table$SeasonTime[nearest_table$Season == "2024-25" & nearest_table$Stage == "S4"] <- 55

pheno2 <- merge(
  pheno1,
  nearest_table[, c("Season", "SeasonTime", "Stage")],
  by = c("Season", "SeasonTime"),
  all.x = TRUE
)
head(pheno2)
table(pheno2$Season, pheno2$Stage)
table(pheno2$SeasonTime)


length(unique(pheno$Genotype))
length(unique(pheno1$Genotype))
length(unique(pheno1$Genotype))

pheno3<- pheno2
pheno3$Stage <- factor(pheno3$Stage,levels = paste0("S", 1:12))

pheno3 <- pheno3[order(pheno3$Season, pheno3$Stage),]
head(pheno3)

stage_col <- which(names(pheno3) == "Stage")

pheno3 <- pheno3[, c(1:2, stage_col,setdiff(seq_along(pheno3), c(1:2, stage_col)))]
head(pheno3)
table(pheno3$Season, pheno3$Stage)

pheno3 <- pheno3[!is.na(pheno3$Stage), ]

#####time summary
time_summary <- pheno3 %>%
  group_by(Season, Stage) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(Season, Stage)

time_summary[1:14,]

library(ggplot2)
ggplot(time_summary, aes(x = Stage, y = Season, size = n)) +geom_point() +theme_bw()

head(pheno3)

pheno3 <- pheno3[,-10]
head(pheno3)


##############################################
pheno <- pheno[
  order(pheno$Season,
        pheno$Genotype,
        pheno$SeasonTime),
]

pheno1 <- pheno1[
  order(pheno1$Season,
        pheno1$Genotype,
        pheno1$SeasonTime),
]
pheno2 <- pheno2[
  order(pheno2$Season,
        pheno2$Genotype,
        pheno2$SeasonTime),
]
pheno3 <- pheno3[
  order(pheno3$Season,
        pheno3$Genotype,
        pheno3$SeasonTime),
]
pheno[1:5, 1:5]

pheno1[1:5, 1:5]

tiff("NA_patterns.tif",width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))

image(t(is.na(pheno3)),main = "Pheno3")
image(t(is.na(pheno2)),main = "Pheno2")
image(t(is.na(pheno1)),main = "Pheno1")
image(t(is.na(pheno)), main = "Pheno")

dev.off()

#image(t(is.na(pheno0)),main = "Pheno0")
na_before <- colMeans(is.na(pheno2))
na_after <- colMeans(is.na(pheno3))


head(na_compare)

dim(pheno1)
dim(pheno3)

############################################################


struct_traits <- c("CanopyArea","CanopyVolume", "Height")

plot_data <- pheno3 %>%
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
index_traits <- c("R", "G", "B","NIR", "EXG","RRI2")

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

head(pheno3)
write.csv(pheno3, file= 'Phenomics.csv', row.names= FALSE)
write.csv(geno1, file= 'Genomics.csv', row.names= FALSE)
write.csv(nearest_table, file= 'Phenomic_Stages1.csv', row.names= FALSE)
write.csv(res2, file= 'Response.csv', row.names= FALSE)
