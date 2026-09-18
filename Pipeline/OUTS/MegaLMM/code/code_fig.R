rm(list=ls())

getwd()

setwd('./output')

M2021 <- read.csv('MegaLMM_accuracy_2020-21.csv')
M2122 <- read.csv('MegaLMM_accuracy_2021-22.csv')
M2223 <- read.csv('MegaLMM_accuracy_2022-23.csv')
M2324 <- read.csv('MegaLMM_accuracy_2023-24.csv')
M2425 <- read.csv('MegaLMM_accuracy_2024_25.csv')
M2526 <- read.csv('MegaLMM_accuracy_2025_26.csv')


head(M2526)
setwd('../../../OUTS/E+L+G/output')


G2021 <- read.csv('ELG_2020-21_correlations.csv')
G2122 <- read.csv('ELG_2021-22_correlations.csv')
G2223 <- read.csv('ELG_2022-23_correlations.csv')
G2324 <- read.csv('ELG_2023-24_correlations.csv')
G2425 <- read.csv('ELG_2024-25_correlations.csv')
G2526 <- read.csv('ELG_2025-26_correlations.csv')

G2526[, 1:2]
M2526[, 1:2]

library(tidyverse)

M <- bind_rows("2020-21"=M2021,"2021-22"=M2122,"2022-23"=M2223,"2023-24"=M2324,"2024-25"=M2425,"2025-26"=M2526,.id="Year") %>% transmute(Year,Trait,Accuracy=Correlation,Model="MegaLMM")
G <- bind_rows("2020-21"=G2021,"2021-22"=G2122,"2022-23"=G2223,"2023-24"=G2324,"2024-25"=G2425,"2025-26"=G2526,.id="Year") %>% transmute(Year,Trait,Accuracy=Cor,Model="E+L+G")
df <- bind_rows(M,G) %>% filter(Trait %in% paste0("X",1:10,"_cum"))

ggplot(df, aes(Trait, Accuracy, color = Model, group = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~Year) +
  scale_x_discrete(
    limits = paste0("X", 1:10, "_cum"),
    labels = paste0("W", 1:10)
  ) +
  labs(x = "Week", y = "Predictive ability", color = "Model") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.91, 0.09),
    legend.background = element_rect(fill = "white", color = "black")
  )
