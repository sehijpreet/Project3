rm(list=ls())

setwd( "C:/Users/sidhu/OneDrive - University of Florida/Desktop/PhD data/P3/Pipeline/Fit_models/P_pred/Within_S/2025-26/output/")

df <- read.csv('All_model_accuracy.csv')
head(df)
df[1:15, 1:5]

df2 <- read.csv('Stage_specific_accuracy.csv')
df2[1:15, 1:6]

library(tidyverse)

df_test <- df %>%
  filter(Set == "Test") %>%
  mutate(
    Time_num = parse_number(TimePoint),
    TimePoint = factor(
      TimePoint,
      levels = unique(TimePoint[order(Time_num)])
    )
  )

ggplot(df_test, aes(
  x = TimePoint,
  y = Correlation,
  color = Model,
  group = Model
)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_discrete(labels = \(x) sub("_cum$", "", x)) +
  labs(
    x = "Time point",
    y = "Predictive ability",
    color = "Model"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )



##Df2
library(tidyverse)

stages <- paste0("S", 1:12)

# G results from the previous dataset, repeated for every stage
G_test <- df %>%
  filter(Set == "Test", Model == "G") %>%
  select(TimePoint, Set, Model, Correlation, RMSE) %>%
  crossing(Stage = stages)

# Combine with test results from df2
plot_df <- df2 %>%
  filter(Set == "Test") %>%
  bind_rows(G_test) %>%
  mutate(
    Stage = factor(Stage, levels = stages),
    Time_num = parse_number(TimePoint)
  )

time_levels <- plot_df %>%
  distinct(TimePoint, Time_num) %>%
  arrange(Time_num) %>%
  pull(TimePoint)

plot_df <- plot_df %>%
  mutate(TimePoint = factor(TimePoint, levels = time_levels))

# Plot
ggplot(
  plot_df,
  aes(
    x = TimePoint,
    y = Correlation,
    color = Model,
    linetype = Model,
    group = Model
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~Stage, ncol = 4) +
  scale_x_discrete(labels = \(x) sub("_cum$", "", x)) +
  scale_color_manual(values = c(
    "G" = "black",
    "Phat" = "#0072B2",
    "G+Phat" = "#D55E00"
  )) +
  scale_linetype_manual(values = c(
    "G" = "dashed",
    "Phat" = "solid",
    "G+Phat" = "solid"
  )) +
  labs(
    x = "Time point",
    y = "Predictive ability",
    color = "Model",
    linetype = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



##

library(tidyverse)

model_order <- c(
  "G", "P", "G+P",
  "Phat_stage", "G+Phat_stage",
  "Phat_all", "G+Phat_all",
  "Phat_traits", "G+Phat_traits"
)

panel_order <- c(
  "G, P, and G+P",
  "G, Phat_stage, and G+Phat_stage",
  "G, Phat_all, and G+Phat_all",
  "G, Phat_traits, and G+Phat_traits"
)

# Assign non-G models to their comparison panels
model_data <- df %>%
  filter(Set == "Test", Model != "G") %>%
  mutate(
    Panel = case_when(
      Model %in% c("P", "G+P") ~ panel_order[1],
      Model %in% c("Phat_stage", "G+Phat_stage") ~ panel_order[2],
      Model %in% c("Phat_all", "G+Phat_all") ~ panel_order[3],
      Model %in% c("Phat_traits", "G+Phat_traits") ~ panel_order[4]
    )
  ) %>%
  filter(!is.na(Panel))

# Repeat G in all four panels
G_data <- df %>%
  filter(Set == "Test", Model == "G") %>%
  crossing(Panel = panel_order)

# Combine and set plotting order
plot_df <- bind_rows(G_data, model_data) %>%
  mutate(
    Model = factor(Model, levels = model_order),
    Panel = factor(Panel, levels = panel_order),
    Time_num = parse_number(TimePoint)
  )

time_order <- plot_df %>%
  distinct(TimePoint, Time_num) %>%
  arrange(Time_num) %>%
  pull(TimePoint)

plot_df <- plot_df %>%
  mutate(TimePoint = factor(TimePoint, levels = time_order))

ggplot(
  plot_df,
  aes(
    x = TimePoint,
    y = Correlation,
    color = Model,
    group = Model
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  facet_wrap(~Panel, ncol = 2) +
  scale_x_discrete(labels = \(x) sub("_cum$", "", x)) +
  scale_color_manual(
    breaks = model_order,
    values = c(
      "G" = "black",
      "P" = "#00A6D6",
      "G+P" = "#D89000",
      "Phat_stage" = "#CC66E3",
      "G+Phat_stage" = "#00A63C",
      "Phat_all" = "#619CFF",
      "G+Phat_all" = "#7A9E00",
      "Phat_traits" = "#F45BB5",
      "G+Phat_traits" = "#00B89C"
    ),
    guide = guide_legend(nrow = 3, byrow = TRUE)
  ) +
  labs(
    x = "Time point",
    y = "Predictive ability",
    color = "Model"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

