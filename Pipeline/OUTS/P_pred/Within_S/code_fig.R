rm(list=ls())
getwd()
library(dplyr)
library(tidyverse)

setwd("./PhD data/P3/Pipeline/Fit_models/P_pred/Within_S/")

seasons <- c("2020-21","2021-22","2022-23","2023-24","2024-25","2025-26")

all_results <- bind_rows(lapply(seasons,function(s){
  
  P <- read.csv(file.path(s,"output","Basic_models_accuracy.csv"))
  R <- read.csv(file.path(s,"output","Res","Residual_models_accuracy.csv"))
  
  P$Season <- s
  R$Season <- s
  
  bind_rows(P,R)
}))

write.csv(all_results,"All_seasons_model_accuracy.csv",row.names=FALSE)

d <- all_results %>%
  filter(
    Set=="Test",
    Model %in% c("G","P","G+P", "Residual","G+Residual")
  ) %>%
  mutate(
    Time=as.numeric(gsub("X|_cum","",Trait)),
    Trait=factor(Trait,levels=paste0("X",sort(unique(Time)),"_cum"))
  )

ggplot(d,aes(Trait,Correlation,color=Model,group=Model))+
  geom_line(linewidth=1)+
  geom_point(size=2.2)+
  facet_wrap(~Season,ncol=2)+
  labs(x="Yield time point",y="Correlation",color="Model")+
  theme_classic(base_size=12)+
  theme(
    axis.text.x=element_text(angle=45,hjust=1),
    strip.text=element_text(face="bold"),
    legend.position="bottom"
  )

ggsave("Within_season_selected_models.png",
       width=11,height=8,dpi=600)