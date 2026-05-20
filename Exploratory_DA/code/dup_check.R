
rm(list=ls())

getwd()

setwd('../Desktop/PhD data/P3/Exploratory_DA/output/')

VIs <- read.csv('../../UAV Database/Data2/VIs_Database_SpATS.csv')
Biomass <- read.csv('../../UAV Database/Data2/Biomass_Database_SpATS.csv')
head(Biomass)
head(VIs)

table(grepl('/', Biomass$Date))
table(grepl('/', VIs$Date))
VIs[!grepl('/', VIs$Date),]

VIs$Date <- trimws(VIs$Date)
sum(is.na(VIs$Date))

i1 <- grepl("/", VIs$Date)

i2 <- grepl("-", VIs$Date)

#Change format
VIs$Date[i1] <- format(
  as.Date(VIs$Date[i1], format = "%m/%d/%Y"),
  "%m/%d/%Y"
)

VIs$Date[i2] <- format(
  as.Date(VIs$Date[i2], format = "%Y-%m-%d"),
  "%m/%d/%Y"
)

table(grepl('/', VIs$Date))

#table(grepl('\\.', VIs$Genotype))
#VIs[!grepl('\\.', VIs$Genotype),]

#table(grepl('-', VIs$Genotype))
#c<-VIs[!grepl('-', VIs$Genotype),]

#table(grepl('_', VIs$Genotype))
#a <-VIs[grepl('_', VIs$Genotype),]

#Clean geno names
VIs$Genotype <- trimws(VIs$Genotype)
Biomass$Genotype <- trimws(Biomass$Genotype)

#Check duplicates
dup_VIs <- duplicated(VIs[, setdiff(names(VIs), "BLUE")])
sum(dup_VIs)
VIs[dup_VIs,]

a <- VIs[dup_VIs,]
VIs[VIs$Season=='2020-21' & VIs$SeasonTime==43 & VIs$Rep==5 & VIs$Genotype=='Radiance',]
VIs[VIs$Season=='2020-21' & VIs$SeasonTime==43 & VIs$Rep==5 & VIs$Genotype=='Festival',]
table(VIs$SeasonTime[dup_VIs], VIs$Season[dup_VIs])


dup_Bio <- duplicated(Biomass[, setdiff(names(Biomass), "BLUE")])
sum(dup_Bio)
b <- Biomass[dup_Bio,]

Biomass[Biomass$Season=='2020-21' & Biomass$SeasonTime==43 & Biomass$Rep==5 & Biomass$Genotype=='Radiance',]

