rm(list=ls())

getwd()

setwd('../../../../../../Desktop/PhD data/P3/Exploratory_DA/output/')

VIs <- read.csv('../../UAV Database/VIs_Database_SpATS.csv')
Biomass <- read.csv('../../UAV Database/Biomass_Database_SpATS.csv')

head(VIs)
head(Biomass)

table(VIs$Trait)
length(unique(VIs$Genotype))
dim(VIs)
table(VIs$Season)
table(VIs$Rep)
table(Biomass$Rep)

VIs[VIs$Genotype=='19.55-42'& VIs$Season=='2020-21'& VIs$SeasonTime=='43'& VIs$Trait=='NDVIrededge',]

VIs[VIs$Genotype=='19.16-273'& VIs$Season=='2020-21'& VIs$SeasonTime=='56'& VIs$Trait=='NDVIrededge',]


#Biomass[Biomass$Genotype=='19.55-42'& Biomass$Season=='2020-21'& Biomass$SeasonTime == '56' & Biomass$Trait=='CanopyArea',]

a <- unique(VIs$Genotype[VIs$Season == "2020-21"])
b <- unique(Biomass$Genotype[Biomass$Season == "2020-21"])
c<- intersect(a,b)
setdiff(a,b)


Biomass[Biomass$Genotype=='19.55-42'& Biomass$Season=='2020-21'& Biomass$SeasonTime=='56' & Biomass$Trait=='CanopyArea',]

Biomass[Biomass$Genotype=='19.16-273'& Biomass$Season=='2020-21'& Biomass$SeasonTime=='56' & Biomass$Trait=='CanopyArea',]



head(Biomass)
head(VIs)

Biomass1 <- reshape(Biomass,idvar = c("Season", "Date", "SeasonTime", "Rep", "Genotype"), timevar = "Trait", direction = "wide")
names(Biomass1) <- sub("^BLUE\\.", "", names(Biomass1))
head(Biomass1)

VIs1 <- reshape(VIs,idvar = c("Season", "Date", "SeasonTime", "Rep", "Genotype"), timevar = "Trait", direction = "wide")
names(VIs1) <- sub("^BLUE\\.", "", names(VIs1))
head(VIs1)

any(is.na(VIs1))
any(is.na(Biomass1))

sum(is.na(Biomass1$CanopyArea))
sum(is.na(Biomass1$CanopyVolume))
sum(is.na(Biomass1$Height))
sum(is.na(Biomass1$HeightSD))

Biomass1[is.na(Biomass1$CanopyArea), ]
Biomass1[is.na(Biomass1$CanopyVolume), ]
Biomass1[is.na(Biomass1$Height), ]
Biomass1[is.na(Biomass1$HeightSD), ]

sum(is.na(VIs1$NDVIrededge))
sum(is.na(VIs1$NDVI))
sum(is.na(VIs1$OSAVI))
sum(is.na(VIs1$BNDVI))

VIs1[is.na(VIs1$NDVIrededge), ]
VIs1[is.na(VIs1$NDVI), ]
VIs1[is.na(VIs1$OSAVI), ]
VIs1[is.na(VIs1$BNDVI), ]


df <- merge(Biomass1, VIs1, by = c("Season", "Date", "SeasonTime", "Rep", "Genotype"), all = TRUE)
head(df)

write.csv(df, file = 'BVI.csv', row.names= F)

table(df$Season, df$SeasonTime, df$Rep)
str(df)
table(df$SeasonTime[df$Season=='2020-21'], df$Rep[df$Season=='2020-21'])
table(df$SeasonTime[df$Season=='2021-22'])
table(df$SeasonTime[df$Season=='2022-23'])
table(df$SeasonTime[df$Season=='2023-24'])
table(df$SeasonTime[df$Season=='2024-25'])
table(df$SeasonTime[df$Season=='2025-26'])

id_cols <- c("Season","Date","SeasonTime","Rep","Genotype")
trait_cols <- setdiff(names(df), id_cols)
#image(is.na(df[, trait_cols]),xlab = "Traits",ylab = "Observations", main = "NA pattern")


seasons <- unique(df$Season)

par(mfrow = c(1,1))

for(s in seasons){
  sub <- df[df$Season == s, ]

 image(t(is.na(sub[, trait_cols])),
              axes = FALSE,
              main = paste("NA pattern -", s))
  
  readline("Press Enter for next season...")
}

par(mfrow=c(1,1))  
image(t(is.na(df[, trait_cols])),
      axes = FALSE,
      main = paste("NA pattern -", 'full'))

image((is.na(df[, trait_cols])),
      axes = FALSE,
      main = paste("NA pattern -", "full"))


colSums(is.na(df[, trait_cols])) == nrow(df)
row_all_na <- apply(is.na(df[, trait_cols]), 1, all)
sum(row_all_na)

par(mfrow=c(3,3))  

for(i in trait_cols[1:9]){
  hist(df[[i]],
       main = i,
       xlab = "",
       col = "lightblue")
}

for(i in trait_cols[10:19]){
  hist(df[[i]],
       main = i,
       xlab = "",
       col = "lightblue")
}

for(i in trait_cols[20:29]){
  hist(df[[i]],
       main = i,
       xlab = "",
       col = "lightblue")
}

for(i in trait_cols[30:39]){
  hist(df[[i]],
       main = i,
       xlab = "",
       col = "lightblue")
}

for(i in trait_cols[40:49]){
  hist(df[[i]],
       main = i,
       xlab = "",
       col = "lightblue")
}


pheno <- read.csv('../../../P1/P1/1.data.combining/output/Phenos.csv')
