rm(list=ls())

getwd()

setwd('../Desktop/PhD data/P3/Exploratory_DA/output/')

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

#VIs[VIs$Genotype=='19.55-42'& VIs$Season=='2020-21'& VIs$SeasonTime=='43'& VIs$Trait=='NDVIrededge',]

#VIs[VIs$Genotype=='19.16-273'& VIs$Season=='2020-21'& VIs$SeasonTime=='56'& VIs$Trait=='NDVIrededge',]


#Biomass[Biomass$Genotype=='19.55-42'& Biomass$Season=='2020-21'& Biomass$SeasonTime == '56' & Biomass$Trait=='CanopyArea',]


#Biomass[Biomass$Genotype=='19.55-42'& Biomass$Season=='2020-21'& Biomass$SeasonTime=='56' & Biomass$Trait=='CanopyArea',]

#Biomass[Biomass$Genotype=='19.16-273'& Biomass$Season=='2020-21'& Biomass$SeasonTime=='56' & Biomass$Trait=='CanopyArea',]

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


row_na <- rowMeans(is.na(df[, trait_cols]))

# top rows with most NA
sort(row_na, decreasing = TRUE)[1:10]

col_na <- colMeans(is.na(df[, trait_cols]))

sort(col_na, decreasing = TRUE)[1:10]

i <- which.max(row_na)
df[i, trait_cols]



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







seasons <- unique(df$Season)

pdf("All_histograms.pdf", width = 10, height = 10)

for(s in seasons){
  
  sub <- df[df$Season == s, ]
  
  for(start in seq(1, length(trait_cols), by = 9)){
    
    end <- min(start + 8, length(trait_cols))
    
    par(mfrow = c(3,3), mar = c(2,2,2,1))
    
    for(i in trait_cols[start:end]){
      
      x <- sub[[i]]
      x <- x[!is.na(x)]
      
      if(length(x) > 0){
        hist(x,
             main = paste(i, "-", s),
             xlab = "",
             col = "lightblue",
             cex.main = 0.8)
      } else {
        plot.new()
        title(main = paste(i, "-", s))
      }
    }
  }
}

dev.off()



mean(df$RRI2[df$Season=='2020-21'], na.rm= TRUE)

mean(df$RRI2[df$Season=='2024-25'], na.rm= TRUE)


table(df$SeasonTime[df$Season=='2024-25'])
mean(df$RRI2[df$Season=='2024-25' & df$SeasonTime=='48'], na.rm= TRUE)
mean(df$RRI2[df$Season=='2024-25' & df$SeasonTime=='50'], na.rm= TRUE)
mean(df$RRI2[df$Season=='2024-25' & df$SeasonTime=='55'], na.rm= TRUE)
mean(df$RRI2[df$Season=='2024-25' & df$SeasonTime=='62'], na.rm= TRUE)
mean(df$RRI2[df$Season=='2024-25' & df$SeasonTime=='99'], na.rm= TRUE)
mean(df$RRI2[df$Season=='2024-25' & df$SeasonTime=='120'], na.rm= TRUE)
mean(df$RRI2[df$Season=='2024-25' & df$SeasonTime=='134'], na.rm= TRUE)
mean(df$RRI2[df$Season=='2024-25' & df$SeasonTime=='139'], na.rm= TRUE)

SS <- df

SS$RRi2 <- SS$RE/SS$R
mean(SS$RRi2[SS$Season=='2024-25'], na.rm= TRUE)

mean(SS$RRi2[SS$Season=='2024-25'& SS$SeasonTime=='48'], na.rm= TRUE)
mean(SS$RRi2[SS$Season=='2024-25'& SS$SeasonTime=='50'], na.rm= TRUE)
mean(SS$RRi2[SS$Season=='2024-25'& SS$SeasonTime=='55'], na.rm= TRUE)
mean(SS$RRi2[SS$Season=='2024-25'& SS$SeasonTime=='62'], na.rm= TRUE)
mean(SS$RRi2[SS$Season=='2024-25'& SS$SeasonTime=='99'], na.rm= TRUE)
mean(SS$RRi2[SS$Season=='2024-25'& SS$SeasonTime=='120'], na.rm= TRUE)
mean(SS$RRi2[SS$Season=='2024-25'& SS$SeasonTime=='134'], na.rm= TRUE)
mean(SS$RRi2[SS$Season=='2024-25'& SS$SeasonTime=='139'], na.rm= TRUE)




head(VIs1)
VIs2 <- VIs[VIs$Season=='2021-22',]
head(VIs2)
VIs2$RRi2 <- VIs2$RE/VIs2$R
mean(VIs2$RRi2[VIs2$SeasonTime=='48'], na.rm= TRUE)
mean(VIs2$RRi2[VIs2$SeasonTime=='50'], na.rm= TRUE)

VIs3 <- reshape(VIs2,idvar = c("Season", "Date", "SeasonTime", "Rep", "Genotype"), timevar = "Trait", direction = "wide")
head()

VIs2[VIs2$Genotype=='19.47-180',] 
pheno <- read.csv('../../../P1/P1/1.data.combining/output/Phenos.csv')
