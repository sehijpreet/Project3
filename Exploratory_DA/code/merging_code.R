
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

# rows with slash format
i1 <- grepl("/", VIs$Date)

# rows with dash format
i2 <- grepl("-", VIs$Date)

# convert separately
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

VIs$Genotype <- trimws(VIs$Genotype)
Biomass$Genotype <- trimws(Biomass$Genotype)

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

#####Fr Screenshot####

# VIs - keep all duplicate groups together
dup_all_VIs <- duplicated(VIs[, setdiff(names(VIs), "BLUE")]) |
  duplicated(VIs[, setdiff(names(VIs), "BLUE")], fromLast = TRUE)
a1 <- VIs[dup_all_VIs, ]
a1 <- a1[order(a1$Season, a1$SeasonTime, a1$Rep, a1$Genotype, a1$Trait), ]
head(a1, 50)


dup_all_Bio <- duplicated(Biomass[, setdiff(names(Biomass), "BLUE")]) |
  duplicated(Biomass[, setdiff(names(Biomass), "BLUE")], fromLast = TRUE)
b1 <- Biomass[dup_all_Bio, ]
b1 <- b1[order(b1$Season, b1$SeasonTime, b1$Rep, b1$Genotype, b1$Trait), ]
head(b1, 50)

###Mean

# VIs - take mean of duplicate BLUE values and keep one row
VIs_mean <- aggregate(BLUE ~ .,data = VIs,FUN = mean,na.rm = TRUE)
# check duplicates again
sum(duplicated(VIs_mean[, setdiff(names(VIs_mean), "BLUE")]))

# Biomass - take mean of duplicate BLUE values and keep one row
Biomass_mean <- aggregate(BLUE ~ .,data = Biomass,FUN = mean, na.rm = TRUE)
# check duplicates again
sum(duplicated(Biomass_mean[, setdiff(names(Biomass_mean), "BLUE")]))


write.csv(Biomass_mean, file= 'Biomass_mean.csv', row.names=FALSE)
write.csv(VIs_mean, file = 'VIs_mean.csv', row.names= FALSE)

library(dplyr)
library(tidyr)

Biomass1 <- Biomass_mean %>%
  spread(key = Trait, value = BLUE)

VIs1 <- VIs_mean %>% 
  spread(key = Trait, value = BLUE)

head(Biomass1)
head(VIs1)

df <- merge(Biomass1, VIs1, by = c("Season", "Date", "SeasonTime", "Rep", "Genotype"), all = TRUE)
head(df)

write.csv(df, file = 'BVI_mean.csv', row.names= F)

df1 <- read.csv('../../Exploratory_DA/output/BVI.csv')

df2 <- read.csv('../../Exploratory_DA/output/BVI_mean.csv')


# Check NA pattern from column 6 to last column

image(t(is.na(df[, 6:ncol(df)])),
  axes = FALSE,
  col = c("white", "black"),
  main = "NA Pattern")

# Add axis labels
axis(1,at = seq(0, 1, length.out = nrow(df)),labels = FALSE)

axis(2,at = seq(0, 1, length.out = ncol(df[, 6:ncol(df)])),labels = names(df)[6:ncol(df)],
     las = 2,cex.axis = 0.5)



Biomass_mean<- read.csv('Biomass_mean.csv')

VIs_mean <- read.csv('VIs_mean.csv')

#####  Reshape _ picking Biomass_mean.csv#####  Reshape _ picking first value from duplicates###
Biomass2 <- reshape(Biomass_mean, idvar = c("Season", "Date", "SeasonTime", "Rep", "Genotype"), timevar = "Trait", direction = "wide")
names(Biomass2) <- sub("^BLUE\\.", "", names(Biomass2))
head(Biomass2)

VIs2 <- reshape(VIs_mean, idvar = c("Season", "Date", "SeasonTime", "Rep", "Genotype"), timevar = "Trait", direction = "wide")
names(VIs2) <- sub("^BLUE\\.", "", names(VIs2))
head(VIs2)


df2 <- merge(Biomass2, VIs2, by = c("Season", "Date", "SeasonTime", "Rep", "Genotype"), all = TRUE)
head(df2)

write.csv(df2, file = 'BVI.csv', row.names= F)


image(t(is.na(df2[, 6:ncol(df2)])),
      axes = FALSE,
      col = c("white", "black"),
      main = "NA Pattern")

identical(df1, df2)
all.equal(df1, df2)





###################################### DO NOT RUN #####################

# subset RRI2 rows for 2024-25
RRI2_2024_25 <- VIs[VIs$Season == "2024-25" & VIs$Trait == "RRI2",c("Date", "SeasonTime", "Rep", "Genotype", "BLUE")]

# fetch R values
R_vals <- VIs[VIs$Season == "2024-25" & VIs$Trait == "R", c("Date", "SeasonTime", "Rep", "Genotype", "BLUE")]

# fetch REdge values
REdge_vals <- VIs[VIs$Season == "2024-25" & VIs$Trait == "RE", c("Date", "SeasonTime", "Rep", "Genotype", "BLUE")]

# rename BLUE columns
names(R_vals)[5] <- "R"
names(REdge_vals)[5] <- "RE"

# merge R values
RRI2_2024_25 <- merge(RRI2_2024_25, R_vals, by = c("Date", "SeasonTime", "Rep", "Genotype"),all.x = TRUE)

# merge REdge values
RRI2_2024_25 <- merge(RRI2_2024_25,REdge_vals,by = c("Date", "SeasonTime", "Rep", "Genotype"),all.x = TRUE)

RRI2_2024_25$BLUE2 <- RRI2_2024_25$RE/RRI2_2024_25$R

cor(RRI2_2024_25$BLUE2, RRI2_2024_25$BLUE, na.rm= TRUE)
mean(RRI2_2024_25$BLUE2, na.rm= TRUE)
hist(RRI2_2024_25$BLUE, na.rm= TRUE)



##############################################################

0.90 *54

head(df1)
df1[,6]

# proportion of NA per row
row_na <- rowMeans(is.na(df1[, 6:ncol(df1)]))

# Remove rows with >=90% NA
df_up1 <- df1[row_na < 0.90, ]
df_up2 <- df1[row_na < 0.95, ]

image(t(is.na(df_up2[, 6:ncol(df_up2)])),
      axes = FALSE,
      col = c("white", "black"),
      main = "NA Pattern")


image(t(is.na(df_up1[, 6:ncol(df_up1)])),
      axes = FALSE,
      col = c("white", "black"),
      main = "NA Pattern (90% clean)")



