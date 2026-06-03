rm(list=ls())

getwd()
setwd('../Desktop/PhD data/P3/Pipeline/3.Phenomic_profiling/output/')

Pheno <- read.csv('../../2.BLUEs/output/Phenomic_BLUEs_wide.csv')
Res <- read.csv('../../2.BLUEs/output/Response_BLUEs_wide.csv')
head(Res)
table(Res$Season)

head(Pheno)
table(Pheno$Season , Pheno$Stage)

hist(Res$X3)
#pheno0<- read.csv('../../1. Standardizing_phenomics/output/Pheno_all.csv')
#res0<- read.csv('../../../Preparing_data/output/Response.csv')
#Res[Res$Genotype=='19.11-199',]
#length(unique(res0$Genotype[res0$Season=='2020-21']))
#length(unique(Res$Genotype[Res$Season=='2020-21']))

### CUMULATIVE##############

summary(as.matrix(Res[,3:ncol(Res)]))
range(as.matrix(Res[,3:ncol(Res)]))

### pmax= parallel maximum (pick maximum value eg. -5, 0 - it will pick 0)

tmp <- pmax(as.matrix(Res[,3:ncol(Res)]),  0)
head(tmp)
head(Res)
cumY <- t(apply(tmp, 1, cumsum))

colnames(cumY) <- paste0(colnames(Res[,3:ncol(Res)]), "_cum")

Res1 <- cbind(Res, cumY)
head(Res1)

Pheno1 <- Pheno
Pheno1 <- Pheno1[, colMeans(is.na(Pheno1)) <= 0.2]
Pheno1 <- Pheno1[rowMeans(is.na(Pheno1)) <= 0.5, ]


image(t(is.na(Pheno1)))
dim(Pheno1)
dim(Pheno)

############################

Stage <- unique(Pheno1$Stage)
Season <- unique(Pheno1$Season)
head(Pheno1)
Traits <- colnames(Pheno1)[4:ncol(Pheno1)]


##########Correlation bwn Phenomics covariates#####


dir.create("Correlation_Matrices_Pheno", showWarnings = FALSE)
library(corrplot)

pdf("Correlation_Plots.pdf", width = 10, height = 8)

for (s in Season) {
  for (i in Stage) {
    
    tmp <- subset(Pheno1, Season == s & Stage == i)
    if (nrow(tmp) < 2) next
    
    corr_mat1 <- cor(tmp[,4:52], use = "pairwise.complete.obs")
    corr_mat2 <- cor(tmp[,53:99], use = "pairwise.complete.obs")
    
    write.csv(corr_mat1, paste0("Correlation_Matrices_Pheno/Corr1_", s, "_", i, ".csv"))
    write.csv(corr_mat2, paste0("Correlation_Matrices_Pheno/Corr2_", s, "_", i, ".csv"))
    
    corrplot(corr_mat1, method="color", type="upper",
             tl.col="black", main=paste("Season:", s, "Stage:", i, "- Set 1"))
    
    corrplot(corr_mat2, method="color", type="upper",
             tl.col="black", main=paste("Season:", s, "Stage:", i, "- Set 2"))
  }
}

dev.off()


sum(is.na(Pheno1[,4:100]))

table(Pheno1$Season, Pheno1$Stage)
table(Pheno$Season, Pheno$Stage)



#############Response and phenomics


library(corrplot)

dir.create("Corr_Matrices_Pheno_Res", showWarnings = FALSE)

pheno1 <- names(Pheno1)[4:51]
pheno2 <- names(Pheno1)[52:99]

pdf("Corr_Pheno_Res_All.pdf", width = 14, height = 11)

for(s in unique(Pheno1$Season)){
  for(i in unique(Pheno1$Stage)){
    
    tmp <- merge(subset(Pheno1, Season == s & Stage == i),Res1, by = c("Season","Genotype"))
    
    if(nrow(tmp) < 2) next
    
    resp <- grep("^X", names(tmp), value = TRUE)
    
    corr_mat1 <- cor(tmp[, pheno1], tmp[, resp], use = "pairwise.complete.obs")
    
    corr_mat2 <- cor(tmp[, pheno2], tmp[, resp], use = "pairwise.complete.obs")
    
    
    par(mfrow = c(1,2), mar = c(4,4,8,4), oma = c(0,0,0,0))    
    corrplot(corr_mat1, is.corr = FALSE, method = "color", tl.col = "black", tl.cex = 0.5)
    title(main = paste("Season", s, "Stage", i),  line = 7,  cex.main = 1    )
    
    corrplot(corr_mat2, is.corr = FALSE, method = "color", tl.col = "black", tl.cex = 0.5)
    title(main = paste("Season", s, "Stage", i),  line = 7,  cex.main = 1    )
    
    cat("Done:", s, i, "\n")
  }
}

dev.off()

 #Pheno[Pheno$Season=='2021-22'& Pheno$Stage == 'S12',]


#####  Check pheno and pheno1....missing information ###########

head(Pheno)
dim(Pheno)
dim(Pheno1)
table(Pheno$Season , Pheno$Stage)
table(Pheno1$Season, Pheno1$Stage)
sum(is.na(Pheno))
pdf("na_plot_Pheno1.pdf", width = 8, height = 6)
image(t(is.na(Pheno1)))
dev.off()


a <- Pheno[Pheno$Season=='2020-21' & Pheno$Stage=='S11',]
nearest_table <- read.csv('../../../Preparing_data/output/Phenomic_Stages1.csv')
in_pHeno <- read.csv('../../../Preparing_data/data/BVI_mean.csv')
b <- in_pHeno[in_pHeno$Season=='2020-21'& in_pHeno$SeasonTime== '147',]
####it means no data in season time 147,

#S11/S12= 145/155 ~~147 ##no data for VIs


write.csv(Pheno, file='Pheno.csv', row.names= FALSE)
write.csv(Pheno1, file='Pheno_cl.csv', row.names= FALSE)
write.csv(Res1, file= 'Res.csv', row.names= FALSE)
################################# Matching data ####################


load('../../../Preparing_data/data/X4.rda')
X4[1:5, 1:9]
head(Pheno1)
Geno <- X4
Geno[1:5, 1:9]
head(Res1)


cat("Geno:  ", nrow(Geno), "genotypes x", ncol(Geno), "markers\n")
cat("Pheno1:", nrow(Pheno1), "rows |", length(unique(Pheno1$Genotype)), "genotypes\n")
cat("Res1:  ", nrow(Res1), "rows |", length(unique(Res1$Genotype)), "genotypes\n")

#### Pheno
meta_cols  <- 1:3
orig_cols  <- 4:51
std_cols   <- 52:99

orig_names <- colnames(Pheno1)[orig_cols]
orig_names
std_names  <- colnames(Pheno1)[std_cols]
std_names

#####
for(s in unique(Pheno1$Season)){
  sub <- Pheno1[Pheno1$Season == s, ]
  cat(s, ":", paste(sort(unique(sub$Stage)), collapse = " "), "\n")
}
####Check again dup###
dup_check <- duplicated(Pheno1[, c("Season", "Stage", "Genotype")])
Pheno1_clean <- Pheno1[!dup_check, ]
dim(Pheno1)
dim(Pheno1_clean)


### Wide Pheno   ###########################################################################################

col_names_std  <- as.vector(outer(stages_list, std_names,  paste, sep = "_"))
col_names_orig <- as.vector(outer(stages_list, orig_names, paste, sep = "_"))

#  empty dataframes 
Pheno_wide_std  <- data.frame(matrix(NA, nrow = nrow(Pheno1_clean), ncol = 2 + length(col_names_std)))
Pheno_wide_orig <- data.frame(matrix(NA, nrow = nrow(Pheno1_clean), ncol = 2 + length(col_names_orig)))

colnames(Pheno_wide_std)  <- c("Season", "Genotype", col_names_std)
colnames(Pheno_wide_orig) <- c("Season", "Genotype", col_names_orig)

# Get unique Season x Genotype combinations
sea_geno        <- unique(Pheno1_clean[, c("Season", "Genotype")])
sea_geno        <- sea_geno[order(sea_geno$Season, sea_geno$Genotype), ]
rownames(sea_geno) <- NULL

Pheno_wide_std  <- data.frame(sea_geno, matrix(NA, nrow = nrow(sea_geno), ncol = length(col_names_std)))
Pheno_wide_orig <- data.frame(sea_geno, matrix(NA, nrow = nrow(sea_geno), ncol = length(col_names_orig)))

colnames(Pheno_wide_std)  <- c("Season", "Genotype", col_names_std)
colnames(Pheno_wide_orig) <- c("Season", "Genotype", col_names_orig)

# Fill stage by stage
for(stg in stages_list){
  
  stg_data <- Pheno1_clean[Pheno1_clean$Stage == stg, ]
  
  if(nrow(stg_data) == 0) next
  
  std_here  <- paste0(stg, "_", std_names)
  orig_here <- paste0(stg, "_", orig_names)
  
  # Match rows by Season + Genotype
  idx <- match(paste(stg_data$Season, stg_data$Genotype),  paste(Pheno_wide_std$Season, Pheno_wide_std$Genotype))
  
  Pheno_wide_std [idx, std_here]  <- stg_data[, std_cols]
  Pheno_wide_orig[idx, orig_here] <- stg_data[, orig_cols]
}

# Check
dim(Pheno_wide_std)
dim(Pheno_wide_orig)
Pheno_wide_std[1:5, 1:10]

table(Pheno_wide_std$Season)
table(Pheno_wide_orig$Season)



###################################################################
# MATCH PHENO x RES1 x GENO


# Common genotypes across all three
genos_pheno <- unique(Pheno_wide_std$Genotype)
genos_res   <- unique(Res1$Genotype)
genos_geno  <- rownames(Geno)

common_genos <- Reduce(intersect, list(genos_pheno, genos_res, genos_geno))
length(common_genos)

# Subset
Pheno_std_matched  <- Pheno_wide_std[Pheno_wide_std$Genotype  %in% common_genos, ]
Pheno_orig_matched <- Pheno_wide_orig[Pheno_wide_orig$Genotype %in% common_genos, ]
Res1_matched       <- Res1[Res1$Genotype %in% common_genos, ]
Geno_matched       <- Geno[common_genos, ]
dim(Geno_matched)

# Sort all by Season then Genotype
Pheno_std_matched  <- Pheno_std_matched [order(Pheno_std_matched$Season, Pheno_std_matched$Genotype), ]
Pheno_orig_matched <- Pheno_orig_matched[order(Pheno_orig_matched$Season, Pheno_orig_matched$Genotype), ]
Res1_matched       <- Res1_matched[order(Res1_matched$Season, Res1_matched$Genotype), ]

dim(Pheno_orig_matched)
dim(Pheno_std_matched)
dim(Res1_matched)
dim(Geno_matched)

Pheno_orig_matched[1:8, 1:7]
nrow(Pheno_std_matched) == nrow(Res1_matched)


key_pheno <- paste(Pheno_std_matched$Season, Pheno_std_matched$Genotype)
key_res   <- paste(Res1_matched$Season, Res1_matched$Genotype)

# How many in each
length(key_pheno)
length(key_res)

length(setdiff(key_pheno, key_res))
length(setdiff(key_res, key_pheno))

# Common keys
length(intersect(key_pheno, key_res))


# Keep only rows present in BOTH
common_keys <- intersect(key_pheno, key_res)

Pheno_std_matched  <- Pheno_std_matched[key_pheno %in% common_keys, ]
Pheno_orig_matched <- Pheno_orig_matched[paste(Pheno_orig_matched$Season,
                                               Pheno_orig_matched$Genotype) %in% common_keys, ]
Res1_matched       <- Res1_matched[key_res %in% common_keys, ]

# Re-sort
Pheno_std_matched  <- Pheno_std_matched[order(Pheno_std_matched$Season, Pheno_std_matched$Genotype), ]
Pheno_orig_matched <- Pheno_orig_matched[order(Pheno_orig_matched$Season, Pheno_orig_matched$Genotype), ]
Res1_matched       <- Res1_matched[order(Res1_matched$Season, Res1_matched$Genotype), ]


length(unique(Pheno_orig_matched$Genotype))

length(unique(Res1_matched$Genotype))
length(rownames(Geno_matched))

table(Pheno_std_matched$Season)
table(Res1_matched$Season)



####
dir.create('Final_outs',showWarnings = FALSE)

saveRDS(Pheno_std_matched, "Final_outs/Pheno_std_matched.rds")
saveRDS(Pheno_orig_matched, "Final_outs/Pheno_orig_matched.rds")
saveRDS(Res1_matched, "Final_outs/Res1_matched.rds")
saveRDS(Geno_matched, "Final_outs/Geno_matched.rds")

saveRDS(common_genos, "Final_outs/common_genos.rds")


####
# Save as CSV
write.csv(Pheno_std_matched,  "Final_outs/Pheno_std_matched.csv",  row.names = FALSE)
write.csv(Pheno_orig_matched, "Final_outs/Pheno_orig_matched.csv", row.names = FALSE)
write.csv(Res1_matched, "Final_outs/Res1_matched.csv",       row.names = FALSE)
write.csv(Geno_matched, "Final_outs/Geno_matched.csv",       row.names = TRUE)   
write.csv(data.frame(Genotype = common_genos), "Final_outs/common_genos.csv", row.names = FALSE)
