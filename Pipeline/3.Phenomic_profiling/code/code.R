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

### With response ####

library(corrplot)

for(s in unique(Pheno1$Season)){
  for(i in unique(Pheno1$Stage)){
    
    tmp <- merge(subset(Pheno1, Season==s & Stage==i), Res, by=c("Season","Genotype"))
    
    Pheno1 <- colnames(Pheno1)[4:52]
    resp  <- grep("^X", names(tmp), value=TRUE)
    
    corr_mat <- cor(tmp[,Pheno1], tmp[,resp], use="pairwise.complete.obs")
    
    write.csv(corr_mat, paste0("Corr_",s,"_",i,".csv"))
    
    pdf(paste0("Corr_",s,"_",i,".pdf"), width=10, height=8)
    corrplot(corr_mat, is.corr=FALSE, method="color", tl.cex=.6)
    dev.off()
    
    cat("Done:", s, i, "\n")
  }
}
