rm(list=ls())

getwd()
setwd('../../Genotyping_data/output/')

pheno <- read.csv('../data/BVI_mean.csv')
load('../data/X4.rda')
geno <- X4

head(pheno)
geno[1:5, 1:5]


g_pheno <- unique(pheno$Genotype)
length(g_pheno)
 
dim(geno)
colnames(geno)
rownames(geno)


a <- intersect(rownames(geno), g_pheno)
length(a)

pheno1 <- pheno[pheno$Genotype %in% a,]
g_pheno1 <- unique(pheno1$Genotype)

geno1<- geno[rownames(geno) %in% a,]
dim(geno1)
 