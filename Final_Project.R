####################################
##START OF SUBSETTING FOR DISEASES##
##START OF SUBSETTING FOR DISEASES##
##START OF SUBSETTING FOR DISEASES##
####################################
setwd('C:/Users/Igor/Desktop/Bioinformatyka I rok II stopnia/III sem/Modelowanie Zaawansowane/Projekt')

library(tidyverse)
library(data.table)
library(dplyr)


cellosaurus <- read.csv(file = 'Cellosaurus_Information.csv', sep = '\t')

#############################################################################
######Subset and write the expression data on basis of the cancer type#######
#############################################################################
x <- as.character(cellosaurus$Diseases)
y <- strsplit(x,";")
z <- vector()
z <- c(q,length(y))

for(i in 1:length(y)){
  z[i] <- (y[[i]][3])
}
#unique(z)

diseases <-cellosaurus[,c(1,4,6)]
diseases$Diseases_extracted <-as.character(z)

disease_vector <- vector()
disease_vector <- c(disease_vector, length(diseases$Diseases_extracted))

cancer <- c("lung", "Lung", "Bronchioloalveolar", "Bronchogenic")
for(i in 1:length(diseases$Diseases_extracted)){
  disease_vector[i] <- as.logical(any(is.element(strsplit(diseases$Diseases_extracted, " ")[[i]], cancer)))
}

cellosaurus.cancer <- diseases[as.logical(disease_vector),]
cellosaurus.cancer <- cellosaurus.cancer[,-c(2,3)]

pheno <- cellosaurus.cancer

expression <- read.csv(file = 'Expression_split_files/Expression-1.csv', colClasses="character", header = TRUE, sep = ',')
cancer_expression <- left_join(cellosaurus.cancer, expression, by = "CellLineName_Cellosaurus")
cancer_expression <- subset(cancer_expression, Source_Type != '')

for(i in 2:10){
  path <- paste("Expression_split_files/Expression-",i,".csv", sep = '')
  expression_file <- read.csv(file = path, colClasses="character", header = TRUE, sep = ',')
  expression_i_lung <- left_join(cellosaurus.cancer, expression_file, by = "CellLineName_Cellosaurus")
  expression_i_lung <- expression_i_lung[which(expression_i_lung$Source_Type != 'NA'),]
  cancer_expression <- rbind(cancer_expression, expression_i_lung)
  
}

cancer_path <- paste("Expression_split_files/expression_", cancer[1],".csv", sep = '')

write.csv(cancer_expression, file = cancer_path)

####################################
###END OF SUBSETTING FOR DISEASES###
###END OF SUBSETTING FOR DISEASES###
###END OF SUBSETTING FOR DISEASES###
####################################

#########################
###THE CANCER PIPELINE###
#########################
library(plyr)
library(reshape)

###EXPLORATORY ANALYSIS

#Upload the data
setwd('C:/Users/Igor/Desktop/Bioinformatyka I rok II stopnia/III sem/Modelowanie Zaawansowane/Projekt')
cancer_types <- c("lung")

for(i in 1:length(cancer_types)){
  cancer_path <- paste("Expression_split_files/expression_", cancer_types[i], ".csv", sep = "")
  cancer_data <- read.csv(file = cancer_path, sep = ",")
  
  #retrive only those genes that all patients have#
  exp_data <- cancer_data[,-c(1,6,8,9,10)]
  
  source_types <- unique(exp_data$Source_Type)
  for(j in 1:length(source_types)){
    expression_data_by_source <- exp_data[which(exp_data$Source_Type == source_types[j]),]
    source_genes <- as.data.frame(unique(expression_data_by_source$HGNC))
    
    source_casted <- cast(expression_data_by_source, formula = "CellLineName_Cellosaurus + Diseases_extracted + ~ HGNC", fun.aggregate = mean, value = "Exp_Value")
    rownames(source_casted) <- source_casted$CellLineName_Cellosaurus
    source_patients <- as.data.frame(source_casted$CellLineName_Cellosaurus)
    colnames(source_patients) <- "patient"
    
    source_path <- paste("Expression_split_files/Data_analysis/", cancer_types[i], "_", source_types[j],".csv", sep = '')
    write.csv(source_casted, file = source_path)
  }
}

#Thyroid
CCLE <- read.csv(file = "Expression_split_files/Data_analysis/thyroid_CCLE.csv", sep = ",")
COSMIC <- read.csv(file = "Expression_split_files/Data_analysis/thyroid_COSMIC.csv", sep = ",")

thyroid_gene_CCLE <- colnames(CCLE)
thyroid_gene_COSMIC <- colnames(COSMIC)

genes_to_remain <- intersect(thyroid_gene_CCLE, thyroid_gene_COSMIC)
length(genes_to_remain)
CCLE.thyroid <- CCLE[,genes_to_remain]
COSMIC.thyroid <- COSMIC[,genes_to_remain]

#Kidney
CCLE <- read.csv(file = "Expression_split_files/Data_analysis/kidney_CCLE.csv", sep = ",")
COSMIC <- read.csv(file = "Expression_split_files/Data_analysis/kidney_COSMIC.csv", sep = ",")
NCI60 <- read.csv(file = "Expression_split_files/Data_analysis/kidney_NCI60.csv", sep = ",")

kidney_gene_CCLE <- colnames(CCLE)
kidney_gene_COSMIC <- colnames(COSMIC)
kidney_gene_NCI60 <- colnames(NCI60)

genes_to_remain <- intersect(kidney_gene_CCLE, kidney_gene_COSMIC)
genes_to_remain <- intersect(genes_to_remain, kidney_gene_NCI60)
length(genes_to_remain)
CCLE.kidney <- CCLE[,genes_to_remain]
COSMIC.kidney <- COSMIC[,genes_to_remain]
NCI60.kidney <- NCI60[,genes_to_remain]

#Pancreas
CCLE <- read.csv(file = "Expression_split_files/Data_analysis/pancreas_CCLE.csv", sep = ",")
COSMIC <- read.csv(file = "Expression_split_files/Data_analysis/pancreas_COSMIC.csv", sep = ",")

pancreas_gene_CCLE <- colnames(CCLE)
pancreas_gene_COSMIC <- colnames(COSMIC)

genes_to_remain <- intersect(pancreas_gene_CCLE, pancreas_gene_COSMIC)
length(genes_to_remain)
CCLE.pancreas <- CCLE[,genes_to_remain]
COSMIC.pancreas <- COSMIC[,genes_to_remain]

#Lung
CCLE <- read.csv(file = "Expression_split_files/Data_analysis/lung_CCLE.csv", sep = ",")
COSMIC <- read.csv(file = "Expression_split_files/Data_analysis/lung_COSMIC.csv", sep = ",")
NCI60 <- read.csv(file = "Expression_split_files/Data_analysis/lung_NCI60.csv", sep = ",")

lung_gene_CCLE <- colnames(CCLE)
lung_gene_COSMIC <- colnames(COSMIC)
lung_gene_NCI60 <- colnames(NCI60)

genes_to_remain <- intersect(lung_gene_CCLE, lung_gene_COSMIC)
genes_to_remain <- intersect(genes_to_remain, lung_gene_NCI60)
length(genes_to_remain)
CCLE.lung <- CCLE[,genes_to_remain]
COSMIC.lung <- COSMIC[,genes_to_remain]
NCI60.lung <- NCI60[,genes_to_remain]


write.csv(CCLE.thyroid, file = "Expression_split_files/Data_analysis/Same_genes/CCLE_thyroid.csv")
write.csv(COSMIC.thyroid, file = "Expression_split_files/Data_analysis/Same_genes/COSMIC_thyroid.csv")

write.csv(CCLE.pancreas, file = "Expression_split_files/Data_analysis/Same_genes/CCLE_pancreas.csv")
write.csv(COSMIC.pancreas, file = "Expression_split_files/Data_analysis/Same_genes/COSMIC_pancreas.csv")

write.csv(CCLE.kidney, file = "Expression_split_files/Data_analysis/Same_genes/CCLE_kidney.csv")
write.csv(COSMIC.kidney, file = "Expression_split_files/Data_analysis/Same_genes/COSMIC_kidney.csv")
write.csv(NCI60.kidney, file = "Expression_split_files/Data_analysis/Same_genes/NCI60_kidney.csv")

write.csv(CCLE.lung, file = "Expression_split_files/Data_analysis/Same_genes/CCLE_lung.csv")
write.csv(COSMIC.lung, file = "Expression_split_files/Data_analysis/Same_genes/COSMIC_lung.csv")
write.csv(NCI60.lung, file = "Expression_split_files/Data_analysis/Same_genes/NCI60_lung.csv")

################
###CLUSTERING###
################

library(Biobase)
library(GEOquery)
library(data.table)
library(broom)
library(tidyverse)
library(ggplot2)
library(BBmisc)
library(irlba)
library(Rtsne)
library(jackstraw)
library(stringi)
library(meanShiftR)
library(dbscan)

?svd

setwd('C:/Users/Igor/Desktop/Bioinformatyka I rok II stopnia/III sem/Modelowanie Zaawansowane/Projekt')

cancer <- c("lung")
source <- c("COSMIC")

for(i in 1:length(cancer)){
  for(j in 1:length(source)){
    edata_path <- paste("Expression_split_files/Data_analysis/Same_genes/",source[j],"_",cancer[i], ".csv", sep = "")
    edata <- read.csv(file = edata_path)
    rownames(edata) <- edata$CellLineName_Cellosaurus
    edata_pheno <- edata$Diseases_extracted
    edata.raw <- edata[,-c(1,2,3,4)]
    edata.norm <- scale(x = edata.raw, center = TRUE, scale = TRUE)
    
    edata.svd <- svd(t(edata.norm))
    #edata.tsvd <- irlba(t(edata.norm), nv = 12)
    #edata.tsne <- Rtsne(t(edata.norm),pca=TRUE,perplexity=30)
    
    
    edata.jackstraw = jackstraw_pca(t(edata.norm), r=2, s=round(ncol(edata.norm)*.1), B=10)
    
    
  }
}

par(mfrow=c(1,2))
plot(edata.svd$d, pch=20, ylab="Singular values")
plot(edata.svd$d^2/sum(edata.svd$d^2)*100, pch=20, ylab="% variance explained")
abline(v=40, col = "red", lty=2)
v.limit <- 50
length(unique(edata$Diseases_extracted))

raw.kmeans <- kmeans(edata.raw, centers = 13, algorithm = "Lloyd")
norm.kmeans <- kmeans(edata.norm, centers = 13, algorithm = "Lloyd")
svd.kmeans <- kmeans(edata.svd$v[,1:v.limit], centers = 13, algorithm = "Lloyd")

#raw.meanshift <- meanShift(edata.raw,edata.raw)
norm.meanshift <- meanShift(edata.norm,edata.norm)
svd.meanshift <- meanShift(edata.svd$v[,1:v.limit],edata.svd$v[,1:v.limit])

dim(edata.svd$v)
kNNdistplot(edata.svd$v[,1:v.limit], k = 20)
abline(h=0.8, col = "red", lty=2)
View(edata.svd$v)
raw.dbscan <- dbscan(edata.raw, eps = 147, minPts = 5)
norm.dbscan <- dbscan(edata.norm, eps = 192, minPts = 5)
svd.dbscan <- dbscan(edata.svd$v[,1:v.limit], eps = 0.8, minPts = 5)

svd.mclust <- Mclust(edata.svd$v[,1:v.limit])
Mclust()
pheno_data <- edata[,c(3,4)]
pheno_data$raw.kmeans <- raw.kmeans$cluster
pheno_data$norm.kmeans <- norm.kmeans$cluster
pheno_data$svd.kmeans <- svd.kmeans$cluster

pheno_data$norm.meanshift <- norm.meanshift$assignment
pheno_data$svd.meanshift <- svd.meanshift$assignment

pheno_data$raw.dbscan <- raw.dbscan$cluster
pheno_data$norm.dbscan <- norm.dbscan$cluster
pheno_data$svd.dbscan <- svd.dbscan$cluster

pheno_data$svd.mclust <- svd.mclust$classification

cluster_path <- paste("Expression_split_files/Data_analysis/Same_genes/clustering/",cancer[1],"_",source[1],".csv",sep = "")
cluster_path
write.csv(x = pheno_data, file = cluster_path)

#Read summary tables
colnames(kidney.CCLE)
kidney.CCLE <- read.csv("Expression_split_files/Data_analysis/Same_genes/clustering/kidney_CCLE.csv")
kidney.COSMIC <- read.csv("Expression_split_files/Data_analysis/Same_genes/clustering/kidney_COSMIC.csv")
kidney.NCI60 <- read.csv("Expression_split_files/Data_analysis/Same_genes/clustering/kidney_NCI60.csv")

kidney.CCLE <- kidney.CCLE[,-c(2)]
kidney.COSMIC <- kidney.COSMIC[,-c(2)]
kidney.NCI60 <- kidney.NCI60[,-c(2)]

pancreas.CCLE <- read.csv("Expression_split_files/Data_analysis/Same_genes/clustering/pancreas_CCLE.csv")
pancreas.COSMIC <- read.csv("Expression_split_files/Data_analysis/Same_genes/clustering/pancreas_COSMIC.csv")

pancreas.CCLE <- pancreas.CCLE[,-c(2)]
pancreas.COSMIC <- pancreas.COSMIC[,-c(2)]


thyroid.CCLE <- read.csv("Expression_split_files/Data_analysis/Same_genes/clustering/thyroid_CCLE.csv")
thyroid.COSMIC <- read.csv("Expression_split_files/Data_analysis/Same_genes/clustering/thyroid_COSMIC.csv")

thyroid.CCLE <- thyroid.CCLE[,-c(2)]
thyroid.COSMIC <- thyroid.COSMIC[,-c(2)]

lung.CCLE <- read.csv("Expression_split_files/Data_analysis/Same_genes/clustering/lung_CCLE.csv")
lung.COSMIC <- read.csv("Expression_split_files/Data_analysis/Same_genes/clustering/lung_COSMIC.csv")
lung.NCI60 <- read.csv("Expression_split_files/Data_analysis/Same_genes/clustering/lung_NCI60.csv")

lung.CCLE <- lung.CCLE[,-c(2)]
lung.COSMIC <- lung.COSMIC[,-c(2)]
lung.NCI60 <- lung.NCI60[,-c(2)]

lung.cancers <- unique(lung.CCLE$Diseases_extracted)
freq.lung.CCLE <- as.matrix(table(lung.CCLE$Diseases_extracted, lung.CCLE$norm.kmeans))

freq.lung.COSMIC <- as.matrix(table(lung.COSMIC$Diseases_extracted, lung.COSMIC$norm.kmeans))
freq.lung.COSMIC <- as.matrix(table(lung.COSMIC$Diseases_extracted, lung.COSMIC$norm.kmeans))
freq.lung.NCI60 <- as.matrix(table(lung.NCI60$Diseases_extracted, lung.NCI60$norm.kmeans))

freq.lung.CCLE
freq.lung.COSMIC
freq.lung.NCI60

freq.lung.CCLE.svd <- as.matrix(table(lung.CCLE$Diseases_extracted, lung.CCLE$svd.kmeans))
freq.lung.COSMIC.svd <- as.matrix(table(lung.COSMIC$Diseases_extracted, lung.COSMIC$svd.kmeans))
freq.lung.NCI60.svd <- as.matrix(table(lung.NCI60$Diseases_extracted, lung.NCI60$svd.kmeans))

freq.lung.CCLE.svd
freq.lung.COSMIC.svd
freq.lung.NCI60.svd






