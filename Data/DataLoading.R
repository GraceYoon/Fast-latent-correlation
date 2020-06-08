rm(list=ls())
###Load the data.  Each file can be downloaded from the TCGA Data Portal at link "https://gdc.cancer.gov/about-data/publications/brca_2012/" (accessed 06/07/2020)

library(data.table)
rna_file_name <- "BRCA.exp.547.med.txt"
GE <- fread(paste("Data/", rna_file_name, sep='')) # 547
# Accessed (04/17/2018), dim(GE) # 17814 by 548
###### To add Gene names to rownames of GE data. Changing the format to data.frame can lose the rownames. Keep it.
GErownames <- GE$NAME

mirna_file_name <- "BRCA.780.precursor.txt"
miRNA <- fread(paste("Data/", mirna_file_name, sep='')) #780
# Accessed (04/17/2018), dim(miRNA) # 1046 by 781
##### To add miRNA names to rownames of miRNA data. Changing the format to data.frame can lose the rownames. Keep it.
miRNArownames <- miRNA$Gene

# Do matching between different subjects/samples.
# Important: get rid of normal/metastatic, use only primary tumor (01).
# TCGA code Tables can be foudn at link https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables (accessed 04/19/2018).
# For GE
subset_cancer <- (1:ncol(GE))[grep("-01.-", colnames(GE))] 
length(subset_cancer) # agrees with the numbers of PT posted on the website # 522
GE <- data.frame(GE)
GE <- GE[, subset_cancer]
# For MiRNA
subset_cancer <- (1:ncol(miRNA))[grep("-01.-", colnames(miRNA))] 
length(subset_cancer) # 694
miRNA <- data.frame(miRNA)
miRNA <- miRNA[, subset_cancer]

subjectGE <- substr(colnames(GE),1,12) # There is no subject repetition. only 12 letters are considered to distinguish subjects, for example, "TCGA.A1.A0SD".
subjectmiRNA <- substr(colnames(miRNA), 1, 12)
subject_both = intersect(subjectGE, subjectmiRNA)
length(subject_both) # 500 subjects overlap

# Select only the ones with overlapping subjects
GE <- GE[, subjectGE %in% subject_both] # dim(GE) # 17814 by 500
miRNA <- miRNA[, subjectmiRNA %in% subject_both] # dim(miRNA) # 1046 by 500

GE.mat = as.matrix(GE); rownames(GE.mat) <- GErownames
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("impute")

GE.mat = impute::impute.knn(GE.mat) ##Impute missing values via KNN (K=10)
GE.mat = GE.mat$data
miRNA.mat = as.matrix(miRNA); rownames(miRNA.mat) <- miRNArownames

# Nota that miRNA has tons and tons of zeroes, and highly skewed
hist(miRNA.mat, breaks = 100)

# Filter miRNA only the features with less that a certain percent of zero values
threshold = 0.5
hist(rowSums(miRNA.mat==0), breaks = 100)
miRNA.mat = miRNA.mat[rowSums(miRNA.mat==0) < ncol(miRNA.mat)*threshold,] # dim(miRNA.mat) # 431 by 500

# Filter GE to select only the features that have a large enough standard deviation
GEsd = apply(GE.mat,1,sd)
GE.mat = GE.mat[GEsd>quantile(GEsd, 0.95),] # dim(GE.mat) # 891 by 500

# Match the order of the subjects for GE and miRNA data.
sum(substr(colnames(GE.mat), 1, 12) == subject_both)
ind.miRNA <- match(subject_both, substr(colnames(miRNA.mat), 1, 12))

miRNA.mat.subjectmatched <- miRNA.mat[, ind.miRNA]
sum(substr(colnames(GE.mat), 1, 12) == substr(colnames(miRNA.mat.subjectmatched), 1, 12))

# Combine the data and subjects
matchedTCGA <- list(X1 = t(GE.mat), X2 = t(miRNA.mat.subjectmatched), subject = subject_both)
save(matchedTCGA, file = "Data/matchedTCGA.Rdata")

# dim(matchedTCGA$X1) = 500 by 891 and dim(matchedTCGA$X2) = 500 by 431. (Checked on 06/07/2020)
# Final check
# sum(substr(rownames(matchedTCGA$X1), 1, 12) == substr(rownames(matchedTCGA$X2), 1, 12))
# sum(subject_both == substr(rownames(matchedTCGA$X2), 1, 12))
# sum(substr(rownames(matchedTCGA$X1), 1, 12) == subject_both)


