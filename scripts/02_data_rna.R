###########################################################################
# PhD
# 02_data_rna.R
###########################################################################

# Note: As file size of "MSKCC_PCa_mRNA_data.txt" is too large for Github use
# "clean_rna_data.csv" available instead

#Load mRNA Data
rna_data <- read.table("data/MSKCC_PCa_mRNA_data.txt", header = T)
rna_data <- t(rna_data)
rna_data <- rna_data[-1,]
colnames(rna_data) <- rna_data[1, ]
rna_data <- rna_data[-1, ]
colnames(rna_data)[16923] <- 'NA_var'
write.csv(rna_data, "data/clean_rna_data.csv")