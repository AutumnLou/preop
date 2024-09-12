###########################################################################
# PhD
# 02_data_comb_preop.R
###########################################################################

# Load the dataset
source("scripts/02_data_preop.R")
source("scripts/02_data_rna.R")

#Combine the Clinical and mRNA data
comb_data <- merge(preop_data, rna_data, by = 0)
rownames(comb_data) <- comb_data[, 1]
comb_data <- comb_data[, -1]
comb_data[,8:ncol(comb_data)] <- lapply(comb_data[,8:ncol(comb_data)], function(x) {
  if(is.character(x)) as.numeric(x)
})

comb_notherapy_data <- merge(preop_notherapy_data, rna_data, by = 0)
rownames(comb_notherapy_data) <- comb_notherapy_data[, 1]
comb_notherapy_data <- comb_notherapy_data[, -1]
comb_notherapy_data[,8:ncol(comb_notherapy_data)] <- lapply(comb_notherapy_data[,8:ncol(comb_notherapy_data)], function(x) {
  if(is.character(x)) as.numeric(x)
})

# Remove intermediary datasets
rm(rna_data, preop_data)

