###########################################################################
# PhD
# 02_data_comb_preop_matrix.R
###########################################################################

# Load the dataset
source("scripts/02_data_preop.R")
source("scripts/02_data_rna.R")

# Prep for matrix format
mat_clin <- model.matrix(~ .-1, preop_data[,1:5])
mat_clin <- as.data.frame(mat_clin)
mat_clin$BCR_FreeTime <- preop_data$BCR_FreeTime
mat_clin$BCR_Event <- preop_data$BCR_Event


#Combine the Clinical and mRNA data
comb_data <- merge(mat_clin, rna_data, by = 0)
rownames(comb_data ) <- comb_data [, 1]
comb_data <- comb_data [, -1]
comb_data[,14:26460] <- lapply(comb_data[,14:26460], function(x) {
  if(is.character(x)) as.numeric(x)
})

# Remove intermediary datasets
rm(rna_data, preop_data, mat_clin)