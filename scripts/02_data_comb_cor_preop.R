###########################################################################
# PhD
# 02_data_comb_cor_preop.R
###########################################################################

# Load the dataset
source("scripts/02_data_comb_preop.R")

# Pre-filter Correlation between mRNA variables
cor1_df <- cor(comb_data[8:8823], comb_data[8:8823])
fc1 <- findCorrelation(cor1_df, cutoff = 0.6, names = TRUE)
cor2_df <- cor(comb_data[8:8823], comb_data[8824:17638])
fc2 <- findCorrelation(cor2_df, cutoff = 0.6, names = TRUE)
cor3_df <- cor(comb_data[8:8823], comb_data[17639:26454])
fc3 <- findCorrelation(cor3_df, cutoff = 0.6, names = TRUE)
cor4_df <- cor(comb_data[8824:17638], comb_data[8824:17638])
fc4 <- findCorrelation(cor4_df, cutoff = 0.6, names = TRUE)
cor5_df <- cor(comb_data[8824:17638], comb_data[17639:26454])
fc5 <- findCorrelation(cor5_df, cutoff = 0.6, names = TRUE)
cor6_df <- cor(comb_data[17639:26454], comb_data[17639:26454])
fc6 <- findCorrelation(cor6_df, cutoff = 0.6, names = TRUE)

fc <- c(fc1, fc2, fc3, fc4, fc5, fc6)
fc <- fc[!duplicated(fc)]

cor_df <- comb_data %>%
  select(!(fc))

# Remove intermediary datasets
rm(comb_data, cor1_df, cor2_df, cor3_df, cor4_df, cor5_df,
   cor6_df, fc1, fc2, fc3, fc4, fc5, fc6, fc)


