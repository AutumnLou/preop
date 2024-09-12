###########################################################################
# MSKCC
# 02_data_msk_preop.R
###########################################################################

# Read in the required packages ------------------------------------------

source("scripts/01_data_clean.R")
rm(notherapy_data, final_data)

# Create a dataset for use with the MSKCC Nomogram
msk_nomo_data <- preop_data

# Remove all unused columns
msk_nomo_data <- msk_nomo_data %>%
  
  # Remove Race
  select(!Race)

# Create new variables
msk_nomo_data <- msk_nomo_data %>%
  mutate(biop_grp = ifelse(BxGGS == 6,
                           1,
                           ifelse(BxGG1 == 3 & BxGG2 == 4,
                                  2,
                                  ifelse(BxGG1 == 4 & BxGG2 == 3,
                                         3,
                                         ifelse(BxGG1 == 4 & BxGG2 == 4,
                                                4, 5)))),
         stage = mapvalues(ClinT_Stage,
                           from = c("T1C", "T2", "T2A",
                                    "T2B", "T2C", "T3",
                                    "T3A", "T3B", "T3C"),
                           to = c("1C", "2A", "2A",
                                  "2B", "2C", "3p",
                                  "3p", "3p", "3p"))
  )

# Remove the original columns
msk_nomo_data <- msk_nomo_data %>%
  select(!c(BxGG1, BxGG2, BxGGS, ClinT_Stage))

# Widen the dataset
msk_wide_data <- msk_nomo_data %>%
  dummy_cols(select_columns = c("biop_grp", "stage")) %>%
  
  # Remove original columns
  select(!c(biop_grp, stage))

# Remove intermediary datasets
rm(preop_data, msk_nomo_data)

# Make a .csv file
write.csv(msk_wide_data, 'data/msk_wide_preop_data.csv')
