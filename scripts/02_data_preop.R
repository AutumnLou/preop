###########################################################################
# PhD
# 02_data_preop.R
###########################################################################

# Read in the required packages ------------------------------------------
source("scripts/01_data_clean.R")

# Remove intermediary datasets
rm(final_data, notherapy_data)

# Create new variables
preop_data <- preop_data %>%
  mutate(ClinT_Stage = ifelse(ClinT_Stage == 'T1C',
                             'T1',
                             ifelse(ClinT_Stage == 'T2' |
                                      ClinT_Stage == 'T2A' |
                                      ClinT_Stage == 'T2B'|
                                      ClinT_Stage == 'T2C',
                                    'T2', 'T3'))) %>%
  select(!c(BxGG1, BxGG2))
preop_data$ClinT_Stage <- as.factor(preop_data$ClinT_Stage)

# Create new variables
preop_notherapy_data <- preop_notherapy_data %>%
  mutate(ClinT_Stage = ifelse(ClinT_Stage == 'T1C',
                              'T1',
                              ifelse(ClinT_Stage == 'T2' |
                                       ClinT_Stage == 'T2A' |
                                       ClinT_Stage == 'T2B'|
                                       ClinT_Stage == 'T2C',
                                     'T2', 'T3'))) %>%
  select(!c(BxGG1, BxGG2))
preop_notherapy_data$ClinT_Stage <- as.factor(preop_notherapy_data$ClinT_Stage)
