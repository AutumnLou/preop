###########################################################################

# 01_data_clean.R
###########################################################################

# Read in the required packages ------------------------------------------

source("scripts/00_packages.R")

# Load the data
data <- read.csv("data/raw_data.csv", header = TRUE)

# Attach the data
attach(data)

# Make the row names the Sample ID
rownames(data) <- data$Sample.ID
data <- data %>% select(!Sample.ID)

# Reduce number of variables
bcr_data <- data %>%
  # Remove all rows without BCR_Free Time
  filter(!is.na(BCR_FreeTime) & !is.na(BCR_Event)) %>%
  
  # Remove MetSite and Type
  select(!c(MetSite, Type)) %>%
  
  # Remove ERG columns
  select(!c(ERG.fusion.aCGH, ERG.fusion.gex)) %>%
  
  # Removed Num_Nodes_Removed and Num_Nodes_Positive
  select(!c(Num_Nodes_Removed, Num_Nodes_Positive)) %>%
  
  # Removed Nomogram and Copy Number Cluster Columns
  select(!c(Nomogram.NomoPred_ECE, Nomogram.NomoPred_LNI,
            Nomogram.NomoPred_OCD, Nomogram.NomoPred_SVI,
            Nomogram.PFP_PostRP, Copy.number.Cluster)) %>%
  
  # Remove RP type
  select(!c(RP_Type)) %>%
  
  # Remove MetsEvent, SurvTime and Event
  ## as they are post BCR or censored
  select(!c(MetsEvent, SurvTime, Event))

# New dataset for removal of therapy patients
notherapy_data <- bcr_data %>%
  filter(is.na(NeoAdjRadTx) &
           (is.na(ChemoTx) | ChemoTx == "PostCHEMO") &
           (is.na(HormTx) | HormTx == "PostHORM"))

# Remove therapy columns
bcr_data <- bcr_data %>%
  select(!c(NeoAdjRadTx, ChemoTx, HormTx, RadTxType))
notherapy_data <- notherapy_data %>%
  select(!c(NeoAdjRadTx, ChemoTx, HormTx, RadTxType))

# Create Pre-op data
preop_data <- bcr_data %>%
  select(!c(PreTxPSA, SMS, ECE, SVI, LNI, PathStage,
            PathGG1, PathGG2, PathGGS))
preop_notherapy_data <- notherapy_data %>%
  select(!c(PreTxPSA, SMS, ECE, SVI, LNI, PathStage,
            PathGG1, PathGG2, PathGGS))

# Remove incomplete observations
preop_data <- na.omit(preop_data)
bcr_comp_data <- na.omit(bcr_data)
notherapy_data <- na.omit(notherapy_data)
preop_notherapy_data <- na.omit(preop_notherapy_data)

# Correct variable types
preop_data$ClinT_Stage <- as.factor(preop_data$ClinT_Stage)
bcr_comp_data$ClinT_Stage <- as.factor(bcr_comp_data$ClinT_Stage)
bcr_comp_data$SMS <- as.factor(bcr_comp_data$SMS)
bcr_comp_data$SVI <- as.factor(bcr_comp_data$SVI)
bcr_comp_data$LNI <- as.factor(bcr_comp_data$LNI)
bcr_comp_data$PathStage <- as.factor(bcr_comp_data$PathStage)
bcr_comp_data$PathGG1 <- as.factor(bcr_comp_data$PathGG1)
bcr_comp_data$PathGG2 <- as.factor(bcr_comp_data$PathGG2)
bcr_comp_data$PathGGS <- as.factor(bcr_comp_data$PathGGS)

preop_notherapy_data$ClinT_Stage <- as.factor(preop_notherapy_data$ClinT_Stage)
notherapy_data$ClinT_Stage <- as.factor(notherapy_data$ClinT_Stage)
notherapy_data$SMS <- as.factor(notherapy_data$SMS)
notherapy_data$SVI <- as.factor(notherapy_data$SVI)
notherapy_data$LNI <- as.factor(notherapy_data$LNI)
notherapy_data$PathStage <- as.factor(notherapy_data$PathStage)
notherapy_data$PathGG1 <- as.factor(notherapy_data$PathGG1)
notherapy_data$PathGG2 <- as.factor(notherapy_data$PathGG2)
notherapy_data$PathGGS <- as.factor(notherapy_data$PathGGS)

# Change the categories of ECE column
bcr_comp_data$ECE <- ifelse(bcr_comp_data$ECE == "NONE",
                            "Absent","Present")
bcr_comp_data$ECE <- as.factor(bcr_comp_data$ECE)
notherapy_data$ECE <- ifelse(notherapy_data$ECE == "NONE",
                            "Absent","Present")
notherapy_data$ECE <- as.factor(notherapy_data$ECE)


# Change the categories of the Race column
preop_data$Race <-
  mapvalues(preop_data$Race, from = c("Black Non Hispanic", "White Non Hispanic",
                                         "Black Hispanic", "White Hispanic", "Asian",
                                         "Unknown"),
            to = c("2Black", "1White", "2Black", "1White", "3Asian", "4Unknown"))
preop_data$Race <- as.factor(preop_data$Race)
bcr_comp_data$Race <-
  mapvalues(bcr_comp_data$Race, from = c("Black Non Hispanic", "White Non Hispanic",
                                    "Black Hispanic", "White Hispanic", "Asian",
                                    "Unknown"),
            to = c("2Black", "1White", "2Black", "1White", "3Asian", "4Unknown"))
bcr_comp_data$Race <- as.factor(bcr_comp_data$Race)
notherapy_data$Race <- mapvalues(notherapy_data$Race,
                                 from = c("Black Non Hispanic",
                                          "White Non Hispanic",
                                          "Black Hispanic",
                                          "White Hispanic",
                                          "Asian", "Unknown"),
                                 to = c("2Black", "1White", "2Black", "1White",
                                        "3Asian", "4Unknown"))
notherapy_data$Race <- as.factor(notherapy_data$Race)
preop_notherapy_data$Race <- mapvalues(preop_notherapy_data$Race,
                                 from = c("Black Non Hispanic",
                                          "White Non Hispanic",
                                          "Black Hispanic",
                                          "White Hispanic",
                                          "Asian", "Unknown"),
                                 to = c("2Black", "1White", "2Black", "1White",
                                        "3Asian", "4Unknown"))
preop_notherapy_data$Race <- as.factor(preop_notherapy_data$Race)

# Change levels of Biopsy Gleason score lowest == 3
## and sum == 6 
preop_data$BxGG1 <- mapvalues(preop_data$BxGG1,
                                 from = c(2, 3, 4, 5),
                                 to = c(3, 3, 4, 5))
preop_data$BxGG2 <- mapvalues(preop_data$BxGG2,
                                 from = c(2, 3, 4, 5),
                                 to = c(3, 3, 4, 5))
preop_data$BxGGS <- mapvalues(preop_data$BxGGS,
                                 from = c(5, 6, 7, 8, 9),
                                 to = c(6, 6, 7, 8, 9))
preop_data$BxGG1 <- as.factor(preop_data$BxGG1)
preop_data$BxGG2 <- as.factor(preop_data$BxGG2)
preop_data$BxGGS <- as.factor(preop_data$BxGGS)

bcr_comp_data$BxGG1 <- mapvalues(bcr_comp_data$BxGG1,
                                 from = c(2, 3, 4, 5),
                                 to = c(3, 3, 4, 5))
bcr_comp_data$BxGG2 <- mapvalues(bcr_comp_data$BxGG2,
                                 from = c(2, 3, 4, 5),
                                 to = c(3, 3, 4, 5))
bcr_comp_data$BxGGS <- mapvalues(bcr_comp_data$BxGGS,
                                 from = c(5, 6, 7, 8, 9),
                                 to = c(6, 6, 7, 8, 9))
bcr_comp_data$BxGG1 <- as.factor(bcr_comp_data$BxGG1)
bcr_comp_data$BxGG2 <- as.factor(bcr_comp_data$BxGG2)
bcr_comp_data$BxGGS <- as.factor(bcr_comp_data$BxGGS)

notherapy_data$BxGG1 <- mapvalues(notherapy_data$BxGG1,
                                 from = c(2, 3, 4, 5),
                                 to = c(3, 3, 4, 5))
notherapy_data$BxGG2 <- mapvalues(notherapy_data$BxGG2,
                                 from = c(2, 3, 4, 5),
                                 to = c(3, 3, 4, 5))
notherapy_data$BxGGS <- mapvalues(notherapy_data$BxGGS,
                                 from = c(5, 6, 7, 8, 9),
                                 to = c(6, 6, 7, 8, 9))
notherapy_data$BxGG1 <- as.factor(notherapy_data$BxGG1)
notherapy_data$BxGG2 <- as.factor(notherapy_data$BxGG2)
notherapy_data$BxGGS <- as.factor(notherapy_data$BxGGS)

preop_notherapy_data$BxGG1 <- mapvalues(preop_notherapy_data$BxGG1,
                              from = c(2, 3, 4, 5),
                              to = c(3, 3, 4, 5))
preop_notherapy_data$BxGG2 <- mapvalues(preop_notherapy_data$BxGG2,
                              from = c(2, 3, 4, 5),
                              to = c(3, 3, 4, 5))
preop_notherapy_data$BxGGS <- mapvalues(preop_notherapy_data$BxGGS,
                              from = c(5, 6, 7, 8, 9),
                              to = c(6, 6, 7, 8, 9))
preop_notherapy_data$BxGG1 <- as.factor(preop_notherapy_data$BxGG1)
preop_notherapy_data$BxGG2 <- as.factor(preop_notherapy_data$BxGG2)
preop_notherapy_data$BxGGS <- as.factor(preop_notherapy_data$BxGGS)

# Change the BCR event to a numeric variable
preop_data$BCR_Event <- ifelse(preop_data$BCR_Event == "NO",
                                  0, 1)
bcr_comp_data$BCR_Event <- ifelse(bcr_comp_data$BCR_Event == "NO",
                                  0, 1)
notherapy_data$BCR_Event <- ifelse(notherapy_data$BCR_Event == "NO",
                                  0, 1)
preop_notherapy_data$BCR_Event <- ifelse(preop_notherapy_data$BCR_Event == "NO",
                                   0, 1)

## Drop Observation with PSA 506 (error) and Drop MET duplicate
preop_data <- preop_data[!(row.names(preop_data)
                              %in% c('PCA0045','PCA0187', 'PCA0207')), ]
final_data <- bcr_comp_data[!(row.names(bcr_comp_data)
                                 %in% c('PCA0045','PCA0187', 'PCA0207')), ]
notherapy_data <- notherapy_data[!(row.names(notherapy_data)
                              %in% c('PCA0045','PCA0187', 'PCA0207')), ]
preop_notherapy_data <- preop_notherapy_data[!(row.names(preop_notherapy_data)
                                   %in% c('PCA0045','PCA0187', 'PCA0207')), ]

# Remove intermediary datasets
rm(data, bcr_data, bcr_comp_data)

# Make a .csv file
write.csv(final_data, 'data/clean_bcr_data.csv')
write.csv(notherapy_data, 'data/notherapy_data.csv')
write.csv(preop_data, 'data/clean_preop_data.csv')
write.csv(preop_notherapy_data, 'data/preop_notherapy_data.csv')
