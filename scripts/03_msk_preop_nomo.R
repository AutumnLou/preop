###########################################################################
# MSKCC
# 03_msk_preop_nomo.R
###########################################################################

# Load the dataset
source("scripts/02_data_msk_preop.R")
attach(msk_wide_data)

# Load all coefficients
intercept <- 6.52947562
age <- -0.00296141
psa <- -0.39711402
psa_s1 <- 0.0030316
psa_s2 <- -0.00828079
bgg2 <- -1.07031938
bgg3 <- -2.12762264
bgg4 <- -2.62959495
bgg5 <- -3.00134586
cs_2a <- -0.45324021
cs_2b <- -0.88923348
cs_2c <- -0.89903664
cs_3 <- -1.31807837

# Load scaling parameter
scaling_param <- 1.11186766

# c_index <- 0.79051447
# model_no <- 13640

# Load Spline Knots
PSAPreopKnot1 <- 0.2
PSAPreopKnot2 <- 4.71
PSAPreopKnot3 <- 7.22
PSAPreopKnot4 <- 96.53

# Last Updated: 14 April 2021 

# Create the spline variables
msk_wide_spline_data <- msk_wide_data %>%
  mutate(psa_k1 = PreDxBxPSA - PSAPreopKnot1,
         psa_k2 = PreDxBxPSA - PSAPreopKnot2,
         psa_k3 = PreDxBxPSA - PSAPreopKnot3,
         psa_k4 = PreDxBxPSA - PSAPreopKnot4) %>%
  mutate(sp1var = pmax(psa_k1, 0)**3 -
           (pmax(psa_k3,0)**3)*(PSAPreopKnot4 - PSAPreopKnot1)/(PSAPreopKnot4 - PSAPreopKnot3) +
           (pmax(psa_k4,0)**3)*(PSAPreopKnot3 - PSAPreopKnot1)/(PSAPreopKnot4 - PSAPreopKnot3),
         sp2var = (pmax((psa_k2),0))**3 -
           ((pmax((psa_k3),0))**3)*((PSAPreopKnot4 - PSAPreopKnot2)/(PSAPreopKnot4 - PSAPreopKnot3)) +
           ((pmax((psa_k4),0))**3)*((PSAPreopKnot3 - PSAPreopKnot2)/(PSAPreopKnot4 - PSAPreopKnot3))
  )

# Reattach data with additional variables
attach(msk_wide_spline_data)

# Build the survival model
surv_model <- intercept + age*DxAge + psa*PreDxBxPSA + psa_s1*sp1var +
  psa_s2*sp2var + bgg2*biop_grp_2 + bgg3*biop_grp_3 + bgg4*biop_grp_4 +
  bgg5*biop_grp_5 + cs_2a*stage_2A + cs_2b*stage_2B + cs_2c*stage_2C +
  cs_3*stage_3p

# Calculate the prediction
pred_prob_5 <- 1/(1 + (exp(-surv_model)*5)**(1/scaling_param))

# Add predictions to dataset
msk_pred_data <- msk_wide_data %>%
  select(BCR_FreeTime, BCR_Event)
msk_pred_data$surv5 <- pred_prob_5
msk_pred_data$risk <- surv_model


# c-index
c_index <- survcomp::concordance.index(x = 1 - msk_pred_data$surv5,
                             surv.time = msk_pred_data$BCR_FreeTime,
                             surv.event = msk_pred_data$BCR_Event)
c_index$c.index
## 0.7245467

# Make a .csv file
write.csv(msk_pred_data, 'data/msk_preop_pred_data.csv', row.names = FALSE)

