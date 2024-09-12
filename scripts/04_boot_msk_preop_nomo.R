###########################################################################
# MSKCC
## Bootstrapped
# 04_boot_msk_preop_nomo.R
###########################################################################

source("scripts/03_msk_preop_nomo.R")
source("scripts/mycalibration_function2.R")
attach(msk_wide_data)

# Model w/ all 187 obs (app)
## Discrimination
c_app <- c_index$c.index

# 100 Bootstrap resamples
set.seed(123)
B <- 100
bootIndex <- createResample(factor(msk_wide_data$BCR_Event), times = B, list = TRUE)

c_oob <- numeric(B)
c_boot <- numeric(B)

oob_har_res <- data.frame(group = numeric(B*5), mean.p = numeric(B*5),
                          km = numeric(B*5), lower = numeric(B*5),
                          upper = numeric(B*5))

oob_res <- data.frame(pred = numeric(B*5), obs = numeric(B*5),
                      mean.obs = numeric(B*5), low = numeric(B*5),
                      up = numeric(B*5), means.low = numeric(B*5),
                      means.up = numeric(B*5), group = numeric(B*5))

for(i in 1:length(bootIndex)) {
  print(i)
  # define Training / Validation data
  bag <- msk_wide_data[bootIndex[[i]],]
  oob <- msk_wide_data[-bootIndex[[i]],]
  
  # Create the spline variables (bag)
  bag <- bag %>%
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
  
  # Build the survival model (bag)
  bag_surv_model <- intercept + age*bag$DxAge + psa*bag$PreDxBxPSA +
    psa_s1*bag$sp1var + psa_s2*bag$sp2var + bgg2*bag$biop_grp_2 +
    bgg3*bag$biop_grp_3 + bgg4*bag$biop_grp_4 + bgg5*bag$biop_grp_5 +
    cs_2a*bag$stage_2A + cs_2b*bag$stage_2B + cs_2c*bag$stage_2C +
    cs_3*bag$stage_3p
  
  # Calculate the prediction (bag)
  bag_pred_prob_5 <- 1/(1 + (exp(-bag_surv_model)*5)**(1/scaling_param))
  
  # Add predictions to dataset (bag)
  bag_pred_data <- bag %>%
    select(BCR_FreeTime, BCR_Event)
  bag_pred_data$surv5 <- bag_pred_prob_5
  bag_pred_data$risk <- bag_surv_model
  
  # Create the spline variables (bag)
  oob <- oob %>%
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
  
  # Build the survival model (oob)
  oob_surv_model <- intercept + age*oob$DxAge + psa*oob$PreDxBxPSA +
    psa_s1*oob$sp1var + psa_s2*oob$sp2var + bgg2*oob$biop_grp_2 +
    bgg3*oob$biop_grp_3 + bgg4*oob$biop_grp_4 + bgg5*oob$biop_grp_5 +
    cs_2a*oob$stage_2A + cs_2b*oob$stage_2B + cs_2c*oob$stage_2C +
    cs_3*oob$stage_3p
  
  # Calculate the prediction (oob)
  oob_pred_prob_5 <- 1/(1 + (exp(-oob_surv_model)*5)**(1/scaling_param))
  
  # Add predictions to dataset (oob)
  oob_pred_data <- oob %>%
    select(BCR_FreeTime, BCR_Event)
  oob_pred_data$surv5 <- oob_pred_prob_5
  oob_pred_data$risk <- oob_surv_model
  
  # c-index
  bag_c_index <- survcomp::concordance.index(x = 1 - bag_pred_data$surv5,
                               surv.time = bag_pred_data$BCR_FreeTime,
                               surv.event = bag_pred_data$BCR_Event)
  oob_c_index <- survcomp::concordance.index(x = 1 - oob_pred_data$surv5,
                                   surv.time = oob_pred_data$BCR_FreeTime,
                                   surv.event = oob_pred_data$BCR_Event)
  
  ## OOB
  oob_cal_dat <- data.frame(BCR_FreeTime = oob_pred_data$BCR_FreeTime,
                            BCR_Event = oob_pred_data$BCR_Event,
                            surv5yr = oob_pred_data$surv5)
  oob_cal_dat <- oob_cal_dat[order(oob_cal_dat$surv5yr),]
  oob_cal_dat$group <- cut(1:nrow(oob_cal_dat), 5, labels = FALSE)
  oob_cal_dat <- oob_cal_dat %>% 
    group_by(group) %>% 
    mutate(mean.p = mean(surv5yr),
           km = summary(survfit(Surv(BCR_FreeTime,BCR_Event)~1), times = 60, extend = TRUE)$surv,
           lower = summary(survfit(Surv(BCR_FreeTime,BCR_Event)~1), times = 60, extend = TRUE)$lower,
           upper = summary(survfit(Surv(BCR_FreeTime,BCR_Event)~1), times = 60, extend = TRUE)$upper) %>%
    summarise(mean.p = mean(mean.p),
              km = mean(km),
              lower = mean(lower),
              upper = mean(upper))
  
  oob_har_res[(((i-1)*5)+1):(i*5),] <- oob_cal_dat
  
  ## OOB
  oob_km <- survfit(Surv(oob$BCR_FreeTime,oob$BCR_Event)~1)
  oob_Surv.output = Surv(oob$BCR_FreeTime,oob$BCR_Event)
  # Calculate prediction for each event time (bag)
  oob_times <- unique(sort(oob$BCR_FreeTime))
  oob_overall_surv <- data.frame(matrix(nrow = nrow(oob), ncol = length(oob_times)))
  for (j in 1:length(oob_times)){
    oob_overall_surv[,j] <- 1/(1 + (exp(-oob_surv_model)*oob_times[j])**(1/scaling_param))
  }
  oob_sc <- data.frame(time = oob_times,
                       surv = colMeans(oob_overall_surv))
  oob_cal <- as.data.frame(my.calibration2(oob_sc, oob_km, oob_Surv.output))
  oob_cal <- oob_cal %>% mutate(group = rownames(oob_cal))
  
  oob_res[(((i-1)*5)+1):(i*5),] <- oob_cal
  
  # store results
  c_boot[i] <- bag_c_index$c.index
  c_oob[i] <- oob_c_index$c.index
}

# OOB Discrimination
mean(c_oob)
quantile(c_oob, probs = c(0.025, 0.975))

# Bootstrapped Discrimination
mean(c_boot)
quantile(c_boot, probs = c(0.025, 0.975))

# Calibration
oob_har_res <- oob_har_res %>%
  group_by(group) %>%
  summarise(mean.p = mean(mean.p),
            km = mean(km),
            lower = mean(lower),
            upper = mean(upper))

my.calibration.plot(oob_har_res$mean.p,
                    oob_har_res$km,
                    oob_har_res$lower,
                    oob_har_res$upper,
                    xlab="Predicted 5-year Survival",
                    ylab="Observed 5-year Survival")

oob_res <- oob_res %>%
  group_by(group) %>%
  summarise(pred = mean(pred, na.rm = T),
            obs = mean(obs, na.rm = T),
            mean.obs = mean(mean.obs, na.rm = T),
            low = mean(low, na.rm = T),
            up = mean(up, na.rm = T),
            means.low = mean(means.low, na.rm = T),
            means.up = mean(means.up, na.rm = T))

my.calibration.plot2(oob_res,
                     xlab="Predicted",
                     ylab="Observed")

dis_dat <- data.frame(dis = c_boot, oob = c_oob)

