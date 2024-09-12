###########################################################################
# PhD - Pre Op Journal Article
## RSF Bootstrapped
# 03_boot_preop_co_rsf.R
###########################################################################

# Load the dataset
source("scripts/02_data_preop.R")
source("scripts/mycalibration_function2.R")
source("scripts/mypredictSurvProb_functions.R")
source("scripts/RFE_for_RSF_function.R")

# Prep date for rfe function
colnames(preop_data)[6] <- "time"
colnames(preop_data)[7] <- "status"

# Model w/ all 190 obs (app)
set.seed(123)
rsf_app <- rsf_rfe(preop_data)
pred_app <- predict.rfsrc(rsf_app)
round(rsf_app$importance, 3)
# Discrimination
c_app <- 1 - pred_app$err.rate[500]

###########################################################################
# Resampling Framework
###########################################################################
# 100 Bootstrap resamples
set.seed(123)
B <- 100
bootIndex <- createResample(factor(preop_data$status), times = B, list = TRUE)

c_oob <- numeric(B)
o_boot <- numeric(B)
oob_har_res <- data.frame(group = numeric(B*5), mean.p = numeric(B*5),
                          km = numeric(B*5), lower = numeric(B*5),
                          upper = numeric(B*5))

oob_res <- data.frame(pred = numeric(B*5), obs = numeric(B*5),
                      mean.obs = numeric(B*5), low = numeric(B*5),
                      up = numeric(B*5), means.low = numeric(B*5),
                      means.up = numeric(B*5), group = numeric(B*5))

stab.v <- rep.int(0,length(colnames(preop_data[,1:5])))
names(stab.v) <- colnames(preop_data[,1:5])
selection.v <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(selection.v) <- colnames(preop_data[,1:5])

for(i in 1:length(bootIndex)) {
  print(i)
  # define Training / Validation data
  bag <- preop_data[bootIndex[[i]],]
  oob <- preop_data[-bootIndex[[i]],]
  
  # RSF Selection
  set.seed(123)
  rsf_boot <- rsf_rfe(bag)
  pred_boot <- predict.rfsrc(rsf_boot)
  
  # Stability
  terms <- names(rsf_boot$importance)
  for(j in 1:length(terms)) {
    stab.v[which(names(stab.v) == terms[j])] <- stab.v[which(names(stab.v) == terms[j])] + 1
  }
  for(j in 1:length(terms)) {
    selection.v[i,which(colnames(selection.v) == terms[j])] <- 1
  }
  
  # Validation
  oob_val <- predict.rfsrc(rsf_boot, newdata = oob)
  orig_val <- predict.rfsrc(rsf_boot, newdata = preop_data)
  
  ## OOB
  oob_surv5yr <- pec::predictSurvProb(rsf_boot, newdata = oob, times = 60)
  oob_cal_dat <- data.frame(BCR_FreeTime = oob$time,
                            BCR_Event = oob$status,
                            surv5yr = oob_surv5yr,
                            row.names = rownames(oob))
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
  oob_km <- survfit(Surv(oob$time,oob$status)~1)
  oob_Surv.output = Surv(oob$time,oob$status)
  oob_sc = aggregate.predictSurvProb(pec::predictSurvProb(rsf_boot, newdata = oob, times = oob_km$time), oob_km$time)
  oob_cal <- as.data.frame(my.calibration2(oob_sc, oob_km, oob_Surv.output))
  oob_cal <- oob_cal %>% mutate(group = rownames(oob_cal))
  
  oob_res[(((i-1)*5)+1):(i*5),] <- oob_cal
  
  # store results
  c_oob[i] <- 1 - oob_val$err.rate[500]
  o_boot[i] <- (1 - pred_boot$err.rate[500]) - (1 - orig_val$err.rate[500])
}

###########################################################################
# Results
###########################################################################

# Stability
app_var <- names(rsf_app$importance)
stab.v <- sort(stab.v, decreasing = TRUE)
barchart(head(stab.v, 10),
         xlab = "Number of uses across resamples",
         xlim = 0:100)

# Selection
selection.v <- ifelse(is.na(selection.v), 0, 1)
selection.v <- as.data.frame(selection.v)

# OOB Discrimination
mean(c_oob)
quantile(c_oob, probs = c(0.025, 0.975))

# Harrels Discrimination
o <- mean(o_boot)
c_app - o
quantile(c_app - o_boot, probs = c(0.025, 0.975))

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

# save discrimination data
dis_dat <- data.frame(oob = c_oob,
                      har = c_app - o_boot)