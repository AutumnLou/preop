###########################################################################
# PhD - Pre Op Journal Article
## Traditional Cox Bootstrapped
# 03_boot_preop_co_cox.R
###########################################################################

# Load the dataset
source("scripts/02_data_preop.R")
source("scripts/mycalibration_function2.R")

# Model w/ all 190 obs (app)
app_start_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ 1, data = preop_data, x = TRUE)
app_full_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ ., data = preop_data, x = TRUE)
app_fit_step <- step(app_start_cox, direction = "both", scope = app_full_cox$formula, trace = 0)
cox.zph(app_fit_step)
app_fit_step$coef
# Discrimination
c_app <- concordance(app_fit_step)$concordance

# 100 Bootstrap resamples
set.seed(123)
B <- 100
bootIndex <- createResample(factor(preop_data$BCR_Event), times = B, list = TRUE)

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

###########################################################################
# Resampling Framework
###########################################################################
for(i in 1:length(bootIndex)) {
  print(i)
  # define Training / Validation data
  bag <- preop_data[bootIndex[[i]],]
  oob <- preop_data[-bootIndex[[i]],]
  
  # Stepwise Selection
  start_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ 1, data = bag, x = TRUE)
  
  full_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ ., data = bag, x = TRUE)
  
  fit_step <- step(start_cox, direction = "both", scope = full_cox$formula, trace = 0)
  
  # Stability
  terms <- attr(terms(fit_step), "term.labels")
  for(j in 1:length(terms)) {
    stab.v[which(names(stab.v) == terms[j])] <- stab.v[which(names(stab.v) == terms[j])] + 1
  }
  for(j in 1:length(terms)) {
    selection.v[i,which(colnames(selection.v) == terms[j])] <- 1
  }

  # Validation
  oob_val <- concordance(fit_step, newdata = oob)
  orig_val <- concordance(fit_step, newdata = preop_data)
  
  oob_surv <- survfit(fit_step, newdata = oob)
  oob_surv5yr <- summary(oob_surv, times = 60)
  oob_cal_dat <- data.frame(BCR_FreeTime = oob$BCR_FreeTime,
                            BCR_Event = oob$BCR_Event,
                            surv5yr = t(oob_surv5yr$surv),
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
  oob_km <- survfit(Surv(oob$BCR_FreeTime,oob$BCR_Event)~1)
  oob_Surv.output = Surv(oob$BCR_FreeTime,oob$BCR_Event)
  oob_sc = aggregate(survfit(fit_step, newdata = oob))
  oob_cal <- as.data.frame(my.calibration2(oob_sc, oob_km, oob_Surv.output))
  oob_cal <- oob_cal %>% mutate(group = rownames(oob_cal))
  
  oob_res[(((i-1)*5)+1):(i*5),] <- oob_cal
  
  # store results
  c_oob[i] <- oob_val$concordance
  o_boot[i] <- concordance(fit_step)$concordance - orig_val$concordance
}

###########################################################################
# Results
###########################################################################
# Stability
app_var <- attr(terms(app_fit_step), "term.labels")
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
