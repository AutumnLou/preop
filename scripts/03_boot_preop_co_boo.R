###########################################################################
# PhD - Pre Op Journal Article
## Boosted Bootstrapped
# 03_boot_preop_co_boo.R
###########################################################################

# Load the dataset
source("scripts/02_data_preop.R")

# Load functions/libraries
source("scripts/mycalibration_function2.R")
source("scripts/mysurvFit_functions.R")
require(mboost)

# Boost
## Data preparation
x_full <- model.matrix( ~ .-1, preop_data[,1:5])
y_full <- Surv(time = preop_data$BCR_FreeTime,
               event = preop_data$BCR_Event)

# Model w/ all 190 obs (app)
boost_cox_app = glmboost(x = x_full, y = y_full,
                     family = CoxPH())
pred_app <- predict.glmboost(boost_cox_app, newdata = x_full)

# Discrimination
c_app <- 1-compareC::estC(preop_data$BCR_FreeTime, preop_data$BCR_Event, pred_app)


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


stab.v <- rep.int(0,length(colnames(x_full)))
names(stab.v) <- colnames(x_full)

for(i in 1:length(bootIndex)) {
  print(i)
  # define Training / Validation data
  bag <- preop_data[bootIndex[[i]],]
  oob <- preop_data[-bootIndex[[i]],]
  
  # Boost data prep
  x <- model.matrix( ~ .-1, bag[,1:5])
  y <- Surv(time = bag$BCR_FreeTime,
            event = bag$BCR_Event)
  test_x <- model.matrix( ~ .-1, oob[,1:5])
  test_y <- Surv(time = oob$BCR_FreeTime,
                 event = oob$BCR_Event)
  
  # Boost Model
  boost_cox = glmboost(x = x, y = y, family = CoxPH())
  train_pred <- predict.glmboost(boost_cox, newdata = x)
  
  # Stability
  terms <- names(boost_cox$coef()[1:length(boost_cox$coef())])
  for(j in 1:length(terms)) {
    stab.v[which(names(stab.v) == terms[j])] <- stab.v[which(names(stab.v) == terms[j])] + 1
  }
  
  # Validation
  oob_val <- predict.glmboost(boost_cox, newdata = test_x)
  orig_val <- predict.glmboost(boost_cox, newdata = x_full)
  
  ## OOB
  oob_surv <- survFit(boost_cox, test_x)
  oob_surv_time <- which(abs(oob_surv$time-60)==min(abs(oob_surv$time-60)))
  oob_surv5yr <- oob_surv$surv[oob_surv_time,]
  oob_cal_dat <- data.frame(BCR_FreeTime = oob$BCR_FreeTime,
                            BCR_Event = oob$BCR_Event,
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
  oob_km <- survfit(Surv(oob$BCR_FreeTime,oob$BCR_Event)~1)
  oob_Surv.output = Surv(oob$BCR_FreeTime,oob$BCR_Event)
  oob_sc = aggregate.survFit(survFit(boost_cox, test_x))
  oob_cal <- as.data.frame(my.calibration2(oob_sc, oob_km, oob_Surv.output))
  oob_cal <- oob_cal %>% mutate(group = rownames(oob_cal))
  
  oob_res[(((i-1)*5)+1):(i*5),] <- oob_cal
  
  # store results
  c_oob[i] <- 1-compareC::estC(oob$BCR_FreeTime, oob$BCR_Event, oob_val)
  o_boot[i] <- (1-compareC::estC(bag$BCR_FreeTime, bag$BCR_Event, train_pred)) - (1-compareC::estC(preop_data$BCR_FreeTime, preop_data$BCR_Event, orig_val))
  
}


# Stability
app_var <- names(boost_cox_app$coef()[1:length(boost_cox_app$coef())])
stab.v <- sort(stab.v, decreasing = TRUE)
barchart(head(stab.v, 10),
         xlab = "Number of uses accross resamples",
         xlim = 0:100)

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

