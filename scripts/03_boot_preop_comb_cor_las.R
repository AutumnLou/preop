###########################################################################
# PhD - Pre Op Journal Article
## LASSO Cox Bootstrapped
## mRNA data correlation pre-filter
# 03_boot_preop_comb_cor_las.R
###########################################################################

# Load the dataset
source("scripts/02_data_comb_cor_preop.R")
source("scripts/mycalibration_function2.R")

###########################################################################
# Prep and calculation of Apparent C-Index
###########################################################################

# Lasso
## Data preparation
x_full <- model.matrix( ~ .-1, cor_df[,c(1:5,8:ncol(cor_df))])
y_full <- Surv(time = cor_df$BCR_FreeTime,
               event = cor_df$BCR_Event)

# Model w/ all 135 obs (app)
## Lasso Selection
set.seed(123)
cv_lasso_cox_app <- cv.glmnet(x = x_full,
                              y = y_full,
                              family="cox")
plot(cv_lasso_cox_app)

lasso_cox_app <- glmnet(x = x_full,
                        y = y_full,
                        family="cox",
                        lambda = cv_lasso_cox_app$lambda.1se)
pred_app <- predict(lasso_cox_app, newx = x_full)


# Discrimination
c_app <- glmnet::Cindex(pred_app, y = y_full)


###########################################################################
# Resampling Framework
###########################################################################

# 100 Bootstrap resamples
set.seed(123)
B <- 100
bootIndex <- createResample(factor(cor_df$BCR_Event), times = B, list = TRUE)

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
  bag <- cor_df[bootIndex[[i]],]
  oob <- cor_df[-bootIndex[[i]],]
  
  # Lasso data prep
  x <- model.matrix( ~ .-1, bag[,c(1:5,8:ncol(bag))])
  y <- Surv(time = bag$BCR_FreeTime,
            event = bag$BCR_Event)
  test_x <- model.matrix( ~ .-1, oob[,c(1:5,8:ncol(oob))])
  test_y <- Surv(time = oob$BCR_FreeTime,
                 event = oob$BCR_Event)
  
  # Lasso Modelling
  set.seed(123)
  cv_lasso_cox <- cv.glmnet(x = x,
                            y = y,
                            family="cox")
  lasso_cox <- glmnet(x = x,
                      y = y,
                      family="cox",
                      lambda = cv_lasso_cox$lambda.1se)
  train_pred <- predict(lasso_cox, newx = x)
  
  # Stability
  terms <- rownames(lasso_cox$beta)[lasso_cox$beta@i +1]
  for(j in 1:length(terms)) {
    stab.v[which(names(stab.v) == terms[j])] <- stab.v[which(names(stab.v) == terms[j])] + 1
  }
  
  # Validation
  oob_val <- predict(lasso_cox, newx = test_x)
  orig_val <- predict(lasso_cox, newx = x_full)
  
  ## OOB
  oob_surv <- survfit(lasso_cox, x = x, y = y, newx = test_x)
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
  oob_sc = aggregate(survfit(lasso_cox, x = x, y = y, newx = test_x))
  oob_cal <- as.data.frame(my.calibration2(oob_sc, oob_km, oob_Surv.output))
  oob_cal <- oob_cal %>% mutate(group = rownames(oob_cal))
  
  oob_res[(((i-1)*5)+1):(i*5),] <- oob_cal
  
  # store results
  c_oob[i] <- glmnet::Cindex(oob_val, y = test_y)
  o_boot[i] <- glmnet::Cindex(train_pred, y = y) - glmnet::Cindex(orig_val, y = y_full)
}

###########################################################################
# Results
###########################################################################

# Stability
app_var <- rownames(lasso_cox_app$beta)[lasso_cox_app$beta@i +1]
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
