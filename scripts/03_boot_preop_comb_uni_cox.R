###########################################################################
# PhD - Pre Op Journal Article
## Traditional Cox Bootstrapped
## mRNA data univariate feature selection
# 03_boot_preop_comb_uni_cox.R
###########################################################################

# Load the dataset
source("scripts/02_data_comb_preop.R")

# Load Libraries/Functions
source("scripts/mycalibration_function2.R")
require(RegParallel)

###########################################################################
# Prep and calculation of Apparent C-Index
###########################################################################
#Feature Selection
options(scipen=10)
options(digits=6)
data <- comb_data[,8:ncol(comb_data)]
colnames(data) <- paste0('gene', 1:ncol(data))
rownames(data) <- paste0('sample', 1:nrow(data))
variables = colnames(data)[1:ncol(data)]
data$time <- comb_data[,6]
data$alive <- comb_data[,7]

res4 <- RegParallel(
  data = data,
  formula = 'Surv(time, as.integer(alive)) ~  [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = variables,
  blocksize = 2000,
  p.adjust = "BH"
)

res5 <- res4[!is.na(res4$P),]
res5 <- res5[order(res5$LogRank.adjust, decreasing = FALSE),]
final <- head(res5, 50)
final <- subset(final, LogRank.adjust < 0.05)
probes <- gsub('[^0-9.-]', '', final$Variable)

app_genes_columns <- sort(as.numeric(probes) + 7)

gene_data <- comb_data[, c(1:7, app_genes_columns)]

# Model w/ all 133 obs (app)
app_start_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ 1, data = gene_data, x = TRUE)
app_full_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ ., data = gene_data, x = TRUE)
app_fit_step <- MASS::stepAIC(app_start_cox, direction = "both",
                              scope = app_full_cox$formula,
                              steps = 10)

# Discrimination
c_app <- concordance(app_fit_step)$concordance


###########################################################################
# Resampling Framework
###########################################################################
# 100 Bootstrap resamples
set.seed(123)
B <- 100
bootIndex <- createResample(factor(comb_data$BCR_Event), times = B, list = TRUE)

c_oob <- numeric(B)
o_boot <- numeric(B)

oob_har_res <- data.frame(group = numeric(B*5), mean.p = numeric(B*5),
                          km = numeric(B*5), lower = numeric(B*5),
                          upper = numeric(B*5))

oob_res <- data.frame(pred = numeric(B*5), obs = numeric(B*5),
                      mean.obs = numeric(B*5), low = numeric(B*5),
                      up = numeric(B*5), means.low = numeric(B*5),
                      means.up = numeric(B*5), group = numeric(B*5))

stab.v <- rep.int(0,length(colnames(comb_data[,c(1:5,8:ncol(comb_data))])))
names(stab.v) <- colnames(comb_data[,c(1:5,8:ncol(comb_data))])
stab.uni.v <- rep.int(0,length(colnames(comb_data[,c(1:5,8:ncol(comb_data))])))
names(stab.uni.v) <- colnames(comb_data[,c(1:5,8:ncol(comb_data))])
haz <- data.frame(matrix(ncol = 135, nrow = 0))
colnames(haz) <- rownames(comb_data)

for(i in 1:length(bootIndex)) {
  print(i)
  # define Training / Validation data
  bag <- comb_data[bootIndex[[i]],]
  oob <- comb_data[-bootIndex[[i]],]
  
  #Feature Selection
  options(scipen=10)
  options(digits=6)
  data <- bag[,8:ncol(bag)]
  colnames(data) <- paste0('gene', 1:ncol(data))
  rownames(data) <- paste0('sample', 1:nrow(data))
  variables = colnames(data)[1:ncol(data)]
  data$time <- bag[,6]
  data$alive <- bag[,7]
  
  res4 <- RegParallel(
    data = data,
    formula = 'Surv(time, as.integer(alive)) ~  [*]',
    FUN = function(formula, data)
      coxph(formula = formula,
            data = data,
            ties = 'breslow',
            singular.ok = TRUE),
    FUNtype = 'coxph',
    variables = variables,
    blocksize = 2000,
    p.adjust = "BH"
  )
  
  
  res5 <- res4[!is.na(res4$P),]
  res5 <- res5[order(res5$LogRank.adjust, decreasing = FALSE),]
  final <- head(res5, 50)
  final <- subset(final, LogRank.adjust < 0.05)
  probes <- gsub('[^0-9.-]', '', final$Variable)
  
  genes_columns <- sort(as.numeric(probes) + 7)
  
  bag_gene_data <- bag[, c(1:7, genes_columns)]
  
  # Univariate Stability
  genes <- colnames(comb_data[,c(genes_columns)])
  for(k in 1:length(genes)) {
    stab.uni.v[which(names(stab.uni.v) == genes[k])] <- stab.uni.v[which(names(stab.uni.v) == genes[k])] + 1
  }
  
  # Stepwise Selection
  start_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ 1, data = bag_gene_data, x = TRUE)
  
  full_cox <- coxph(Surv(BCR_FreeTime, BCR_Event) ~ ., data = bag_gene_data, x = TRUE)
  
  fit_step <- MASS::stepAIC(start_cox, direction = "both", scope = full_cox$formula,
                            trace = 0, steps = 10)
  
  # Stability
  terms <- attr(terms(fit_step), "term.labels")
  for(j in 1:length(terms)) {
    stab.v[which(names(stab.v) == terms[j])] <- stab.v[which(names(stab.v) == terms[j])] + 1
  }
  
  # Validation
  oob_val <- concordance(fit_step, newdata = oob)
  orig_val <- concordance(fit_step, newdata = comb_data)
  
  # For ROC and DCA
  preds <- predict(fit_step, newdata = oob)
  preds <- (preds - min(preds))/(max(preds)-min(preds))
  for(k in 1:length(preds)) {
    haz[i,which(colnames(haz) == names(preds)[k])] <- preds[k]
  }
  
  ## OOB
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

# Univariate Stability
uni_var <- colnames(comb_data[,c(app_genes_columns)])
stab.uni.v <- sort(stab.uni.v, decreasing = TRUE)
barchart(head(stab.uni.v, 10),
         xlab = "Number of selections accross resamples",
         xlim = 0:100)

# Stability
app_var <- attr(terms(app_fit_step), "term.labels")
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
roc_dat <- data.frame(surv = comb_data$BCR_FreeTime,
                      time = comb_data$BCR_Event,
                      marker = colMeans(haz, na.rm = TRUE))
