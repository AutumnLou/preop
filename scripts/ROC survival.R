roc_dat <- read.csv("data/roc.csv", header = TRUE, row.names = 1)
colnames(roc_dat)[2] = "status"

####### 3 years
roc_dat$status3yr <- 0
roc_dat$status3yr[is.na(roc_dat$status)] <- NA
roc_dat$status3yr[roc_dat$surv<=36 & roc_dat$status==0] <- NA 
roc_dat$status3yr[roc_dat$surv<=36 & roc_dat$status==1] <- 1 

############# 5 years
roc_dat$status5yr <- 0
roc_dat$status5yr[is.na(roc_dat$status)] <- NA
roc_dat$status5yr[roc_dat$surv<=60 & roc_dat$status==0] <- NA 
roc_dat$status5yr[roc_dat$surv<=60 & roc_dat$status==1] <- 1 

#### ROC curve
library(pROC)

# 3 year
auc.3.cox <- pROC::roc(roc_dat$status3yr ~ roc_dat$Cox, ci = TRUE)
auc.3.lasso <- pROC::roc(roc_dat$status3yr ~ roc_dat$LASSO, ci = TRUE)
auc.3.boost <- pROC::roc(roc_dat$status3yr ~ roc_dat$Boosted, ci = TRUE)
auc.3.rsf <- pROC::roc(roc_dat$status3yr ~ roc_dat$RSF, ci = TRUE)
auc.3.cox$auc; ci(auc.3.cox)
auc.3.lasso$auc; ci(auc.3.lasso)
auc.3.boost$auc;ci(auc.3.boost)
auc.3.rsf$auc; ci(auc.3.rsf)

# 5 year
auc.5.cox <- pROC::roc(roc_dat$status5yr ~ roc_dat$Cox, ci = TRUE)
auc.5.lasso <- pROC::roc(roc_dat$status5yr ~ roc_dat$LASSO, ci = TRUE)
auc.5.boost <- pROC::roc(roc_dat$status5yr ~ roc_dat$Boosted, ci = TRUE)
auc.5.rsf <- pROC::roc(roc_dat$status5yr ~ roc_dat$RSF, ci = TRUE)
auc.5.cox$auc; ci(auc.5.cox, )
auc.5.lasso$auc; ci(auc.5.lasso)
auc.5.boost$auc;ci(auc.5.boost)
auc.5.rsf$auc; ci(auc.5.rsf)

## Plot
par(mar=c(8,5,2,2),mfrow = c(1,2), cex.lab = 1.2, font.lab = 2, cex.axis = 1.2, font.axis = 2)
plot(1-auc.3.cox$specificities, auc.3.cox$sensitivities, type="l", col = "gold", lwd = 2, xlim=c(0,1), ylim=c(0,1),
     xlab="",
     ylab="Sensitivity",
     main="ROC: Year = 3")
abline(0,1)
mtext(paste("1 - Specificity", "\n", "AUC: Cox = ",round(auc.3.cox$auc,3),
            " LASSO = ",round(auc.3.lasso$auc,3),"\n",
            " Boosted = ",round(auc.3.boost$auc,3),
            " RSF = ",round(auc.3.rsf$auc,3)), side = 1, line = 5, font = 2, cex = 1.2)
points(1-auc.3.lasso$specificities, auc.3.lasso$sensitivities, type="l", col="red", lwd = 2)
points(1-auc.3.boost$specificities, auc.3.boost$sensitivities, type="l", col="blue", lwd = 2)
points(1-auc.3.rsf$specificities, auc.3.rsf$sensitivities, type="l", col="black", lwd = 2)
legend("bottomright", legend=c("Cox", "LASSO", "Boosted", "RSF"),
       col=c("gold", "red", "blue", "black"), lty=1, lwd = 2, cex=2, bty = "n")

plot(1-auc.5.cox$specificities, auc.5.cox$sensitivities, type="l", col = "gold", lwd = 2, xlim=c(0,1), ylim=c(0,1),
     xlab="",
     ylab="Sensitivity",main="ROC: Year = 5")
abline(0,1)
mtext(paste("1 - Specificity", "\n", "AUC: Cox = ",round(auc.5.cox$auc,3),
            " LASSO = ",round(auc.5.lasso$auc,3),"\n",
            " Boosted = ",round(auc.5.boost$auc,3),
            " RSF = ",round(auc.5.rsf$auc,3)), side = 1, line = 5, font = 2, cex = 1.2)
points(1-auc.5.lasso$specificities, auc.5.lasso$sensitivities, type="l", col="red", lwd = 2)
points(1-auc.5.boost$specificities, auc.5.boost$sensitivities, type="l", col="blue", lwd = 2)
points(1-auc.5.rsf$specificities, auc.5.rsf$sensitivities, type="l", col="black", lwd = 2)
legend("bottomright", legend=c("Cox", "LASSO", "Boosted", "RSF"),
       col=c("gold", "red", "blue", "black"), lty=1, lwd = 2, cex=2, bty = "n")

# Pairwise comparisons
roc.test(auc.3.cox, auc.3.lasso, paired = TRUE)
roc.test(auc.3.cox, auc.3.boost, paired = TRUE)
roc.test(auc.3.cox, auc.3.rsf, paired = TRUE)
roc.test(auc.3.lasso, auc.3.boost, paired = TRUE)
roc.test(auc.3.lasso, auc.3.rsf, paired = TRUE)
roc.test(auc.3.boost, auc.3.rsf, paired = TRUE)

roc.test(auc.5.cox, auc.5.lasso, paired = TRUE)
roc.test(auc.5.cox, auc.5.boost, paired = TRUE)
roc.test(auc.5.cox, auc.5.rsf, paired = TRUE)
roc.test(auc.5.lasso, auc.5.boost, paired = TRUE)
roc.test(auc.5.lasso, auc.5.rsf, paired = TRUE)
roc.test(auc.5.boost, auc.5.rsf, paired = TRUE)
