# Calibration function for overall survival
my.calibration2 <- function(sc,km,
                           Surv.output,
                           q=seq(0,1,by=.2)){
  # Input:
  # 	Surv.output is the output of Surv(), survival object, times with censorship information
  # 	sc = survfit(cox), creates survival curves for the fitted data from the fitted model
  # 	km is a Kaplan-Meier fit obtained from survfit(). (estimate the observed survival function)
  #   Group into quintiles by default (q).
  #   observed, lower, and upper are averages from KM estimates, 
  #   this could be reconsidered.
  # 
  # groups:
  # ... per prediction quantiles:
  # qs = quantile(predicted,q)
  # gp = cut(predicted,qs,labels=F)
  # ... per time:
  events = as.numeric(Surv.output)[1:(length(as.numeric(Surv.output))/2)] # observed times of each observation
  times = as.numeric(quantile(events,q)) # quintiles of the the observed times
  gp = cut(events,times,labels=F) # quintile group each individual belongs to
  # observed data:
  observed = km$surv # estimate of the observed survival probabilities - is this needed?
  refs = approx(x=km$time,y=km$surv,xout=times)$y # estimate of the observed probability of survival at each quintile - is this needed?
  lowers = approx(x=km$time,y=km$lower,xout=times)$y # estimate of the observed probability of survival (lower bound of 95% CI) at each quintile - is this needed?
  uppers = approx(x=km$time,y=km$upper,xout=times)$y # estimate of the observed probability of survival (upper bound of 95% CI) at each quintile - is this needed?
  L = length(refs) # for quintiles L = 6, 0%, 20%, 40%, 60%, 80%, 100% - is this needed?
  refs = refs[-1] # remove the first value related to 0% - is this needed?
  lowers = lowers[-1] # remove the first value related to 0% - is this needed?
  uppers = uppers[-1] # remove the first value related to 0% - is this needed?
  # predicted data:
  sh = approx(x=sc$time,y=sc$surv,xout=events)$y # estimate of the predicted probability of survival (mean) for each observation
  kh = approx(x=km$time,y=km$surv,xout=events)$y # estimate of the observed probability of survival (mean) for each observation
  kh.cihi = approx(x=km$time,y=km$upper,xout=events)$y # estimate of the observed probability of survival (upper bound of 95% CI) for each observation
  kh.cilo = approx(x=km$time,y=km$lower,xout=events)$y # estimate of the observed probability of survival (lower bound of 95% CI) for each observation
  G = max(gp,na.rm=T) # number of groups (quintile = 5) 
  means.p = means.k = numeric(G)
  cihi.k = cilo.k = numeric(G)
  bxp = matrix(NA,ncol=G,nrow=length(sh))
  for(g in 1:G){
    ig = which(gp==g) # setting group = g
    means.p[g] = mean(sh[ig],na.rm=TRUE)  # mean of estimate of the predicted probability of survival (mean) for group g
    means.k[g] = mean(kh[ig],na.rm=TRUE) # mean of estimate of the observed probability of survival (mean) for group g
    cihi.k[g] = mean(kh.cihi[ig],na.rm=TRUE) # mean of estimate of the observed probability of survival (upper bound of 95% CI) for group g
    cilo.k[g] = mean(kh.cilo[ig],na.rm=TRUE) # mean of estimate of the observed probability of survival (lower bound of 95% CI) for group g
    bxp[ig,g] = sh[ig] # estimate of the predicted probability of survival (mean) for each observation in group g - is this needed?
  }
  mylist <- list('pred' = means.p,
                 'obs' = refs, # is this needed?
                 'mean.obs' = means.k,
                 'low' = lowers, # is this needed?
                 'up' = uppers, # is this needed?
                 'means.low' = cilo.k,
                 'means.up' = cihi.k)
  return(mylist)
}

# Function for plotting calibration
my.calibration.plot <- function(mean.p, km,
                                lower, upper,
                                xlab = "Predict",
                                ylab = "Observed",
                                font=1,
                                cex=1){
  plot(c(0:1),c(0:1),t='n',
       xlab=xlab,ylab=ylab,
       font.lab=font,font.axis=font,cex.lab=cex)
  abline(a=0,b=1,lwd=.5)
  points(mean.p,km,pch=20,t='b',cex=cex)
  segments(x0=mean.p,y0=lower,y1=upper)
}

# Function for plotting calibration
my.calibration.plot2 <- function(cal,xlab = "Predict", ylab = "Observed"){
  mean.p=cal$pred
  km=cal$mean.obs
  lower=cal$low
  upper=cal$up 
  font=1
  cex=1
  plot(c(0:1),c(0:1),t='n',
       xlab=xlab,ylab=ylab,
       font.lab=font,font.axis=font,cex.lab=cex)
  abline(a=0,b=1,lwd=.5)
  points(mean.p,km,pch=20,t='b',cex=cex)
  segments(x0=mean.p,y0=lower,y1=upper)
}

####################################################
# # Examples
# 
# # Load library and data set
# library(survival)
# library(dplyr)
# zx <- na.omit(lung)
# zx$status = as.integer(zx$status==2)
# # Fit a model
# cox = coxph(Surv(time,status)~., data=zx)
# # Fitted model survival curve
# surv <- survfit(cox, newdata = zx)
# 
# ## Single time-point example
# 
# # Summary at time = 100 (100 days)
# surv100 <- summary(surv, times = 100)
# # Dataframe of calibration data
# cal.df <- data.frame(time = zx$time,
#                       status = zx$status,
#                       surv100 = t(surv100$surv),
#                       row.names = rownames(zx))
# # Order by predicted survival
# cal.df <- cal.df[order(cal.df$surv100),]
# # Group into quitiles
# cal.df$group <- cut(1:nrow(cal.df), 5, labels = FALSE)
# cal.df$groupq <- cut(cal.df$surv100, as.numeric(quantile(cal.df$surv100, probs = seq(0, 1, 0.2))), labels = FALSE)
# # For each quintile find the mean predicted survival,
# # KM estimate for the quintile and it's upper and lower limit
# cal.df <- cal.df %>% 
#   group_by(group) %>% 
#   mutate(mean.p = mean(surv100),
#          km = summary(survfit(Surv(time,status)~1), times = 100)$surv,
#          lower = summary(survfit(Surv(time,status)~1), times = 100)$lower,
#          upper = summary(survfit(Surv(time,status)~1), times = 100)$upper) %>%
#   summarise(mean.p = mean(mean.p),
#             km = mean(km),
#             lower = mean(lower),
#             upper = mean(upper))
# 
# # Plot calibration
# my.calibration.plot(cal.df$mean.p,
#                     cal.df$km,
#                     cal.df$lower,
#                     cal.df$upper,
#                     xlab="Predicted 100 day Surv",
#                     ylab="Observed 100 day Surv")
# ## Overall survival
# km <- survfit(Surv(zx$time,zx$status)~1)
# Surv.output = Surv(zx$time,zx$status)
# sc = aggregate(survfit(cox, newdata = zx))
# cal <- my.calibration2(sc, km, Surv.output)
# # Plot calibration
# ### WHAT WAS USED FOR PAPER:
# my.calibration.plot(cal$pred,
#                     cal$obs,
#                     cal$low,
#                     cal$up,
#                     xlab="Predicted",
#                     ylab="Observed")
# ### POSSIBLE UPDATED CURVE:
# my.calibration.plot2(cal,
#                     xlab="Predicted",
#                     ylab="Observed")

