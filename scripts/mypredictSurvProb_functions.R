aggregate.predictSurvProb <- function(predictSurvProb, times){
  agg_surv = colMeans(predictSurvProb)
  mylist <- list('time' = times, 'surv' = agg_surv)
  return(mylist)
}

## Example
#library(survival)
#library(randomForestSRC)

#zx = na.omit(lung)
# ... make sure 1 codes for 'death'
#zx$status = as.integer(zx$status==2)

#rsf = rfsrc(Surv(time,status)~., data=zx)
#km = survfit(Surv(zx$time,zx$status)~1)
#overall_surv <- pec::predictSurvProb(rsf, newdata = zx, times = km$time)
#sc = aggregate.predictSurvProb(overall_surv, times = km$time)