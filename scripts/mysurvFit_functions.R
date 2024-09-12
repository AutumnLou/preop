aggregate.survFit <- function(survFit){
  agg_surv <- rowMeans(survFit$surv)
  times <- as.numeric(survFit$time)
  mylist <- list('time' = times, 'surv' = agg_surv)
  return(mylist)
}

# Example
#library(survival)
#library(mboost)
#zx <- na.omit(lung)
#zx$status <- as.integer(zx$status == 2)
#y <- Surv(zx$time, zx$status)
#x <- zx
#x$time <- NULL
#x$status <- NULL
#x <- as.matrix(x)
#boost_cox <- glmboost(x=x, y=y, family = CoxPH())

#km <- survfit(y ~ 1)
#survfit <- survFit(boost_cox, x)
#aggregate.survFit(survfit)
#index_60 <- which(abs(aggregate.survFit(survfit)$time-60)==min(abs(aggregate.survFit(survfit)$time-60)))
#surv_60 <- aggregate.survFit(survfit)$surv[index_60]
