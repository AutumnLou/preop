###########################################################################
# RFE_for_RSF_function.R
###########################################################################
# Function
rsf_rfe <- function(data){
  
  require(randomForestSRC)
  # Full model
  out.rf = rfsrc(formula = Surv(time,status) ~ .,
                 data = data, importance = TRUE)
  # C-Index container
  rfe.cindex.rf <- 1 - mean(out.rf$err.rate, na.rm = TRUE)
  
  # Input set (excluding Y variables)
  input_set = data
  input_set$time = NULL
  input_set$status = NULL
  
  # Remaining X variables
  rfe.set = input_set
  # minP = minimum number of features to keep in model:
  minP = 3 # (you can go to, say, 5 to speed things up)
  out.rf.rfe = out.rf # the full model output
  #
  # RFE:
  #
  rfe.auc.rf = NULL
  ejects = NULL
  while(ncol(rfe.set)>=(minP+1)){
    #cat("\t> rf-rfe: ",ncol(rfe.set),"features left...\n")
    # Find the least important variable
    imp = out.rf.rfe$importance
    vnms = names(imp)
    imin = which.min(imp)
    least.imp = vnms[imin]
    # Make new data set with only the remaining important variables
    ejects = c(ejects, least.imp)
    rfe.nms = setdiff(names(input_set), ejects)
    rfe.set = input_set[, rfe.nms]
    use.set = rfe.set
    use.set$status = data$status
    use.set$time = data$time
    
    # Internal RSF
    out.rf.rfe = rfsrc(formula = Surv(time,status) ~ .,
                       data = use.set, importance = TRUE)
    # Internal c-index
    c.index <- 1 - mean(out.rf.rfe$err.rate, na.rm = TRUE)
    rfe.cindex.rf = c(rfe.cindex.rf, c.index)
  }
  # pick the smallest model out of those with highest c-index:
  opti = rev(which(rfe.cindex.rf==max(rfe.cindex.rf)))[1]
  # reproduce backward elimination up to optimum:
  fset = setdiff(names(input_set),ejects[1:opti]) 
  
  # ... and then re-evaluate your final model
  zx <- data[,fset]
  zx$time <-data$time
  zx$status <-data$status
  rsf = rfsrc(formula = Surv(time,status) ~ .,
             data = zx, importance = TRUE)
  return(rsf)
}

# Example
#data = na.omit(lung)
#data$status = as.integer(data$status==2)
#rsf<-rsf_rfe(data)


# Example
#data = na.omit(lung)
#data$status = as.integer(data$status==2)
#Y = Surv(data$time,data$status)
#input_set = data[,-c(2:3)]
#out.rf = rfsrc(formula = Surv(time,status) ~ .,
#               data = data, importance = TRUE)
#rfe.cindex.rf <- 1 - mean(out.rf$err.rate, na.rm = TRUE)

#input_set = your_input_dataset_excluding_Y # we need this copy
#
# initialisations:
#
#rfe.set = input_set
# minP = minimum number of features to keep in model:
#minP = 3 # (you can go to, say, 5 to speed things up)
#out.rf.rfe = out.rf # the full model output
#
# RFE:
#
#rfe.auc.rf = NULL
#ejects = NULL
#while(ncol(rfe.set)>=(minP+1)){
#  cat("\t> rf-rfe: ",ncol(rfe.set),"features left...\n")
  # you have to adapt these lines to your model output format:
#  imp = out.rf.rfe$importance
#  vnms = names(imp)
#  imin = which.min(imp)
#  least.imp = vnms[imin]
  # 
#  ejects = c(ejects, least.imp)
#  rfe.nms = setdiff(names(input_set), ejects)
#  rfe.set = input_set[, rfe.nms]
#  use.set = rfe.set
#  use.set$status = data$status
#  use.set$time = data$time
  # rfe.vvb = test_set[, rfe.nms] # if required
  # do not re-tune RF at this point	
  # you have to adapt these lines to your model output format:
#  out.rf.rfe = rfsrc(formula = Surv(time,status) ~ .,
#                     data = use.set, importance = TRUE)
  # here you must below lines change to use the corresponding deviance:
  # - I compute accuracy:
#  c.index <- 1 - mean(out.rf.rfe$err.rate, na.rm = TRUE)
#  rfe.cindex.rf = c(rfe.cindex.rf, c.index)
#}
# pick the smallest model out of those with highest Accuracy:
#opti = rev(which(rfe.cindex.rf==max(rfe.cindex.rf)))[1]
# reproduce backward elimination up to optimum:
#fset = setdiff(names(input_set),ejects[1:opti]) 

# ... and then re-evaluate your final model
#zx <- data[,fset]
#zx$time <-data$time
#zx$status <-data$status
#rf = rfsrc(formula = Surv(time,status) ~ .,
#               data = data, importance = TRUE)
#cindex.rf <- 1 - mean(rf$err.rate, na.rm = TRUE)

