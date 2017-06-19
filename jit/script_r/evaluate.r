subFunc <- function(data.train=NULL, data.valid=NULL, data.eval=NULL,
                    method=NULL, meta=NULL, base=NULL, 
                    yname.numeric=NULL, yname.nominal=NULL, xnames=NULL, xname=NULL, model_control=NULL) {
  if (is.null(method)) {
    if (!is.null(meta)) { method <- paste(meta, base, sep="+") } else { method <- base }
  } else {
    base <- method
  }

  pred.train <- pred.valid <- NULL
  
  if (grepl("UP", method)) {
    tryCatch({
      data.train[, xname]

    }, error = function(err){
      print(err)
      print(xname)
      browser()
    }
    )
    pred.train <- 10000/(data.train[, xname]+1) # why 10000?
    pred.valid <- 10000/(data.valid[, xname]+1)
  } else if (grepl("EALR",method)) {
    formula <- as.formula(paste(yname.numeric, paste(xnames, collapse="+"), sep="~"))
    # browser()
    model <- glm(formula=formula, data=data.train, family=gaussian)
    try(model <- step(model, k=log(nrow(data.train)), trace=FALSE), silent=TRUE)
    
    if (length(attr(model$terms, "term.labels"))>=2) {
      err <- try(model.vif <- max(car::vif(model)), silent=TRUE)
      if (class(err)!="try-error") {
        if (!is.na(model.vif)) {
          if (model.vif > 10) { cat("VIF of EALR model is larger than 10.\n") }
        }
      }
    }
    
    pred.train <- predict(model, new=data.train)
    tryCatch({
      pred.valid <- predict(model, new=data.valid)
    }, error = function(err){
      print(err)
      browser()
    }
    )

    
    step.xnames <- attr(model$terms, "term.labels")
    for (aname in step.xnames) {
      data.train[, aname] <- (data.train[, aname]-min(data.train[, aname]))/(max(data.train[, aname])-min(data.train[, aname]))
      data.valid[, aname] <- (data.valid[, aname]-min(data.valid[, aname]))/(max(data.valid[, aname])-min(data.valid[, aname]))
    }
  } else if (!is.null(meta)) {
    data.train[, yname.nominal] <- as.factor(data.train[, yname.nominal])
    formula <- as.formula(paste(yname.nominal, paste(xnames, collapse="+"), sep="~"))
    ########################################
    ### Ensemble methods
    ########################################
    WC <- getWekaClassifier(meta)
    BC <- getWekaClassifier(base)
    model <- WC(formula=formula, data=data.train, control=RWeka::Weka_control(W=BC))
    
    if (is.null(yname.numeric)) {
      pred.train <- predict(model, new=data.train, type="probability")[, "1"]
      pred.valid <- predict(model, new=data.valid, type="probability")[, "1"]
    } else {
      pred.train <- predict(model, new=data.train)
      pred.valid <- predict(model, new=data.valid)
    }
  } else {
    data.train[, yname.nominal] <- as.factor(data.train[, yname.nominal])
    formula <- as.formula(paste(yname.nominal, paste(xnames, collapse="+"), sep="~"))
    
    WC <- getWekaClassifier(base)
    model <- NULL
    if (method=="IBk") {
      model <- WC(formula=formula, data=data.train, control=RWeka::Weka_control(K=8))
    } else {
      model <- WC(formula=formula, data=data.train, control=model_control)
    }
    
    if (is.null(yname.numeric)) {
      pred.train <- predict(model, new=data.train, type="probability")[, "1"]
      pred.valid <- predict(model, new=data.valid, type="probability")[, "1"]
      # pred.train <- predict(model, new=data.train)
      # pred.valid <- predict(model, new=data.valid)
      # levels(pred.valid) <- c(0,1)
    } else {
      pred.train <- predict(model, new=data.train)
      pred.valid <- predict(model, new=data.valid)
    }
  }
  
  pred.train <- as.vector(pred.train)
  pred.valid <- as.vector(pred.valid)
  train.dt <- valid.dt <- NULL
  if (grepl("UP", method)) { # why "UP" method use different datasets than others.
    train.dt <- data.frame(NUM=data.eval$fbug.old, REL=data.eval$fbug.old, LOC=data.eval$churn.fit.old, PRE=pred.train)
    valid.dt <- data.frame(NUM=data.eval$ebug.old, REL=data.eval$ebug.old, LOC=data.eval$churn.est.old, PRE=pred.valid)
  } else {
    train.dt <- data.frame(NUM=data.eval$fbug, REL=data.eval$fbug, LOC=data.eval$churn.fit, PRE=pred.train)
    valid.dt <- data.frame(NUM=data.eval$ebug, REL=data.eval$ebug, LOC=data.eval$churn.est, PRE=pred.valid)
  }
  
  sorted           <- FALSE
  worstcase        <- TRUE  ### compute the worst performance for unsupervised models
  bestcase         <- FALSE
  LOCUP            <- FALSE
  allpercentcutoff <- TRUE
  sub.Popt <- sub.ACC <-sub.Prec.F1<- NULL
  if (grepl("UP", method)) {
    sub.Popt <- ComputePopt(sorted=sorted, data=valid.dt, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sub.ACC <- ComputeACC(sorted=sorted, data=valid.dt, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sub.Prec.F1 <- ComputePrecF1(sorted=sorted, data=valid.dt, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
  } else {

    sub.Popt <- ComputePopt(sorted=sorted, data=valid.dt)
    sub.ACC <- ComputeACC(sorted=sorted, data=valid.dt)
    sub.Prec.F1 <- ComputePrecF1(sorted=sorted, data=valid.dt)
  }
  
  sub.Prec <- sub.Prec.F1$Prec
  sub.F1 <- sub.Prec.F1$F1
  if(grepl("_Sup",method)){
    method <- "OneR.UP"
  }
  names(sub.Popt) <- paste(method, "Popt", sep=".")
  names(sub.ACC) <- paste(method, "ACC", sep=".")
  names(sub.Prec) <-paste(method,"Prec", sep=".")
  names(sub.F1) <-paste(method,"F1", sep=".")
  return (list(Popt= sub.Popt, ACC=sub.ACC, Prec=sub.Prec, F1=sub.F1))
  
}