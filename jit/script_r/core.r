# setwd("~/Github/Research/TuneFSE2016/jit/script_r")
# RWeka::WPM("install-package", "./packages/linearForwardSelection1.0.1.zip")
# RWeka::WPM("install-package", "./packages/ridor1.0.1.zip")
# RWeka::WPM("install-package", "./packages/EMImputation1.0.1.zip")
# RWeka::WPM("install-package", "./packages/citationKNN1.0.1.zip")
# RWeka::WPM("install-package", "./packages/isotonicRegression1.0.1.zip")
# RWeka::WPM("install-package", "./packages/paceRegression1.0.1.zip")
# RWeka::WPM("install-package", "./packages/leastMedSquared1.0.1.zip")
# RWeka::WPM("install-package", "./packages/RBFNetwork1.0.8.zip")
# RWeka::WPM("install-package", "./packages/conjunctiveRule1.0.4.zip")
# RWeka::WPM("install-package", "./packages/rotationForest1.0.2.zip")

source("evaluate.r")
getWekaClassifier <- function(name) {
  WC <- NULL
  
  if (name %in% c("LinearRegression", "Logistic", "SMO", "IBk", "LBR", "AdaBoostM1", "Bagging", "LogitBoost", "MultiBoostAB", "Stacking", "CostSensitiveClassifier", "JRip", "M5Rules", "OneR", "PART", "J48", "LMT", "M5P", "DecisionStump", "SimpleKMeans")) {
    WC <- get(name, asNamespace("RWeka"))
  } else {
    if (name=="NaiveBayes") {
      WC <- RWeka::make_Weka_classifier("weka/classifiers/bayes/NaiveBayes")
    } else if (name %in% c("GaussianProcesses", "SimpleLogistic", "SimpleLinearRegression", "IsotonicRegression", "LeastMedSq", "MultilayerPerceptron", "PaceRegression", "PLSClassifier", "RBFNetwork", "SMOreg")) {
      WC <- RWeka::make_Weka_classifier(paste("weka/classifiers/functions/", name, sep=""))
    } else if (name %in% c("KStar", "LWL")) {
      WC <- RWeka::make_Weka_classifier(paste("weka/classifiers/lazy/", name, sep=""))
    } else if (name %in% c("Bagging", "AdaBoostM1", "RotationForest", "RandomSubSpace", "AdditiveRegression", "AttributeSelectedClassifier", "CVParameterSelection", "Grading", "GridSearch", "MultiScheme", "RegressionByDiscretization", "StackingC", "Vote")) {
      WC <- RWeka::make_Weka_classifier(paste("weka/classifiers/meta/", name, sep=""))
    } else if (name %in% c("ConjunctiveRule", "DecisionTable", "ZeroR", "Ridor")) {
      WC <- RWeka::make_Weka_classifier(paste("weka/classifiers/rules/", name, sep=""))
    } else if (name %in% c("REPTree", "UserClassifier", "RandomForest")) {
      WC <- RWeka::make_Weka_classifier(paste("weka/classifiers/trees/", name, sep=""))
    }
  }
  return(WC)
}

getTuneSapce <- function(){
  
  ## J48 tuning space
  C <- c(min=0.0001, max=0.5)
  M <- c(min=2L, max=20L)
  J48 <- data.frame(C,M)
  names(J48) <- c("C","M")
  
  ## IBk tuning space
  K <- c(min=1L, max=20L)
  W <- c(min=0L, max=20L)
  IBk <- data.frame(K,W)
  names(IBk) <- c("K","W")
  
  ## RandomForests tuning space
  K <- c(min=0L,max=10L)
  M <- c(min=1L, max=20L)
  V <- c(min=0.0001, max=0.01)
  depth <- c(min=0L, max=10L)
  RandomForest <- data.frame(K,M,V,depth)
  names(RandomForest) <- c("K","M","V","depth")
  
  
  temp <- list(J48=J48, IBk=IBk, RandomForest=RandomForest)
}
coreExperiment <- function(fit, est, sampling=TRUE, seed=0) {
  ###### evaluation criterias ######
  Popt <- ACC <- F1 <- Prec<- NULL
  
  # Remove the commit date
  fit$commitdate <- est$commitdate <- NULL
  
  fit.old <- fit
  est.old <- est
  ### sampling
  if (sampling) {
    fit <- doSampling(fit, "bug", seed=seed)
  }
  # cat("training data size:",nrow(fit),"testing data size:", nrow(est),"\n")

  
  ### calc churn
  churn.fit.old <- decChurn(fit.old); churn.fit <- decChurn(fit)
  churn.est.old <- decChurn(est.old); churn.est <- decChurn(est)
  
  fbug.old <- fit.old$bug; fbug <- fit$bug
  ebug.old <- est.old$bug; ebug <- est$bug

  # log transformation for each measure except fix
  idx <- charmatch(c("ns","nm","la","ld","lt","exp","rexp","sexp","ndev","pd","npt","entropy"), colnames(fit))
  fit[,idx] <- fit[,idx] + 1
  idx <- charmatch(c("fix","bug"), colnames(fit))
  fit[,-c(idx)] <- apply(fit[,-c(idx)], 2, log)

  
  ###### remove the variable that has high correlation to other variables
  # idx <- charmatch(c("nm","rexp"), colnames(fit))
  # fit <- fit[,-c(idx)]
  
  ###### remove la and ld from the prediction model
  idx <- charmatch(c("la","ld"), colnames(fit))
  fit <- fit[,-c(idx)]
  
  # ### remove those variables with 0 variance
  va <- apply(fit, 2, var)
  fit <- fit[,((!is.na(va))&(va!=0))] ### fit <- fit[,!(va==0)]
  
  # 
  err <- try(vifstep.names <- as.character(usdm::vifstep(fit[, setdiff(colnames(fit), "bug")], th=10)@results$Variables), silent=TRUE)
  if (class(err)!="try-error") {
    fit <- fit[, c(vifstep.names, "bug")]
  }

  ### VIF analysis
  idx <- charmatch(c("bug"), colnames(fit))
  x <- fit[,-c(idx)]
  err <- try(x <- remVarsByVIF(x), silent=TRUE)
  if (class(err)!="try-error") {
    fit <- cbind(x, fit["bug"])
  }
  
  ### log transformation for est
  idx <- charmatch(c("ns","nm","la","ld","lt","exp","rexp","sexp","ndev","pd","npt","entropy"), colnames(est))
  est[,idx] <- est[,idx] + 1
  idx <- charmatch(c("fix","bug"), colnames(est))
  est[,-c(idx)] <- apply(est[,-c(idx)], 2, log)
  
  
  ###### remove la and ld from the prediction model
  idx <- charmatch(c("la","ld"), colnames(est))
  est <- est[,-c(idx)]
  
  # va <- apply(est, 2, var)
  # est <- est[,((!is.na(va))&(va!=0))] ### est <- est[,!(va==0)]

  ### evaluate the prediction performance when considering effort (effort=churn)
  fit$bugdensity <- fit$bug/(churn.fit+1)
  est$bugdensity <- est$bug/(churn.est+1)
  
  fit.old$bugdensity <- fit.old$bug/(churn.fit.old+1)
  est.old$bugdensity <- est.old$bug/(churn.est.old+1)
  
  # browser()
  # cat("neg/pos: ",(length(est.old$bug)-sum(est.old$bug))/sum(est.old$bug),"\n")
  
  ### obtain original values for each variable
  est.old$lt  <- est.old$lt  * est.old$nf
  est.old$la  <- est.old$la  * est.old$lt
  est.old$ld  <- est.old$ld  * est.old$lt
  est.old$nuc <- est.old$npt * est.old$nf
  
  fit.old$lt  <- fit.old$lt  * fit.old$nf
  fit.old$la  <- fit.old$la  * fit.old$lt
  fit.old$ld  <- fit.old$ld  * fit.old$lt
  fit.old$nuc <- fit.old$npt * fit.old$nf
  
  
  ### remove la, ld from fit.old and est.old

  fit.old <-fit.old[,-c(idx)]
  est.old <-est.old[,-c(idx)]
  
  
  
  ### make a list of data for evaluation FOO..
  data.eval <- list(fbug=fbug, fbug.old=fbug.old, ebug=ebug, ebug.old=ebug.old,
                    churn.fit.old= churn.fit.old, churn.fit=churn.fit, churn.est.old=churn.est.old,
                    churn.est=churn.est)
  
  
  yname.numeric <- "bugdensity"; yname.nominal <- "bug"
  xnames <- setdiff(colnames(fit), c("bug", "bugdensity"))
  
  
  T.idx <- sample(which(fit$bug==1),round(0.75*length(which(fit$bug==1))))
  F.idx <- sample(which(fit$bug==0),round(0.75*length(which(fit$bug==0))))
  shuffle.idx <- c(T.idx, F.idx)
  train.data <- fit[shuffle.idx,]
  # tune.data <- fit[-shuffle.idx,]
  
  
  ############# EALR, using the original features.....
  result <-subFunc(method="EALR", data.train=fit, data.valid=est, data.eval=data.eval,yname.numeric=yname.numeric, xnames=xnames)
  Popt <-c(Popt, result$Popt)
  ACC <-c(ACC, result$ACC)
  Prec <- c(Prec, result$Prec)
  F1 <- c(F1, result$F1)
  

  
  supervised.methods <- c("RandomForest","J48","IBk")
  tuned.method <-NULL
  for (method in supervised.methods) { ###
    result <-subFunc(data.train=fit, data.valid=est, data.eval=data.eval, method=method, yname.numeric=NULL, yname.nominal=yname.nominal, xnames=xnames)
    Popt <-c(Popt, result$Popt)
    ACC <-c(ACC, result$ACC)
    Prec <- c(Prec, result$Prec)
    F1 <- c(F1, result$F1)
  }

  cnames <- c("ns", "nm", "nf", "entropy", "lt", "fix", "ndev", "pd",  "nuc", "exp", "rexp", "sexp","PC1","PC2","PC3")
  labels <- c("NS", "ND", "NF", "Entropy", "LT", "FIX", "NDEV", "AGE", "NUC", "EXP", "REXP", "SEXP","PC1","PC2","PC3")
  cnames <- setdiff(colnames(fit.old), c("bug", "bugdensity"))
  labels <- setdiff(colnames(fit.old), c("bug", "bugdensity"))
  for (i in seq(cnames)) {
    xname <- cnames[i]; label <- labels[i]
    methods <- paste(label, "UP", sep=".")

    for (method in methods) {
      result <-subFunc(method=method, data.train=fit.old, data.valid=est.old, data.eval=data.eval,xname=xname)
      Popt <-c(Popt, result$Popt)
      ACC <-c(ACC, result$ACC)
      Prec <- c(Prec, result$Prec)
      F1 <- c(F1, result$F1)
    }
  }
  
  #################### Feature Prunining   ####################
  findFeature <-function(theDF, goal=NULL){
    getMean <- function(theDF){
      x <- apply(theDF, c(1), mean)
      method <- names(x[x==max(x)])
      return(method)
    }
    method <- NULL
    if(is.null(goal)){
      ### calculate mean of all the performance,return the highest.
      # x <- apply(theDF, c(1), mean)
      method <- getMean(theDF)
    }else{
      x <- apply(theDF, c(2), max)
      best_ones <- theDF[theDF[,"ACC"] == y["ACC"], ]
      if(nrow(best_ones)>1){
        method <-getMean(best_ones)
      }else{
        method <-row.names(best_ones)
      }
    }
    if(length((method) >1))
    {
      method <- method[1]
    }
    return(method)
  }

  all.attributes <- names(est.old)
  cnames <-  all.attributes[-charmatch(c("bug","bugdensity"),all.attributes)]
  
  
  theDF <-NULL
  data.eval.supervised <-data.eval
  data.eval.supervised$ebug.old <-data.eval$fbug.old
  data.eval.supervised$churn.est.old <-data.eval$churn.fit.old
  for( i in seq(cnames)){
    method <- paste(cnames[i],"UP", "supervised",sep=".")
    result <- subFunc(method=method, data.train=fit.old, data.valid=fit.old, data.eval=data.eval.supervised,xname=cnames[i])
    result <-as.data.frame(result)
    row.names(result) <- cnames[i]
    theDF <-rbind(theDF,result)
  }

  bestName <- findFeature(theDF)
  method <-paste(bestName,"UP", "_Sup",sep=".")
  # browser()
  result <-subFunc(method=method, data.train=fit.old, data.valid=est.old, data.eval=data.eval,xname=bestName)
  Popt <-c(Popt, result$Popt)
  ACC <-c(ACC, result$ACC)
  Prec <- c(Prec, result$Prec)
  F1 <- c(F1, result$F1)

  #######################################################################
  ### -> Output

  return(list(Popt=Popt, ACC=ACC, Prec=Prec, F1=F1))
}