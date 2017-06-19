#setwd
# if(Sys.info()[1]=="Darwin"){
#   setwd("~/Github/R/jit/script_r")
# }else
# {
#   setwd("C:/Users/R/jit/script_r")
# }

### following RWeka packages can be downloaded in site: http://sourceforge.net/projects/weka/files/weka-packages/
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

source("core.r")

### 10 times 10-fold cross-validation
exeExperiment.crossvalidation <- function(project, sampling=TRUE, totalFolds=2, totalRuns=1) {
    Popt <- ACC <- NULL

    all.name <- paste(project, ".csv", sep="")
    file <- file.path("..", "input", all.name, fsep=.Platform$file.sep)
    if (!file.exists(file)) { return(NULL) }
    all.data <- read.csv(file=file, header=TRUE, row.names=1)

    cat("(BEG)Cross-validation", "for", project, "\n")

    theList <- NULL
    index <- 1
    for (run in seq(totalRuns)) {
        sub <- divideToKsubsets(data=all.data, k=totalFolds, seed=run)

        for (fold in seq(totalFolds)) {
            fit.name <- paste(project, "_fit_", fold, sep="")
            est.name <- paste(project, "_est_", fold, sep="")

            fit <- est <- NULL
            for (j in seq(totalFolds)) {
                if (j != fold) {
                    fit <- rbind(fit, sub[[j]])
                }
            }
            est <- sub[[fold]]
            tmp.res <- coreExperiment(fit=fit, est=est, sampling=sampling, seed=fold)
            theList[[index]] <- tmp.res
            index <- index + 1
        }
    }

    cat("(END)Cross-validation", "for", project, "\n")
    
    Popt.lst <- sapply(theList, '[', "Popt")
    ACC.lst <- sapply(theList, '[', "ACC")
    
    for (i in seq(totalRuns*totalFolds)) {
        Popt <- rbind(Popt, Popt.lst[[i]])
        ACC <- rbind(ACC, ACC.lst[[i]])
    }

    return(list(Popt=Popt, ACC=ACC))
}

### timewise cross-validation
exeExperiment.crossvalidation.timewise <- function(project, gap=2,interval=1, sampling=FALSE) {
    Popt <- ACC <- F1<- Prec<-NULL

    fname <- paste(project, ".csv", sep="")
    file <- file.path("..", "input", fname, fsep=.Platform$file.sep)
    if (!file.exists(file)) { return(NULL) }
    data <- read.csv(file=file, header=TRUE, row.names=1, as.is=TRUE)

    cat("(BEG)Timewise cross-validation for", project, "\n")

    data$commitTime <- strptime(data$commitdate, format="%Y/%m/%d %H:%M") #### data$timePOSIX  <- as.POSIXct(data$commitdate, format="%Y/%m/%d %H:%M") ### 
    data$commitTime <- strftime(data$commitTime, format="%Y/%m")
    data <- data[order(data$commitTime), ] ### 

    unimon <- unique(data$commitTime)
    unimon <- unimon[order(unimon)]

    totalFolds <- length(unimon)

    sub <- NULL ### dive data into totalFolds parts, each part corresponding to changes within one month
    for (fold in seq(totalFolds)) {
        sub[[fold]] <- data[which(data$commitTime==unimon[fold]), ]
    }

    cat("\n", project, "has", totalFolds, "folds, length of each fold:\n", unlist(lapply(sub, nrow)), "\n")

    theList <- NULL
    train.valid <- NULL
    list.index <- 1
    neg.pos <-NULL
    for (fold in seq(totalFolds)) {
        if (fold+2+interval+gap > totalFolds) { next }
        fit <- NULL
        for(i in seq(fold,fold+interval)){
          fit <- rbind(fit,sub[[i]])
        }
        est <- rbind(sub[[fold+1+interval+gap]], sub[[fold+2+interval+gap]])
        neg.pos <-c(neg.pos, (length(est$bug)-sum(est$bug))/sum(est$bug))
        fit$commitTime <- est$commitTime <- NULL
        theList[[list.index]] <- coreExperiment(fit=fit, est=est, sampling=sampling, seed=fold)
#         print(theList[[list.index]])
        list.index <- list.index + 1
        train.valid <- rbind(train.valid, c(fold, fold+2+gap))
        
    }
    # browser()
    cat("neg/pos: ",median(neg.pos),"\n")
    Popt.lst <- sapply(theList, '[', "Popt")
    ACC.lst <- sapply(theList, '[', "ACC")
    Prec.lst <- sapply(theList, '[', "Prec")
    F1.lst <- sapply(theList,'[',"F1")
    
    for (i in seq(length(Popt.lst))) {
        Popt <- rbind(Popt, Popt.lst[[i]])
        ACC <- rbind(ACC, ACC.lst[[i]])
        F1 <- rbind(F1, F1.lst[[i]])
        Prec <-rbind(Prec, Prec.lst[[i]])
    }

    colnames(train.valid) <- c("trainfold", "validfold")
    train.valid <- data.frame(train.valid)

    train.valid$trainmonth <- unimon[train.valid$trainfold]
    train.valid$validmonth <- unimon[train.valid$validfold]

    Popt <- cbind(train.valid, Popt)
    ACC <- cbind(train.valid, ACC)
    F1 <- cbind(train.valid, F1)
    Prec <- cbind(train.valid, Prec)
    
    cat("(END)Timewise cross-validation for", project, "\n")
    return(list(Popt=Popt, ACC=ACC, F1=F1, Prec=Prec))
}

### cross project validation
exeExperiment.crossproject <- function(project.fit, project.est, sampling=TRUE) {
    Popt <- ACC <- NULL

    cat("(BEG)Across-project", project.fit, ":", project.est, "\n")

    fit.name <- paste(project.fit,".csv",sep="")
    fit.data <- read.csv(file.path("..", "input", fit.name, fsep=.Platform$file.sep), header=TRUE, row.names=1)

    est.name <- paste(project.est,".csv",sep="")
    est.data <- read.csv(file.path("..", "input", est.name, fsep=.Platform$file.sep), header=TRUE, row.names=1)

    tmp.res <- NULL
    tmp.res <- coreExperiment(fit=fit.data, est=est.data, sampling=sampling, seed=0)
    
    cat("(END)Across-project", project.fit, ":", project.est, "\n")

    Popt <- rbind(Popt, tmp.res$Popt)
    ACC <- rbind(ACC, tmp.res$ACC)

    return(list(Popt=Popt, ACC=ACC))
}

printResult <- function(tmp, path, name) {
    out.fname <- paste(path, "out_", name, ".txt", sep="")
    dput(tmp, file=out.fname)

    out.fname <- paste(path, "out_", name, "_Popt.txt", sep="")
    dput(tmp$Popt, file=out.fname)

    out.fname <- paste(path, "out_", name, "_ACC.txt", sep="")
    dput(tmp$ACC, file=out.fname)
    
    out.fname <- paste(path, "out_", name, "_F1.txt", sep="")
    dput(tmp$F1, file=out.fname)
    
    out.fname <- paste(path, "out_", name, "_Prec.txt", sep="")
    dput(tmp$Prec, file=out.fname)
}

########################################
### output cross validation results
########################################
subFunc.cv <- function(project, sampling=TRUE) {
    path <- "../output/cross-validation/"
    if (!file.exists(file.path(path))) { dir.create(path) }
    file <- paste("../input/", project, ".csv", sep="")
    if (!file.exists(file)) { return(NULL) }
    file <- paste(path, "out_", project, "_ce.txt", sep="")
    if (file.exists(file)) { return(NULL) }
    tmp <- exeExperiment.crossvalidation(project=project, sampling=sampling)
    printResult(tmp, path, project)
}

########################################
### output timewise cross validation results
########################################
subFunc.tw <- function(project, sampling=T) {
    path <- "../output/cross-validation-timewise/"
    if (!file.exists(file.path(path))) { dir.create(path) }
    file <- paste("../input/", project, ".csv", sep="")
    if (!file.exists(file)) { return(NULL) }
    file <- paste(path, "out_", project, "_ce.txt", sep="")
    if (file.exists(file)) { return(NULL) }
    tmp <- exeExperiment.crossvalidation.timewise(project=project, sampling=sampling)
    printResult(tmp, path, project)
    
}

########################################
### output cross-project performance
########################################
subFunc.cp <- function(arg, sampling=TRUE) {
    path <- "../output/cross-project/"
    if (!file.exists(file.path(path))) { dir.create(path) }
    file1 <- paste("../input/", arg[1], ".csv", sep="")
    file2 <- paste("../input/", arg[2], ".csv", sep="")
    if (!(file.exists(file1))&&(file.exists(file2))) { return(NULL) }
    file <- paste(path, "out_", arg[1], "_cross_", arg[2], "_ce.txt", sep="")
    if (file.exists(file)) { return(NULL) }
    tmp <- exeExperiment.crossproject(project.fit=arg[1], project.est=arg[2], sampling=sampling)
    printResult(tmp, path, paste(arg[1], "_cross_", arg[2], sep=""))
}

source("utils.R")
# source("core.r")

projects <- c("bugzilla","columba","jdt","mozilla","postgres","platform")

print.time <- function(){
    this_time <- Sys.time()  ## Foo, print time for tacking.
    cat(as.character(this_time), "\n")
}

# #### cross-validation and time-wise-cross-validation
for (project in projects) {
    # print.time()
    # subFunc.cv(project)
    print.time()
    subFunc.tw(project)
}

# # #### across-project prediction
# args <- NULL; index <- 1
# for (i in seq(projects)) {
#     for (j in seq(projects)) {
#         if (i==j) { next }
#         args[[index]] <- c(projects[i], projects[j])
#         index <- index + 1
#     }
# }
# 
# for (arg in args) {
#     print.time()
#     subFunc.cp(arg)
# }

