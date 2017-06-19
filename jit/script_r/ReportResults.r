
require(effsize)

projects <- c("bugzilla", "platform","mozilla","jdt","columba","postgres")
# projects <- c("columba")   ################################### <<<<<<<<----------------------changed!
GetMagnitude <- function(str) {
    lbl <- NULL
    if (str == "large") {
        lbl <- "L"
    } else if (str == "medium") {
        lbl <- "M"
    } else if (str == "small") {
        lbl <- "S"
    } else if (str == "negligible") {
        lbl <- "T"
    } else {
        stop("error")
    }
    return(lbl)
}

GetOverall <- function(data, xnames, base, fname, pth=0.05, mth=c("negligible")) {

    tt.pvs <- wt.pvs <- NULL
    for (i in seq(xnames)) {
        treat <- data[, xnames[i]]

        tt.pvs <- c(tt.pvs, t.test(treat, base, paired=TRUE)$"p.value")
        apvs <- wilcox.test(treat, base, paired=TRUE)$"p.value"
        if (is.na(apvs)) { apvs <- 1 }
        wt.pvs <- c(wt.pvs, apvs)
    }

    wt.pvs <- p.adjust(wt.pvs, method="BH")
    tt.pvs <- p.adjust(tt.pvs, method="BH")

    dt.out <- data.frame(matrix(0, nrow=length(xnames), ncol=12))
    rownames(dt.out) <- xnames
    colnames(dt.out) <- c("avg", "improve.avg", "cohend", "cohend.magnitude", "tt.p-value", "med", "improve.med", "cliffd", "cliffd.magnitude", "wt.p-value", "overall", "colsall")

    base.out <- c(mean(base), 0, NA, NA, NA, median(base), NA, NA, NA, NA, NA, NA)
    dt.out[, "tt.p-value"] <- tt.pvs
    dt.out[, "wt.p-value"] <- wt.pvs
    dt.out[, "med"] <- apply(data[, xnames], MARGIN=2, FUN=median)
    dt.out[, "avg"] <- apply(data[, xnames], MARGIN=2, FUN=mean)

    colsall <- overall <- NULL
    cohend <- cliffd <- NULL
    cohend.magnitudes <- cliffd.magnitudes <- NULL
    for (i in seq(xnames)) {
        d2 <- data[, xnames[i]]

        cd <- effsize::cohen.d(d2, base, paired=TRUE)
        cohend <- c(cohend, cd$estimate)
        cohend.magnitudes <- c(cohend.magnitudes, cd$magnitude)

        cd <- effsize::cliff.delta(d2, base)
        cliffd <- c(cliffd, cd$estimate)
        cliffd.magnitudes <- c(cliffd.magnitudes, cd$magnitude)

        if (wt.pvs[i]<=pth) {
            if (median(d2, na.rm=TRUE) > median(base, na.rm=TRUE)) {
                if (cd$magnitude %in% mth) {
                    overall <- c(overall, "==")
                    colsall <- c(colsall, "black")
                } else {
                    overall <- c(overall, "vv")
                    colsall <- c(colsall, "blue")
                }
            } else {
                if (cd$magnitude %in% mth) {
                    overall <- c(overall, "==")
                    colsall <- c(colsall, "black")
                } else {
                    overall <- c(overall, "xx")
                    colsall <- c(colsall, "red")
                }
            }
        } else {
            overall <- c(overall, "==")
            colsall <- c(colsall, "black")
        }
    }

    return(data.frame(overall=overall, colsall=colsall, stringsAsFactors=FALSE))
}

Plot_SKtest.test <- function (validation="cv", criteria="ACC", types="UP", cnames=cnames,num_sup=4,all.projects=F) {

    require(ScottKnott)
    smnames <- cnames[1:num_sup]
    smlabels <- substr(cnames[1:num_sup],1,(nchar(cnames)-1-nchar(criteria)))
    uslabels <- substr(cnames[num_sup+1:length(cnames)],1,(nchar(cnames[num_sup+1:length(cnames)])-4-nchar(criteria)))
    # lbls.org <- c(lbnames, lbls.org)
    lbls.org <- c(smlabels, uslabels)
    xlabels <- c("EALR","RF","J48","IBk","NS","ND","NF","Entropy","LT","FIX","NDEV","AGE","EXP","REXP","SEXP","NUC")
    lbls.org <-xlabels
    cols.org <- c(rep("black", length(smnames)), rep("blue", length(uslabels)))

    df.1st <- NULL
    if (validation=="cv") {
        for (i in 1:length(projects)) {
            project <- projects[i]
            out.fname <- sprintf("../output/cross-validation/out_%s_%s.txt", project, criteria)
            rdata <- dget(out.fname)
            rdata <- rdata[, cnames]

            for (cname in cnames) {
                arun <- rdata[, cname]
                df.1st <- rbind(df.1st, data.frame(tr=xlabels, r=i, y=mean(rdata[, cname])))
            }
        }
    } else if (validation=="cv-tw") {
        # if(!all.projects){
          pdf(file=sprintf("results/SK-%s-%s-%s.pdf", validation, criteria,"ALL"), width=3, height=5,paper = "special")
          par(mfrow=c(6, 1), mai=c(0, 0, 0, 0), omi=c(0.15, 0.25, 0.30, 0), mex=1, cex=0.75)
        for (i in seq(length(projects))) {
            project <- projects[i]
            out.fname <- sprintf("../output/cross-validation-timewise/out_%s_%s.txt", project, criteria)
            rdata <- dget(out.fname)
            rdata <- rdata[, cnames]
            if(!all.projects){
              df.1st <-NULL
              for (cname in cnames) {
                arun <- rdata[, cname]
                df.1st <- rbind(df.1st, data.frame(tr=cname, r=i, y=(rdata[, cname])))
              }
              sk.1st <- with(df.1st, SK(x=df.1st, y=y, model='y ~ tr', which='tr'))
              lbls <- lbls.org[sk.1st$ord]
              cols <- cols.org[sk.1st$ord]
              
              plot(sk.1st, rl=FALSE, col=rep("black", length(lbls)), id.lab=rep("", length(lbls)), xlab=NA, ylab=NA, id.col=FALSE,xaxt="n", yaxt="n", main="", title="")###, col=rainbow(max(sk.1st$groups)), xlim=c(2, length(lbls)),

              box(bty="L")
              for (group in unique(sk.1st$groups)) {
                if (group == max(unique(sk.1st$groups))) {
                  next
                }
                abline(v=max(which(sk.1st$groups==group))+0.5, lty=2, col="black", lwd=1.5)
              }
              # mtext("Avg", side=2, line=2.5, adj=0.1, cex=0.6)
              axis(2, mgp=c(0.15, 0.15, 0), las=0, tck=-0.02, cex.axis=0.6, font=1,lwd=1,)
              text(x=seq(length(lbls)), y=par("usr")[3], srt=30, adj=c(1, 1.2), xpd=NA, cex=0.6, labels=lbls, col=cols)
              if(i==1){
                legend(x="top",ncol = 2, legend=expression(italic("Supervised"), italic("Unsupervised")), lty=c(1, 1), lwd=c(2, 2), col=c("black", "blue"), cex=0.8, merge=FALSE, bty="n")
              }
            }
            else{
              for (cname in cnames) {
                arun <- rdata[, cname]
                df.1st <- rbind(df.1st, data.frame(tr=cname, r=i, y=mean(rdata[, cname])))
              }
              
            }
        }
        dev.off()
    } else if (validation=="cp") {
        for(i in seq(length(projects))) {
            project1 <- projects[i]
            rdata <- NULL
            for(j in seq(length(projects))) {
                if(i != j) {
                    project2 <- projects[j]
                    out.fname <- sprintf("../output/cross-project/out_%s_cross_%s_%s.txt", project1, project2, criteria)
                    srdata <- dget(out.fname)

                    rdata <- rbind(rdata, srdata)
                }
            }

            rdata <- rdata[, cnames]
            for (cname in cnames) {
                df.1st <- rbind(df.1st, data.frame(tr=cname, r=i, y=mean(rdata[, cname])))
            }
        }
    }
    if(all.projects){
      sk.1st <- with(df.1st, SK(x=df.1st, y=y, model='y ~ tr', which='tr'))
      browser()
      lbls <- lbls.org[sk.1st$ord]
      cols <- cols.org[sk.1st$ord]
      
      png(file=sprintf("results/SK-%s-%s.png", validation, criteria), width=6, height=1.0, units="in", res=600, pointsize=10)
      par(mai=c(0, 0, 0, 0), omi=c(0.35, 0.3, 0.05, 0), mex=0.4, cex=1.0, cex.axis=0.7)
      
      plot(sk.1st, rl=FALSE, col=rep("black", length(lbls)), id.lab=rep("", length(lbls)), xlab=NA, ylab=NA, id.col=FALSE, xaxt="n", yaxt="n", main="", title="", mex=0.4, cex=0.6, cex.axis=0.7)###, col=rainbow(max(sk.1st$groups)), xlim=c(2, length(lbls)),
      
      box(bty="L")
      for (group in unique(sk.1st$groups)) {
        if (group == max(unique(sk.1st$groups))) {
          next
        }
        abline(v=max(which(sk.1st$groups==group))+0.5, lty=2, col="black", lwd=1.5)
      }
      # mtext("Avg", side=2, line=2.5, adj=0.5, cex=0.8)
      axis(2, mgp=c(0.15, 0.15, 0), las=0, tck=-0.02, cex.axis=0.6, font=1,lwd=1,)
      text(x=seq(length(lbls)), y=par("usr")[3], srt=60, adj=c(1, 1.2), xpd=NA, cex=0.6, labels=lbls, col=cols)
      
      legend(x="top", ncol=2, legend=expression(italic("Supervised"), italic("Unsupervised")), lty=c(1, 1), lwd=c(2, 2), col=c("black", "blue"), cex=0.8, merge=FALSE, bty="n")
      
      dev.off()
      
    }
    
}

plot_hist_single <- function(validation="cv", criteria="ACC", types="UP",all.projects=F,num_sup=12) {
    # smnames <- c("EALR", "NaiveBayes", "SimpleLogistic", "RBFNetwork", "SMO", "J48", "LMT", "RandomForest", "Ridor", "JRip", "IBk", "Bagging+LMT", "Bagging+NaiveBayes", "Bagging+SimpleLogistic", "Bagging+SMO", "Bagging+J48", "RotationForest+LMT", "RotationForest+NaiveBayes", "RotationForest+SimpleLogistic", "RotationForest+SMO", "RotationForest+J48", "AdaBoostM1+LMT", "AdaBoostM1+NaiveBayes", "AdaBoostM1+SimpleLogistic", "AdaBoostM1+SMO", "AdaBoostM1+J48", "RandomSubSpace+LMT", "RandomSubSpace+NaiveBayes", "RandomSubSpace+SimpleLogistic", "RandomSubSpace+SMO", "RandomSubSpace+J48")
    # lbnames <- c("EALR", "NB", "SL", "RBFN", "SMO", "J48", "LMT", "RF", "Ridor", "JRip", "IBk", "BG+LMT", "BG+NB", "BG+SL", "BG+SMO", "BG+J48", "RF+LMT", "RF+NB", "RF+SL", "RF+SMO", "RF+J48", "AB+LMT", "AB+NB", "AB+SL", "AB+SMO", "AB+J48", "RS+LMT", "RS+NB", "RS+SL", "RS+SMO", "RS+J48")
  smnames <- c("EALR","J48","J48.tuned","IBk","IBk.tuned","RandomForest","RandomForest.tuned")
  lbnames <- c("EALR","J48","J48.tuned","IBk","IBk.tuned","RandomForest","RandomForest.tuned")
  smnames <- c("EALR") ################################### <<<<<<<<----------------------changed!
  lbnames <- c("EALR","RF","J48","IBk","","NS","ND","NF","Entropy","LT","FIX","NDEV","AGE","EXP","REXP","SEXP","NUC","OneWay") ################################### <<<<<<<<----------------------changed!
    title <- NULL
    data <- data.all <- NULL
    if (validation=="cv") {
        for (project in projects) {
            out.fname <- sprintf("../output/cross-validation/out_%s_%s.txt", project, criteria)
            if (!file.exists(out.fname)) { next }
            rdata <- dget(out.fname)

            data.all[[project]] <- rdata[, cnames]
            data  <- rbind(data, rdata[, cnames])
        }
    } else if (validation=="cv-tw") {
        for (project in projects) {
            out.fname <- sprintf("../output/cross-validation-timewise/out_%s_%s.txt", project, criteria)
            if (!file.exists(out.fname)) { next }
            rdata <- dget(out.fname)

            data.all[[project]] <- rdata[, 5:length(rdata)]
            data  <- rbind(data, rdata[, 5:length(rdata)])
            # browser()
        }
    } else if (validation=="cp") {
        for (project1 in projects) {
            for (project2 in projects) {
                if (project1!=project2) {
                    out.fname <- sprintf("../output/cross-project/out_%s_cross_%s_%s.txt", project1, project2, criteria)
                    if (!file.exists(out.fname)) { next }
                    rdata <- dget(out.fname)

                    data <- rbind(data, rdata[, cnames])
                }
            }
        }
    }
    umnames <-NULL
    cnames <- names(data)
    cnames <- cnames[c(5,6,7,8,9,10,11,12,14,15,16,17,18)]  ## get rid of npt
    # browser()
    smnames <- cnames[1:num_sup]
    smlabels <- substr(cnames[1:num_sup],1,(nchar(cnames)-1-nchar(criteria)))
    uslabels <- substr(cnames[num_sup+1:length(cnames)],1,(nchar(cnames[num_sup+1:length(cnames)])-4-nchar(criteria)))
    # xlabels <- c(smlabels, uslabels)
    xlabels <- c("EALR","RF","J48","IBk","NS","ND","NF","Entropy","LT","FIX","NDEV","AGE","EXP","REXP","SEXP","NUC","OneWay")
    # xlabels <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")
    xlabels <- c("NS","ND","NF","Entropy","LT","FIX","NDEV","AGE","EXP","REXP","SEXP","NUC","OneWay")
    ################ get colors for boxplot border ####################

    overall <- data.frame(matrix(0, nrow=length(cnames), ncol=length(projects)))
    colnames(overall) <- projects
    overall <- apply(overall, c(1, 2), as.character)
    colsall <- overall
    if (validation!="cp") {
        for (project in projects) {
            # browser()
            s.index <- which.max(apply(data.all[[project]][, cnames[length(cnames)], drop=FALSE], 2, median))
            ares <- GetOverall(data=data.all[[project]], xnames=cnames, base=data.all[[project]][, names(data)[18]], fname=sprintf("results/(%s)%s-all-out-%s.csv", validation, project, criteria))
            overall[, project] <- ares$overall
            colsall[, project] <- ares$colsall
        }
    }

    s.index <- which.max(apply(data[, smnames, drop=FALSE], 2, median))
    ares <- GetOverall(data=data, xnames=cnames, base=data[, smnames[s.index]], fname=sprintf("results/(%s)all-out-%s.csv", validation, criteria))
    overall <- cbind(overall, ares$overall)
    colsall <- cbind(colsall, ares$colsall)
    
    # browser()
    ################ get colors for boxplot filledin, if better than best1 ####################
    overallOO <- data.frame(matrix(0, nrow=length(cnames), ncol=length(projects)))
    colnames(overallOO) <- projects
    overallOO <- apply(overallOO, c(1, 2), as.character)
    colsallOO <- overallOO
    if (validation!="cp") {
      for (project in projects) {
        # browser()
        s.index <- which.max(apply(data.all[[project]][, smnames, drop=FALSE], 2, median))
        ares <- GetOverall(data=data.all[[project]], xnames=cnames, base=data.all[[project]][, cnames[length(cnames)]], fname=sprintf("results/(%s)%s-all-out-%s.csv", validation, project, criteria))
        overallOO[, project] <- ares$overall
        colsallOO[, project] <- ares$colsall
      }
    }

    box.col <- gsub("black|red",NA,colsallOO,fixed=FALSE)
    box.col <- gsub("blue","green",box.col)

    #### parameters for plots
    ##### plots for each project
    par(mfrow = c(3, 2))
    pdf(file=sprintf("results/hist-%s-%s-%s.pdf", validation,"Unsup-OneWay",criteria), width=4, height=4,paper = "special")
    par(mfrow=c(6, 1), mai=c(0, 0, 0, 0), omi=c(0.35, 0.25, 0.30, 0), mex=1, cex=0.75)
    text_count<-1
    if(!all.projects){
      for(project in projects ){
        len <- length(smnames)

        pos <- c(1:len, (len+2):(length(cnames)+1))
        # browser()

        ylab.names <- projects
        
        cols <- colsall[, project]
        if(project !="bugzilla" && (criteria %in% c("Precision","PREC","F1"))){
          if(project == "mozilla"){
            boxplot(data.all[[project]][, cnames], xlab=NA, ylab=NA, xaxt="n", yaxt="n", at=pos, border=cols,boxwex=0.8, frame=FALSE, outline=FALSE,ylim=c(0,0.24)) ### xlim=c(2, len+12),
            axis(2, mgp=c(0.5, 0.15, 0), las=0, tck=-0.02, cex.axis=1, font=1,lwd=1,)
          }else{
            boxplot(data.all[[project]][, cnames], xlab=NA, ylab=NA, xaxt="n", yaxt="n", at=pos, border=cols,boxwex=0.8, frame=FALSE, outline=FALSE,ylim=c(0,0.7)) ### xlim=c(2, len+12),
            axis(2, mgp=c(0.5, 0.15, 0), las=0, tck=-0.02, cex.axis=1, font=1,lwd=1,)
          }
        }else{
          boxplot(data.all[[project]][, cnames], xlab=NA, ylab=NA, xaxt="n", yaxt="n", at=pos, border=cols,boxwex=0.8, frame=FALSE, outline=FALSE,ylim=c(0,1.1)) ### xlim=c(2, len+12),
          axis(2, mgp=c(0.5, 0.15, 0), las=0, tck=-0.02, cex.axis=1, font=1,lwd=1,)
        }

        ### tck: length of tick marks ### las: vertical or horizontal

        box(bty="L", lwd=1)
        
        ### xpd=NA print text outside plot margin

        temp <- max(apply(data.all[[project]][, cnames[len+1], drop=FALSE], 2, median))[1]

        points(x=c(1, length(c(smnames, uslabels))), y=c(temp, temp), type="l", lty="dashed")
        abline(v=len+1,  lty="dashed", col="black")
        if(text_count ==6){
          text(x=pos, y=par("usr")[3], srt=35, adj=c(0.9, 1.6), labels=xlabels, font=1, xpd=NA, cex=1) ### the model names
          # text(x=0,y=sum(par("usr")[3:4])*4/2,str=90,adj=c(0, 1.1),labels=criteria, cex=0.5)
          # browser()
        }
        if(text_count ==1){
          text(x=c((1+len)/2, (length(cnames))/2+0.4+(1+len)/2), y=par("usr")[4], adj=c(0.5, -0.5), font=3, xpd=NA, labels=c("Unsupervised", "Proposed"), cex=1)
          if(criteria == "ACC")
          {criteria <- "Recall"}
          if(criteria == "PREC")
          {criteria <- "Precision"}

          text(x=length(cnames)/2, y=par("usr")[4]+0.25, adj=c(0.5, -0.5), font=3, xpd=NA, labels=criteria, cex=1) 
        }
        text_count <- text_count +1
       
      }  
    }
  else{
    #### plots for all projects in one
    len <- length(smnames)
    pos <- c(1:len, (len+2):(len+13))
    
    png(file=sprintf("results/hist-%s-%s.png", validation, criteria), width=6, height=1.0, units="in", res=600, pointsize=10)
    par(mfrow=c(1, 1), mai=c(0, 0, 0, 0), omi=c(0.35, 0.15, 0.15, 0), mex=0.4, cex=1.0)
    ylab.names <- projects
    
    cols <- colsall[, length(projects)+1]
    boxplot(data[, cnames], xlab=NA, ylab=NA, xaxt="n", yaxt="n", at=pos, border=cols, boxwex=0.6, frame=FALSE, outline=FALSE) ### xlim=c(2, len+12),
    
    ### tck: length of tick marks ### las: vertical or horizontal
    axis(2, mgp=c(0.5, 0.5, 0), las=0, tck=-0.02, cex.axis=0.8, lwd=1.5)
    box(bty="L", lwd=1.5)
    
    ### xpd=NA print text outside plot margin
    temp <- max(apply(data[, cnames[1:len], drop=FALSE], 2, median))[1]
    points(x=c(1, length(c(smnames, umnames))), y=c(temp, temp), type="l", lty="dashed")
    abline(v=len+1,  lty="dashed", col="black")
    
    text(x=pos, y=par("usr")[3], srt=60, adj=c(1, 1.2), labels=xlabels, xpd=NA, cex=0.6) ### the model names
    text(x=c((1+len)/2, len+6.5), y=par("usr")[4], adj=c(0.5, -0.5), font=3, xpd=NA, labels=c("Supervised", "Unsupervised"), cex=0.8)
    
    # dev.off()
    
  }
    
  dev.off()
  return(cnames)
}

All_simple_models.tbl <- function (validation="cv", criteria="ACC", types="UP",cnames=NULL) {
    smname <- cnames[length(cnames)]
    n.row <- NULL
    a.names <- NULL
    if (validation=="cp") {
        n.row <- 1
        a.names <- "ALL"
    } else {
        n.row <- length(projects)
        a.names <- projects
    }

    data.out.r <- matrix(0, nrow=n.row+3, ncol=length(cnames))
    data.out.r <- as.data.frame(data.out.r)
    rownames(data.out.r) <- c(a.names, "AVG", "WTL", "Improve.e")
    colnames(data.out.r) <- cnames

    data.cliffd.r <- matrix(0, nrow=n.row, ncol=length(cnames))
    data.cliffd.r <- as.data.frame(data.cliffd.r)
    rownames(data.cliffd.r) <- paste(a.names, "cliff", sep="-")
    colnames(data.cliffd.r) <- cnames

    data.all <- NULL
    if (validation=="cv") {
        for (i in seq(projects)) {
            project <- projects[i]
            out.fname <- sprintf("../output/cross-validation/out_%s_%s.txt", project, criteria)
            rdata <- dget(out.fname)

            data.all[[i]] <- rdata
            data.out.r[i, ] <- apply(rdata[, cnames], MARGIN=2, FUN=median, na.rm=FALSE)
        }
    } else if (validation=="cv-tw") {
        for (i in seq(length(projects))) {
            project <- projects[i]
            out.fname <- sprintf("../output/cross-validation-timewise/out_%s_%s.txt", project, criteria)
            rdata <- dget(out.fname)

            data.all[[i]] <- rdata
            data.out.r[i, ] <- apply(rdata[, cnames], MARGIN=2, FUN=median, na.rm=FALSE)
        }
    } else if (validation=="cp") {
        rdata <- NULL
        cnamed <- NULL
        for(i in seq(length(projects))) {
            for(j in seq(length(projects))) {
                if(i != j) {
                    project1 <- projects[i]
                    project2 <- projects[j]
                    out.fname <- sprintf("../output/cross-project/out_%s_cross_%s_%s.txt", project1, project2, criteria)
                    srdata <- dget(out.fname)

                    rdata <- rbind(rdata, srdata[, cnames])
                    cnamed <- rbind(cnamed, c(project1, project2))
                }
            }
        }

        cnamed <- rbind(cnamed, c(" ", " "))
        avg <- apply(rdata, MARGIN=2, FUN=mean, na.rm=FALSE)
        write.csv(cbind(cnamed, round(rbind(rdata, avg), 3)), file=sprintf("results/(cp)-(%s).csv", criteria))

        data.out.r[1, ] <- apply(rdata[, cnames], MARGIN=2, FUN=median, na.rm=FALSE)
        data.all[[1]] <- rdata
    }

    data.out.r["AVG", ]       <- colMeans(data.out.r[1:n.row, ], na.rm=TRUE)
    data.out.r["Improve.e", ] <- 100*(data.out.r["AVG", ]-data.out.r["AVG", 1])/data.out.r["AVG", 1]

    data.out.r <- round(data.out.r, 3)
    data.out.r <- apply(data.out.r, c(1, 2), as.character)

    data.out.cliff <- data.out.r
    data.out.cliff <- 0

    PVs <- matrix(0, nrow=n.row, ncol=length(cnames))
    for (i in seq(cnames)) {
        for(j in seq(n.row)) {
            rdata <- data.all[[j]]
            base <- rdata[, smname]
            d1 <- rdata[, cnames[i]]
            err_catch <- try(wt.r <- wilcox.test(d1, base, paired=TRUE))
            if (class(err_catch)=="try-error") {
                PVs[j, i] <- NA
            } else {
                PVs[j, i] <- wt.r$"p.value"
            }
        }
    }

    BH.PVs <- PVs
    for (i in ncol(PVs)) {
        BH.PVs[, i] <- p.adjust(p=PVs[, i], method="BH")
    } ### BH.PVs <- matrix(BH.PVs, nrow=nrow(PVs), ncol=ncol(PVs))

    for (i in seq(cnames)) {
        wtl <- 0
        for(j in seq(n.row)) {
            rdata <- data.all[[j]]
            base <- rdata[, smname]
            d1 <- rdata[, cnames[i]]

            wt.r <- BH.PVs[j, i]
            if (!is.na(wt.r)) {
                if(wt.r <= 0.05) {
                    if(median(d1, na.rm=TRUE) > median(base, na.rm=TRUE)) {
                        wtl <- wtl + 100  ### win
                        data.out.r[j, i] <- paste(data.out.r[j, i], "v", sep="")
                    } else {
                        wtl <- wtl + 1 ### loss
                        data.out.r[j, i] <- paste(data.out.r[j, i], "x", sep="")
                    }
                } else {
                    wtl <- wtl + 10   ### tie
                }
            }

            err_catch <- try(cd.r <- effsize::cliff.delta(d1, base))
            if (class(err_catch)=="try-error") {
                data.cliffd.r[j, i] <- "NA"
            } else {
                if (is.na(cd.r$magnitude) || is.na(cd.r$estimate)) {
                    data.cliffd.r[j, i] <- "NA"
                } else {
                    data.cliffd.r[j, i] <- paste(GetMagnitude(cd.r$magnitude), as.character(round(cd.r$estimate, 3)), sep="")
                }
            }
        }

        data.out.r["WTL", i] <- as.character(wtl)
    }

    if (validation=="cp") {
        data.out.r["AVG", ] <- data.out.r["ALL", ]
    }

    data.out.r <- rbind(data.out.r, data.cliffd.r)

    if (validation!="cp") {
        fname <- sprintf("results/(%s)all-out-%s-models(All).csv", validation, criteria)
        write.table(c("class"), file=fname, row.names=FALSE, col.names=FALSE, append=FALSE, eol=",")
        write.table(data.out.r,   file=fname, row.names=TRUE,  col.names=TRUE,  append=TRUE, sep=",")
    }
}

# validations <- c("cv", "cv-tw", "cp")
validations <- c("cv-tw")
for (validation in validations) {
    criterias <- c("Popt","ACC","F1","PREC")
    # criterias <- c("PREC")
    for (criteria in criterias) {
        cat(validation, criteria, "\n")
        cnames=plot_hist_single(validation=validation, criteria=criteria)
        All_simple_models.tbl(validation=validation, criteria=criteria,cnames=cnames)
        # Plot_SKtest.test(validation=validation, criteria=criteria, cnames=cnames)
    }
}
