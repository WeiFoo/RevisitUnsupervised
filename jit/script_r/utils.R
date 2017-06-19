
divideToKsubsets <- function(data, k=10, seed, rnd=T){
	set.seed(seed=seed)

	if(rnd){data <- data[order((runif(nrow(data)))),]}

	numResidual <- k - nrow(data)%%k
	dummyData<-as.data.frame(matrix(NA,nrow=numResidual,ncol=ncol(data)))
	names(dummyData)<-names(data)

	data<-rbind(data,dummyData)
	splitData<-split(data,1:k)

	for(i in 1:k){
		splitData[[i]] <- na.omit(splitData[[i]])
	}

	return(splitData)
}

decChurn <- function(data){
	nf  <- data$nf             
	lt  <- data$lt * nf
	lt[lt==0] <- lt[lt==0] + 1
	churn <- ((data$la + data$ld) * lt)/2  # LA and LD was normalized by LT

	return (churn)
}

doSampling <- function(data, obj, seed=0){
	set.seed(seed=seed)   #### random seed make the experiment result to be replicated
	dataT <- data[data[, obj]==1, ]
	dataN <- data[data[, obj]==0, ]
	if (nrow(dataT) >= nrow(dataN)) {
		dataT <- dataT[order(runif(nrow(dataT))), ]
		dataT <- dataT[1:nrow(dataN), ]
	} else {
		dataN <- dataN[order(runif(nrow(dataN))), ]
		dataN <- dataN[1:nrow(dataT), ]
	}
	data <- rbind(dataT, dataN)

	return (data)
}

### this function sort the data frame base on the predicted value and second by other policies
SortData <- function(data, effortaware=FALSE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
{
	if (!all(c("NUM", "REL", "LOC", "PRE") %in% colnames(data))) { stop("ERROR: NUM, REL, LOC, or PRE is not colnames of data frame") }

	key.2nd <- "REL"
	if (effortaware) { key.2nd <- "density" }
	if (effortaware) {
		if (!(key.2nd %in% colnames(data))) {
			data[[key.2nd]] <- data$NUM/(data$LOC+1)
		}
		sorted <- FALSE
	}

	if (!sorted) {
		if (worstcase) {
			data <- data[order(-data$PRE, +data[[key.2nd]], -data$LOC), ]
		} else if (bestcase) {
			data <- data[order(-data$PRE, -data[[key.2nd]], +data$LOC), ]
		} else if (LOCUP) {
			data <- data[order(-data$PRE, +data$LOC), ]
		} else {
			data <- data[order(-data$PRE), ]
		}
		sorted <- TRUE
	}

	return(data)
}

### compute area under alberg curve
ComputeArea <- function(data, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
{
	if (!sorted) {
		data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
		sorted <- TRUE
	}

	len <- nrow(data)
	### cumulative sums of LOC or NUM
	cumXs <- cumsum(data$LOC) # x: LOC%
	cumYs <- cumsum(data$NUM) # y: Bug%

	Xs <- cumXs/cumXs[len]
	Ys <- cumYs/cumYs[len]

	fix_subareas <- vector(length=len)
	fix_subareas[1] <- 0.5 * Ys[1] * Xs[1]
	fix_subareas[2:len] <- 0.5 * (Ys[1:(len-1)] + Ys[2:len]) * abs(Xs[1:(len-1)] - Xs[2:len])

	area <- sum(fix_subareas)

	return(area)
}

ComputePopt <- function(data, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
{
	if (!sorted) {
		data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
		sorted <- TRUE
	}

	data.mdl <- data
	data.opt <- data[order(-data$density, +data$LOC), ]
	data.wst <- data[order(+data$density, -data$LOC), ]

	area.mdl <- ComputeArea(data=data.mdl, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
	area.opt <- ComputeArea(data=data.opt, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
	area.wst <- ComputeArea(data=data.wst, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)

	Popt <- 1 - (area.opt - area.mdl)/(area.opt - area.wst)

	return(Popt)
}

### compute recall with 20% effort
ComputeACC <- function(data, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
{
	if (!sorted) {
		data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
		sorted <- TRUE
	}

	len <- nrow(data)
	### cumulative sums of LOC or NUM
	cumXs <- cumsum(data$LOC) # x: LOC%
	cumYs <- cumsum(data$NUM) # y: Bug%

	Xs <- cumXs/cumXs[len]
	Ys <- cumYs/cumYs[len]
	pos <- min(which(Xs >= 0.2))
	ACC <- cumYs[pos] / cumYs[len]

	return(ACC)
}

ComputePrecF1 <- function(data, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE){
  if (!sorted) {
    data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sorted <- TRUE
  }
  
  len <- nrow(data)
  ### cumulative sums of LOC or NUM
  cumXs <- cumsum(data$LOC) # x: LOC%
  cumYs <- cumsum(data$NUM) # y: Bug%
  
  Xs <- cumXs/cumXs[len]
  Ys <- cumYs/cumYs[len]
  
  pos <- min(which(Xs >= 0.2))
  ACC <- cumYs[pos] / cumYs[len]
  Prec <- cumYs[pos] /pos
  F1 <- (ACC * Prec* 2 )/(Prec + ACC+0.000001)
  x <-list(Prec=Prec, F1=F1)
  return(x)
}

ComputeF <- function(data, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE){
  if (!sorted) {
    data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sorted <- TRUE
  }
  
  len <- nrow(data)
  ### cumulative sums of LOC or NUM
  cumXs <- cumsum(data$LOC) # x: LOC%
  cumYs <- cumsum(data$NUM) # y: Bug%
  
  Xs <- cumXs/cumXs[len]
  Ys <- cumYs/cumYs[len]
  
  pos <- min(which(Xs >= 0.2))
  ACC <- cumYs[pos] / cumYs[len]
  Precision <- cumsYs[pos] /pos
  
  F <- ACC * Precision* 2 /(Precision + ACC)
  return(F)
}

ComputePCA <- function(fit, select=3){
  # ### remove those variables with 0 variance
  va <- apply(fit, 2, var)
  fit <- fit[,((!is.na(va))&(va!=0))] ### fit <- fit[,!(va==0)]
  pca <- prcomp(fit, center = TRUE, scale. = TRUE)
  value <- NULL
  if(select==1){
    value <- as.data.frame(pca$x[,1])
    names(value) <- "PC1"
  }else{
    value <-pca$x[,1:select]
  }
  return (value)
  
}

ComputeSquare <- function(fit, original.attr=NULL){
  idx <- charmatch(c("fix","PC1","PC2","PC3","bug"), original.attr)
  fit <- fit[,-c(idx)]
  
  names <- paste(names(fit),"2", sep = "^")
  names(fit) <- names
  fit <-fit^2
  
  va <- apply(fit, 2, var)
  fit <- fit[,((!is.na(va))&(va!=0))] ### fit <- fit[,!(va==0)]
  return (fit)
}


ComputePower <- function(fit, original.attr=NULL, power=NULL){

  idx <- charmatch(c("PC1","PC2","PC3","bug"), original.attr)
  fit <- fit[,-c(idx)]
  
  names <- paste(names(fit),round(power,1), sep = "")
  names(fit) <- names
  fit <-fit^power
  
  # va <- apply(fit, 2, var)
  # fit <- fit[,((!is.na(va))&(va!=0))] ### fit <- fit[,!(va==0)]
  return (fit)
}

ComputeSum <- function(fit,original.attr){

  idx <- charmatch(c("bug","bugdensity"), original.attr)
  fit <- fit[,-c(idx)]
  fit.min <- apply(fit, 2, min)
  fit.max <- apply(fit, 2, max)
  doit <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -
                                                   min(x, na.rm=TRUE)+0.0001)}
  temp <- sapply(fit,doit)
  new.feature <- as.data.frame(rowSums(temp))
  names(new.feature) <- "C"
  
  return(new.feature)
  
}
