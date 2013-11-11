fPCAcapacity <- function(sftData, dimensions, acc.cutoff=.75, OR = TRUE, ratio=TRUE, plotPCs=FALSE, ...) {
  subjects <- sort(unique(sftData$Subject))
  nsubjects <- length(subjects)
  conditions <- sort(unique(sftData$Condition))
  nconditions <- length(conditions)
  subj.out <- character()
  cond.out <- character()

  channels <- grep("Channel", names(sftData), value=T)
  nchannels <- length(channels)

  # Currently only does present versus absent
  #  To be implemented:  separate tests for each factorial
  #  salience condition;  Negative numbers for distractor
  for ( ch in channels ) {
      #sftData[,ch] <- sftData[,ch] != 0
    if(is.factor(sftData[,ch])) {
      sftData[,ch] <- as.numeric(levels(sftData[,ch]))[sftData[,ch]]
    }
    sftData <- subset(sftData, sftData[,ch] >=0)
  }

  tvec <- seq(quantile(sftData$RT,.001), quantile(sftData$RT,.999), 
              length.out=1000)# - median(sftData$RT)

  midpoint <- floor(length(tvec)/2)


  if (OR) {capacity <- capacity.or} else {capacity <- capacity.and }

  capAllMat <- numeric()
  varAllMat <- numeric()
  subjVec <- c()
  condVec <- c()

  allRT <- numeric()
  register <- numeric()
  good <- logical()
  RTlist <- vector("list", nchannels+1)
  CRlist <- vector("list", nchannels+1)

  #ltyvec <- rep(1:nsubjects, nconditions)
  #colvec <- rep(1:nconditions, each=nsubjects)
  ltyvec <- numeric()
  colvec <- numeric()
  condLegend <- levels(conditions)

  # Calculate capacity for each participant in each condition
  for ( cn in 1:nconditions ) {
    if (is.factor(conditions)) {cond <- levels(conditions)[cn]} else {cond <- conditions[cn] }
    condsubjects <- factor(with(sftData, sort(unique(Subject[Condition==cond]))))
    ncondsubjects <- length(condsubjects)
    for ( sn in 1:ncondsubjects ) {
      if (is.factor(condsubjects)) {subj <- levels(condsubjects)[sn]} else {subj <- condsubjects[sn] }

      subjVec <- c(subjVec, subj)
      condVec <- c(condVec, cond)

      ds <- sftData$Subject==subj & sftData$Condition==cond

      # Redundant Target Response Times
      usechannel <- ds & apply(sftData[,channels]>0, 1, all)
      RTlist[[1]] <- sftData$RT[usechannel]
      CRlist[[1]] <- sftData$Correct[usechannel]
      

      # Single Target Response Times
      for ( ch in 1:nchannels ) {
        usechannel <- ds & sftData[,channels[ch]]>0 & 
                      apply(as.matrix(sftData[,channels[-ch]]==0), 1, all)
        RTlist[[ch+1]] <- sftData$RT[usechannel]
        CRlist[[ch+1]] <- sftData$Correct[usechannel]
      }


      # Check to make sure accuracy on each condition is higher than acc.cutoff
      if(any(lapply(CRlist, mean)<acc.cutoff) | any(lapply(RTlist, length) < 10) ) {
        good <- c(good, FALSE)
        capAllMat <- rbind(capAllMat, rep(NA, length(tvec)))
        varAllMat <- rbind(varAllMat, rep(NA, length(tvec)))
        next
      } else{
        good <- c(good, TRUE)
      }

      # Tracks the amount of offset for each capacity function  (the median RT)
      register <- c(register, mean(median(RTlist[[1]], median(c(RTlist[2:nconditions],recursive=TRUE)))) )

      #  Shift response time distributions to register Ct curves
      #for ( ch in 1:(nchannels+1)) {
      #  RTlist[[ch]] <- RTlist[[ch]] - tail(register,1)
      #}

      capout <- capacity(RTlist, CRlist, ratio=TRUE)
      shiftn <- midpoint - max( which(tvec < tail(register,1)))

      subj.out <- c(subj.out, subj)
      cond.out <- c(cond.out, cond)
      ltyvec <- c(ltyvec, sn)
      colvec <- c(colvec, cn)

      if (ratio) {
        tmin <- max( c(lapply(RTlist, quantile, probs=c(.01)), recursive=TRUE), na.rm=TRUE)
        tmax <- min( c(lapply(RTlist, quantile, probs=c(.99)), recursive=TRUE), na.rm=TRUE)
        ct <- capout$Ct(tvec)
        #ct[tvec < tmin] <- mean(ct[tmin:(tmin+10)])
        #ct[tvec > tmax] <- mean(ct[(tmax+10):tmax])
        ct[tvec < tmin] <- NA
        ct[tvec > tmax] <- NA
        capAllMat <- rbind(capAllMat, shift(ct, shiftn))
      } else {
        varAllMat <- rbind(varAllMat, capout$Var(tvec))
        capAllMat <- rbind(capAllMat, shift(capout$Ct(tvec), shiftn))
      }
    }
   }

  tvec <- tvec - midpoint
  
  tmin <- min(tvec[!apply(is.na(capAllMat[good,]), 2, all)])
  tmax <- max(tvec[!apply(is.na(capAllMat[good,]), 2, all)])
  tgood <- tvec[tvec >= tmin & tvec <= tmax]
  capGoodMat <- capAllMat[good,tvec >= tmin & tvec <= tmax]
  k <- dim(capGoodMat)[1]


  if(plotPCs) {
      dev.new()
      par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
      matplot(tgood, t(capGoodMat), type='l', lty=ltyvec, col=colvec, 
              #xlim=xbound,#c(tmin, tmax+50),
              xlim=c(tmin,1300),# ylim=c(-4,4),
              main="Capacity", ylab="C(t)", xlab="Time (Adjusted)")
      if(nconditions <= 5) {
        legend("topright", legend=condLegend, lty=1, col=1:5, cex=.9)
      }
  }

  if(ratio) {
    # Replace NA values in each function with the average capacity across functions.
    capGoodmn <- apply(capGoodMat, 2, mean, na.rm=TRUE)
    for (i in 1:k) {
      capGoodMat[i, is.na(capGoodMat[i,])]  <- capGoodmn[is.na(capGoodMat[i,])]
    }
  } else {
    varGoodMat <- varAllMat[good,tvec >= tmin & tvec <= tmax] 
  }
  

  #  subtract mean (across participants and conditions) capacity function 
  capGoodmn <- apply(capGoodMat, 2, mean)
  capGoodMat <- capGoodMat - matrix(capGoodmn, k, length(tgood), byrow=T)

  
  if(plotPCs) {
      dev.new()
      par(mfrow=c(1,2), mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
      plot(c(tmin-1000, tmax+1000), c(0,0), type='l', 
              xlim=c(tmin, tmax), #ylim=c(-4,4),
              main="Mean C(t)", ylab="C(t)", xlab="Time (Adjusted)")
      lines(tgood, capGoodmn, lwd=2)
      if (ratio) { abline(0,0, lty=1, col=grey(.4)) }
      matplot(tgood, t(capGoodMat), type='l', lty=ltyvec, col=colvec, 
              #xlim=xbound,#c(tmin, tmax+50),
              xlim=c(tmin,tmax),# ylim=c(-4,4),
              main="C(t)-Mean C(t)", ylab="C(t)", xlab="Time (Adjusted)")
      if(nconditions <= 5) {
        legend("topright", legend=condLegend, lty=1, col=1:5, cex=.9)
      }
  }

  if (ratio) {
    wtvec <- rep(1, length(tgood))
  } else {
    wtvec <- apply(varGoodMat, 2, sum, na.rm=TRUE) 
    wtvec[wtvec<1E-4] <- 1E-4
    wtvec <- 1/wtvec
    wtvec[is.na(wtvec)] <- 0
    wtvec <- wtvec / sum(wtvec)

    if (OR) { 
      xbound <- c(min(tgood), min(tgood[which(wtvec < 1E-6)]))
    } else {
      xbound <- c(max(tgood[which(wtvec < 1E-6)]), max(tgood))
    }

    if(plotPCs) {
        dev.new()
        par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
        plot(tgood, wtvec, col='forestgreen', type='l',
            xlim=c(tmin, 1300),
            main="Weighting Function", xlab="Time (Adjusted)", ylab="")
    }
  }

  
  #capGoodMat <- t(capGoodMat)
  wtGoodMat <- t(capGoodMat) * matrix(wtvec, nrow=length(tgood), ncol=sum(good))

  basis <- create.bspline.basis(rangeval=c(min(tgood),max(tgood)), 
              nbasis=sum(good)-1, norder=4)
  capGoodfd <- smooth.basis(tgood, wtGoodMat, basis)
  pcastrGood <- pca.fd(capGoodfd$fd,dimensions)
  pcastrGoodVarmx <- varmx.pca.fd(pcastrGood)
 
  if(plotPCs) {
    values <- pcastrGood$values
    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0))
    plot(1:5, values[1:5]/sum(values),
               xlim=c(1, 5), ylim=c(0,1), pch=19,
               main="Scree Plot", xlab="Eigenfunction", ylab="Variance Accounted For")
    lines(1:5, values[1:5]/sum(values))
  }
  

  harmmat <- eval.fd(tgood, pcastrGood$harmonics)
  harmmat <- harmmat / (wtvec %*% matrix(1, 1, dimensions))
  facmult <- apply(abs(pcastrGood$scores), 2, mean)

  harmmatV <- eval.fd(tgood, pcastrGoodVarmx$harmonics)
  harmmatV <- harmmatV / (wtvec %*% matrix(1, 1, dimensions))
  facmultV <- apply(abs(pcastrGoodVarmx$scores), 2, mean)

  scoreout <- data.frame(subjVec,condVec)
  for ( i in 1:dimensions) {
      scoreout[[i+2]] <- rep(NA, length(scoreout[[1]]))
      scoreout[[i+2]][good] <- pcastrGood$scores[,i]
  }
  names(scoreout) <- c("Subject","Condition",paste("D",1:dimensions,sep=""))

  scoreoutV <- data.frame(subjVec,condVec)
  for ( i in 1:dimensions) {
      scoreoutV[[i+2]] <- rep(NA, length(scoreoutV[[1]]))
      scoreoutV[[i+2]][good] <- pcastrGoodVarmx$scores[,i]
  }
  names(scoreoutV) <- c("Subject","Condition",paste("D",1:dimensions,sep=""))

  pflist <- vector("list", length=dimensions)
  for (ifac in 1:dimensions) {
    pflist[[ifac]] <- approxfun(tgood,harmmatV[,ifac])
  }

  if(plotPCs) {
    if (ratio) { ylim<-c(0,2*mean(capGoodmn)) } else { ylim=c(-1,1) }
    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0), 
            mfrow=c(dimensions,3))
    for ( ifac in 1:dimensions) {
        mainstr <- paste("PC", ifac, "-", floor(100*pcastrGood$varprop[ifac]), "%")

        Wveci <- capGoodmn + facmult[ifac]* harmmat[,ifac]

        plot(tgood, Wveci, type='l', lty=2, main="", xlab="Time (Adjusted)", ylab="",
            ylim=ylim, xlim=c(tmin, tmax))
        lines(tgood, capGoodmn)
        abline(0,0, col=grey(.4))
        mtext(mainstr, side=2, line=1)

        if(ifac==1) {
            mtext("Component Function", side=3, line=.5)
            legend("topright", c("Component", "Mean"), lty=c(2,1))
        }

        plot(tgood, Wveci - capGoodmn, type='l', main="", xlab="Time (Adjusted)", ylab="",
            ylim=ylim, xlim=c(tmin, tmax))
        abline(0,0, col=grey(.4))
        if(ifac==1) {mtext("Component - Mean", side=3, line=.5)}

        plot(scoreout$Subject, scoreout[[ifac+2]], type="n", #ylim=c(-2,2),
          xaxt='n', ylab="", xlab="Subject")
          axis(1,at=1:10, labels=rep("",10), las=0, cex=.1, tck=-.02)
          mtext(side=1, 1:10, at=1:10, line=.05, cex=.7)
        text(scoreout$Subject, scoreout[[ifac+2]], labels=scoreout$Condition,
          col=colvec)
        if(ifac==1) {mtext("Score", side=3, line=.5)}
    }

    dev.new()
    par(mar=c(3.1, 3.1, 2.1, 1.1), mgp=c(1.75, .25,0), mfrow=c(dimensions,3))
    for ( ifac in 1:dimensions) {
        mainstr <- paste("PC", ifac, "-", floor(100*pcastrGoodVarmx$varprop[ifac]), "%")

        Wveci <- capGoodmn + facmultV[ifac]* harmmatV[,ifac]

        plot(tgood, Wveci, type='l', lty=2, main="", xlab="Time (Adjusted)", ylab="",
            ylim=ylim, xlim=c(tmin, tmax))
        lines(tgood, capGoodmn)
        abline(0,0, col=grey(.4))
        mtext(mainstr, side=2, line=1)

        if(ifac==1) {
            mtext("Component Function", side=3, line=.5)
            legend("topright", c("Component", "Mean"), lty=c(2,1))
        }

        plot(tgood, Wveci - capGoodmn, type='l', main="", xlab="Time (Adjusted)", ylab="",
            ylim=ylim, xlim=c(tmin, tmax))
        abline(0,0, col=grey(.4))
        if(ifac==1) {mtext("Component - Mean", side=3, line=.5)}

        plot(scoreout$Subject, scoreoutV[[ifac+2]], type="n", #ylim=c(-2,2),
          xaxt='n', ylab="", xlab="Subject")
          axis(1,at=1:10, labels=rep("",10), las=0, cex=.1, tck=-.02)
          mtext(side=1, 1:10, at=1:10, line=.05, cex=.7)
        text(scoreout$Subject, scoreoutV[[ifac+2]], labels=scoreout$Condition,
          col=colvec)
        if(ifac==1) {mtext("Score", side=3, line=.5)}
    }
  }
  
  return(list(Scores=scoreoutV, MeanCT=approxfun(tgood,capGoodmn), PF=pflist, medianRT=register))
}



shift <- function(x, n, wrap=FALSE) {
  # Shift an array (x) by n
  #  positive n shift right; negative n shfit left

  if (abs(n) > length(x) ) {
    if (!wrap ) { return( rep(NA, length(x))) }
    n <- n %% length(x)
  }

  if ( n >= 0 ) {
    s  <- length(x)-n +1
    if (wrap) {
      xout <- c( x[s:length(x)], x[1:(s-1)])
    } else {
      xout <- c(rep(NA,n), x[1:(s-1)])
    }
  } else {
    s <- abs(n)+1
    if (wrap) {
      xout <- c( x[s:length(x)], x[1:(s-1)])
    } else {
      xout <- c( x[s:length(x)], rep(NA, abs(n)))
    }
  }
  return(xout)
}
