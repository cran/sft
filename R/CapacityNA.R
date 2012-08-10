capacityGroup <- function(inData, acc.cutoff=.9, ratio=TRUE, plotCt=TRUE, ...) {
  subjects <- sort(unique(inData$Subject))
  nsubjects <- length(subjects)

  conditions <- sort(unique(inData$Condition))
  nconditions <- length(conditions)

  channels <- grep("Channel", names(inData), value=T)
  nchannels <- length(channels)

  times <- sort(unique(round(inData$RT)))

  Zor <- numeric()
  Zand <- numeric()
  Zsubj <- character()
  Zcond <- character()
  capORMat <- numeric()
  capANDMat <- numeric()
  if(!ratio) {
    varORMat <- numeric()
    varANDMat <- numeric()
  }

  capORMat <- numeric()
  capANDMat <- numeric()

  devmatOR <- matrix(NA, nconditions, ceiling(nsubjects/9))
  devmatAND <- matrix(NA, nconditions, ceiling(nsubjects/9))

  RTlist <- vector("list", nchannels)
  CRlist <- vector("list", nchannels)

  for ( cn in 1:nconditions ) {
    m <- 0
    if (is.factor(conditions)) {cond <- levels(conditions)[cn]} else {cond <- conditions[cn] }
    #cond <- conditions[cn]
    condsubjects <- with(inData, sort(unique(Subject[Condition==cond])))
    ncondsubjects <- length(condsubjects)
    for ( sn in 1:ncondsubjects ) {
      if (is.factor(condsubjects)) {subj <- levels(condsubjects)[sn]} else {subj <- condsubjects[sn] }
      #subj <- condsubjects[sn]
      if ( sn %% 9 == 1 ) {
        m <- m+1
        if(plotCt) {
          dev.new()
          par(mfrow=c(3,3))
          devmatOR[cn,m] <- dev.cur()

          dev.new()
          par(mfrow=c(3,3))
          devmatAND[cn,m] <- dev.cur()
        }
      }

      ds <- inData$Subject==subj & inData$Condition==cond
      #if ( sum(ds) ==0 ) { next };
      good1 <- TRUE
      usechannel <- ds & apply(inData[,channels]>0, 1, all)

      RTlist[[1]] <- inData$RT[usechannel]
      CRlist[[1]] <- inData$Correct[usechannel]


      if(mean(CRlist[[1]]) < acc.cutoff | sum(CRlist[[1]]) < 10) {
        good1 <- FALSE 
      }

      for ( ch in 1:nchannels ) {
        usechannel <- ds & inData[,channels[ch]]>0 & 
                      apply(as.matrix(inData[,channels[-ch]]==0), 1, all)
        RTlist[[ch+1]] <- inData$RT[usechannel]
        CRlist[[ch+1]] <- inData$Correct[usechannel]
        if(mean(CRlist[[ch+1]]) < acc.cutoff | sum(CRlist[[ch+1]]) < 10) {
            good1 <- FALSE
        }
      }

      if( good1 ) {
        Zsubj <- c(Zsubj, subj)
        Zcond <- c(Zcond, cond)

        cap1or  <- capacity.or(RT=RTlist, CR=CRlist, ratio=ratio)
        capORMat <- rbind(capORMat, cap1or$Ct(times))
        if(!ratio) {
          varORMat <- rbind(varORMat, cap1or$Var(times))
        }
        Zor <- c(Zor,cap1or$statistic)

        cap1and <- capacity.and(RT=RTlist, CR=CRlist, ratio=ratio)
        capANDMat <- rbind(capANDMat, cap1and$Ct(times))
        if(!ratio) {
          varANDMat <- rbind(varANDMat, cap1and$Var(times))
        }
        Zand <- c(Zand,cap1and$statistic)

        if(plotCt) {
          dev.set(devmatOR[cn,m])
          plot(times, tail(capORMat,1), type='l',
              xlab="Time", ylab="C(t)",
              main=paste(cond, "\nParticipant ", subj, sep=""),...)
          if(ratio) {
            abline(1,0, lwd=2)
          } else {
            lines(times, tail(capORMat,1)+1.96*sqrt(tail(varORMat,1)), lty=2)
            lines(times, tail(capORMat,1)-1.96*sqrt(tail(varORMat,1)), lty=2)
            abline(0,0, lwd=2)
          } 
          

          dev.set(devmatAND[cn,m])
          plot(times, tail(capANDMat,1), type='l',
              xlab="Time", ylab="C(t)",
              main=paste(cond, "\nParticipant ", subj, sep=""),...)
          if(ratio) {abline(1,0, lwd=2)} else{
            lines(times, tail(capANDMat,1)+1.96*sqrt(tail(varANDMat,1)), lty=2)
            lines(times, tail(capANDMat,1)-1.96*sqrt(tail(varANDMat,1)), lty=2)
            abline(0,0, lwd=2)
          } 
        }

      } else {
        Zsubj <- c(Zsubj, subj)
        Zcond <- c(Zcond, cond)

        capORMat <- rbind(capORMat, rep(NA, length(times)) )
        if(!ratio){
          varORMat <- rbind(varORMat, rep(NA, length(times)) )
        }
        Zor <- c(Zor,NA)

        capANDMat <- rbind(capANDMat, rep(NA, length(times)) )
        if(!ratio){
          varANDMat <- rbind(varANDMat, rep(NA, length(times)) )
        }
        Zand <- c(Zand,NA)

        if(plotCt) {
          dev.set(devmatOR[cn,m])
          plot(c(min(times),max(times)), c(1,1), type='l', lwd=2,
              xlab="Time", ylab="C(t)",
              main=paste(cond, "\nParticipant ", subj, sep=""),...) 
          text(mean(c(max(times),min(times))),1.2,"Not enough correct.",col='red')

          dev.set(devmatAND[cn,m])
          plot(c(min(times),max(times)), c(1,1), type='l', lwd=2,
              xlab="Time", ylab="C(t)",
              main=paste(cond, "\nParticipant ", subj, sep=""),...) 
          text(mean(c(max(times),min(times))),1.2,"Not enough correct.",col='red')
        }
      }
    }

    if(plotCt) {
      dev.new()
      if(sum(Zcond==cond) > 1) {
        matplot(times, t(capORMat[Zcond==cond,]), type='l', lty=1,
          main=paste(cond,"OR Capacity",sep="\n"), xlab="Time",ylab="C(t)",...)
        if(ratio) {abline(1,0, lwd=2)} else{abline(0,0, lwd=2)} 

        dev.new()
        matplot(times, t(capANDMat[Zcond==cond,]), type='l', lty=1,
          main=paste(cond,"AND Capacity",sep="\n"), xlab="Time",ylab="C(t)",...)
        if(ratio) {abline(1,0, lwd=2)} else{abline(0,0, lwd=2)} 
      } else {
        plot(times, capORMat[Zcond==cond,], type='l', lty=1,
          main=paste(cond,"OR Capacity",sep="\n"), xlab="Time",ylab="C(t)",...)
        if(ratio) {abline(1,0, lwd=2)} else{abline(0,0, lwd=2)} 

        dev.new()
        plot(times, capANDMat[Zcond==cond,], type='l', lty=1,
          main=paste(cond,"AND Capacity",sep="\n"), xlab="Time",ylab="C(t)",...)
      }

    }
  }

  Z <- as.data.frame(list(Subject=Zsubj, Condition=Zcond, 
          Z.or=Zor, p.lower=pnorm(Zor), p.upper=1-pnorm(Zor),
          Z.and=Zand, p.lower=pnorm(Zand), p.upper=1-pnorm(Zand)))
  if(ratio){
    return(list(statistic=Z, Ct.or=capORMat, Ct.and=capANDMat, times=times))
  } else {
    return(list(statistic=Z, Ct.or=capORMat, Var.or=varORMat, Ct.and=capANDMat, Var.and=varANDMat, times=times))
  }
}


capacity.or <- function(RT, CR=NULL, ratio=TRUE) {
    if ( is.null(CR) | (length(CR) != length(RT)) ) {
      CR <- vector("list", length(RT))
      for( i in 1:length(RT) ) {
        CR[[i]] <- rep(1, length(RT[[i]]))
      }
    } 
    times <- sort(unique(c(RT, recursive=TRUE))) 
    ncond <- length(RT) - 1 

    # Find Nelson-Aalen Cumulative Hazard Estimates
    numer <- estimateNAH(RT=RT[[1]], CR=CR[[1]])
    denom <- estimateUCIPor(RT=RT[1+(1:ncond)], CR=CR[1+(1:ncond)])

    rmtest <- ucipTest(RT, CR, OR=TRUE)

    if (ratio) {
      C.or <- numer$H(times) / denom$H(times)

      C.or[is.nan(C.or)] <- NA
      C.or[is.infinite(C.or)] <- NA
      C.or <- approxfun(times, C.or)
      return( list(Ct=C.or, statistic=rmtest$statistic, p.val=rmtest$p.val) )
    } else {
      C.or <- numer$H(times) - denom$H(times)
      C.or <- approxfun(c(0,times), c(0,C.or))
      Var.or <- numer$Var(times) + denom$Var(times)
      Var.or <- approxfun(c(0,times), c(0,Var.or))
      return( list(Ct=C.or, Var=Var.or, statistic=rmtest$statistic, p.val=rmtest$p.val) )
    }
}


capacity.and <- function(RT, CR=NULL, ratio=TRUE) {
    if ( is.null(CR) | (length(CR) != length(RT)) ) {
      CR <- vector("list", length(RT))
      for( i in 1:length(RT) ) {
        CR[[i]] <- rep(1, length(RT[[i]]))
      }
    } 
    times <- sort(unique(c(RT, recursive=TRUE))) 

    ncond <- length(RT) - 1 

    rmtest <- ucipTest(RT, CR, OR=FALSE)

    # Find Nelson-Aalen Reverse Cumulative Hazard Estimates
    denom <- estimateNAK(RT[[1]], CR[[1]])
    numer <- estimateUCIPand(RT=RT[1+(1:ncond)], CR=CR[1+(1:ncond)])

    # Calculate the and capacity coefficient
    if (ratio) {
      C.and <- numer$K(times) / denom$K(times)
      C.and[is.nan(C.and)] <- NA
      C.and[is.infinite(C.and)] <- NA
      C.and <- approxfun(times, C.and)
      return( list(Ct=C.and, statistic=rmtest$statistic, p.val=rmtest$p.val) )
    } else {
      C.and <- denom$K(times) - numer$K(times) 
      C.and <- approxfun(c(times,Inf), c(C.and,0))
      Var.and <- numer$Var(times) + denom$Var(times)
      Var.and <- approxfun(c(times,Inf), c(Var.and,0))
      return( list(Ct=C.and, Var=Var.and, statistic=rmtest$statistic, p.val=rmtest$p.val) )
    }
}
