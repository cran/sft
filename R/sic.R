sicGroup <- function(inData, sictest="ks", domtest="ks", plotSIC=TRUE, ...) {

  subjects <- sort(unique(inData$Subject))
  nsubjects <- length(subjects)

  conditions <- sort(unique(inData$Condition))
  nconditions <- length(conditions)

  SICnames <- c("SerialOR", "ParallelAND", "ParallelOR", "Coactive",
                "ParallelAND", "SerialAND")

  times <- sort(unique(round(inData$RT)))

  n <- 0
  nc1 <- 1
  Dom <- vector("list")

  sicAllMat <- numeric()
  if( sictest=="ks") {
    KSsubj <- character()
    KScond <- character()

    KS <- numeric()
    KSp <- numeric()
    micstat <- numeric()
    micp <- numeric()
    KSwin <- character()
    N <- numeric()
  } else {
    cat("Only KS-SIC test is currently implemented.\n")
    return(NA)
  }
  
  for ( cn in 1:nconditions ) {
    if (is.factor(conditions)) {cond <- levels(conditions)[cn]} else {cond <- conditions[cn] }
    for ( sn in 1:nsubjects ) {
      if (is.factor(subjects)) {subj <- levels(subjects)[sn]} else {subj <- subjects[sn] }
      if (plotSIC & ( sn %% 9 == 1) ) {
        dev.new()
        par(mfrow=c(3,3))
      }

      KScond <- c(KScond, cond)
      KSsubj <- c(KSsubj, subj)

      n <- n+1
      HH <- with(inData, RT[Subject==subj & Condition==cond & Correct & 
                            Channel1==2 & Channel2==2] )
      HL <- with(inData, RT[Subject==subj & Condition==cond & Correct & 
                            Channel1==2 & Channel2==1] )
      LH <- with(inData, RT[Subject==subj & Condition==cond & Correct & 
                            Channel1==1 & Channel2==2] )
      LL <- with(inData, RT[Subject==subj & Condition==cond & Correct & 
                            Channel1==1 & Channel2==1] )

      if ( min( length(HH), length(HL), length(LH), length(LL)) > 10 ) {
        sicn <- sic(HH=HH, HL=HL, LH=LH, LL=LL)
        sicAllMat <- rbind(sicAllMat, sicn$SIC(times))
        N <- rbind(N, sicn$N)
        KS  <- rbind(KS,  sicn$Dvals[,1])
        KSp <- rbind(KSp, sicn$Dvals[,2])
        micstat <- c(micstat, sicn$MIC[[1]])
        micp <- c(micp, sicn$MIC[[2]])
        Dom[[n]] <- sicn$Dominance
        if (sicn$Dvals[1,2] < .05) {
          if (sicn$Dvals[2,2] < .05) {
            if (sicn$MIC[[2]] < .05) {
              KSwin <- c(KSwin, "Coactive")
            } else {
              KSwin <- c(KSwin, "SerialAND")
            }
          } else {
            KSwin <- c(KSwin, "ParallelOR")
          }
        } else {
          if (sicn$Dvals[2,2] < .05) {
            KSwin <- c(KSwin, "ParallelAND")
          } else {
            KSwin <- c(KSwin, "SerialOR")
          }
        }

        if(plotSIC) {
          plot(times, sicn$SIC(times), type='l',
            main=paste(cond, " Condition\nParticipant ", subj, sep=""), 
            xlab="Time",ylab="SIC(t)",...)
        }
      }
    }

    if(plotSIC) {
      dev.new()
      matplot(times, t(sicAllMat[KScond==cond,]),type='l',lty=1,
        main=paste(cond, " Condition", sep=""), 
        xlab="Time",ylab="SIC(t)",...)
    }
    nc1 <- n+1
  }
  colnames(KS) <- c("D+", "D-")
  colnames(KSp) <- c("D+", "D-")
  statistic <- as.data.frame(list(Subject=KSsubj, Condition=KScond, N=N,
      Dpositive= KS[,1], p.val.positive=KSp[,1], Dnegative=KS[,2], p.val.negative=KSp[,2], MIC=micstat, p.val.mic=micp, Model=KSwin))
  return(list(statistic=statistic, SIC=sicAllMat, Dominance=Dom, times=times))
}


sic <- function(HH, HL, LH, LL, sictest="ks", domtest="ks") {
    RTall <- sort(unique(c(HH, HL, LH, LL)))
    HH.ecdf <- ecdf(HH)
    HL.ecdf <- ecdf(HL)
    LH.ecdf <- ecdf(LH)
    LL.ecdf <- ecdf(LL)

    #if (domtest == "ks") {
    dominance<-rbind(c((ks.test(HH,HL,alternative="greater",exact=FALSE)$p<.05),
                      (ks.test(HH,LH,alternative="greater",exact=FALSE)$p<.05),
                      (ks.test(HL,LL,alternative="greater",exact=FALSE)$p<.05),
                      (ks.test(LH,LL,alternative="greater",exact=FALSE)$p<.05)), 
                     c((ks.test(HH,HL,alternative="less",exact=FALSE)$p<.05), 
                      (ks.test(HH,LH,alternative="less",exact=FALSE)$p<.05),
                      (ks.test(HL,LL,alternative="less",exact=FALSE)$p<.05),
                      (ks.test(LH,LL,alternative="less",exact=FALSE)$p<.05)))

    colnames(dominance) <- c("S.hh S.hl", "S.hh S.lh", "S.hl S.ll", "S.lh S.ll")
    rownames(dominance) <- c("<", ">")
    #} else if (domtest=="dp") {
    #  dominance <- c( DPdom(HH,HL)$test, DPdom(HH,LH)$test, 
    #                  DPdom(HL,LL)$test, DPdom(LH,LL)$test)
    #}

    N<-1/length(HH)+1/length(HL)+1/length(LH)+1/length(LL)
    N <- 1/N

    sicall <- LH.ecdf(RTall) + HL.ecdf(RTall) - HH.ecdf(RTall) - LL.ecdf(RTall)
    SIC <- stepfun(RTall, c(0,sicall))

    #if (sictest=="ks") {
      Dplus  <- max(0,sicall)
      Dminus <- abs(min(0,sicall))
      p.Dplus  <- exp(-2 * N * Dplus ^2 )
      p.Dminus <- exp(-2 * N * Dminus ^2 )
      Dvals <- cbind(c(Dplus,Dminus), c(p.Dplus, p.Dminus) )
      colnames(Dvals) <- c("statistic", "p.value")
      rownames(Dvals) <- c("D+", "D-")
      mic <- micTest(HH, HL, LH, LL, ART=FALSE)
      return(list(SIC=SIC, Dominance=dominance, Dvals=Dvals, MIC=mic, N=N))
    #}
}


micTest <- function(HH, HL, LH, LL, ART=TRUE) {
    statistic <- (mean(HL) - mean(HH))  - (mean(LL) - mean(LH))
    if (ART) {
        rtall <- c(HH, HL, LH, LL)
        n1 <- length(HH)
        n2 <- length(HL)
        n3 <- length(LH)
        n4 <- length(LL)
        h1 <- c(rep(1, n1+n2), rep(0,n3+n4))
        h2 <- c(rep(1, n1), rep(0, n2), rep(1,n3), rep(0, n4))
        
        mA0 <- sum( rtall * (1-h1) ) / sum(1-h1)
        mA1 <- sum( rtall * h1 ) / sum(h1)
        mB0 <- sum( rtall * (1-h2) ) / sum(1-h2)
        mB1 <- sum( rtall * h2 ) / sum(h2)
        rtall.m <- rtall - (1-h1)*mA0- h1*mA1 - (1-h2)*mB0 - h2*mB1 
        ranks <- rank(rtall.m, ties.method="average")

        #rval <- anova(lm(ranks ~ h1*h2))[3,4:5]
        rval <- summary(lm(ranks~h1*h2))[[4]][4,4]
        return (list(statistic=statistic, p.val=rval))
    }

    else {
        op <- options(contrasts=c("contr.helmert", "contr.poly"))

        h1vec <- c(rep(1, length(HH)), rep(1, length(HL)), 
                   rep(0, length(LH)), rep(0, length(LL)))
        h2vec <- c(rep(1, length(HH)), rep(0, length(HL)), 
                   rep(1, length(LH)), rep(0, length(LL)))
        allrt <- c(HH, HL, LH, LL)
        rval <- summary(lm(allrt ~h1vec * h2vec))[[4]][4,4]
        options(op)
        return(list(statistic=statistic, p.val=rval))
    }
}
