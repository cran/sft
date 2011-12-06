ucipTest <- function(RT, CR, OR=TRUE) {
    ncond <- length(RT) 
    allRT <- c(RT, recursive=TRUE)
    Nt <- length(allRT)
    if (OR) { 
        Yarr <- rep(0, Nt)
        Y <- vector("list", ncond) 
    } else { 
        Garr <- rep(0, Nt)
        G <- vector("list", ncond) 
    }
    
    S <- vector("list", ncond)

    #  Eliminate ties
    if (length(unique(allRT)) < Nt) {
        for ( i in 1:length(RT) ) {
            RT[[i]] <- RT[[i]] + rnorm(length(RT[[i]]), sd=1E-8)
        }
        allRT <- c(RT, recursive=TRUE)
    }
    tvec <- sort(allRT)
    
    for (j in 1:ncond ) { 
        RTx <- sort(RT[[j]],index.return=TRUE)
        RT[[j]] <- RTx$x
        CR[[j]] <- as.logical(CR[[j]])[RTx$ix]

        if (OR) {
            for (i in 1:Nt ) {Yarr[i] <- sum(RT[[j]] >= tvec[i]) }
            Y[[j]] <- stepfun(tvec, c(0,Yarr) )
        } else {
            for (i in 1:Nt ) {Garr[i] <- sum(RT[[j]] <= tvec[i]) }
            G[[j]] <- stepfun(tvec, c(Garr,0), right=TRUE )
        }
    }

    if (OR) {
        Yst <- rep(0, Nt) 
        for ( i in 2:length(Y) ) {
            Yst <- Yst + Y[[i]](tvec)
        }
        W <- Y[[1]](tvec)*(Yst) / (Y[[1]](tvec)+Yst)
        W[is.nan(W)] <- 0
        W <- stepfun(tvec, c(0,W))


        numer <- sum( W(RT[[1]][ CR[[1]] ]) / Y[[1]](RT[[1]][ CR[[1]] ]),na.rm=TRUE)
        for ( i in 2:ncond) {
           numer<-numer-sum(W(RT[[i]][CR[[i]]])/Y[[i]](RT[[i]][CR[[i]]]),na.rm=TRUE) 
        }
        
        denom <- 0
        for ( i in 1:ncond ) {
            denom <- denom + sum((W(RT[[i]][CR[[i]]])/Y[[i]](RT[[i]][CR[[i]]]))^2,na.rm=TRUE)
        }
        denom <- sqrt(denom)


    } else {
        Gst <- rep(0, Nt) 
        for ( i in 2:length(G) ) {
            Gst <- Gst + G[[i]](tvec)
        }
        W <- G[[1]](tvec)*(Gst) / (G[[1]](tvec)+Gst)
        W <- stepfun(tvec, c(W,0), right=TRUE)

        numer <- sum( W(RT[[1]][ CR[[1]] ]) / G[[1]](RT[[1]][ CR[[1]] ]),na.rm=TRUE)

        for ( i in 2:ncond) {
           numer<-numer-sum(W(RT[[i]][CR[[i]]])/G[[i]](RT[[i]][CR[[i]]]),na.rm=TRUE) 
        }

        
        denom <- 0
        for ( i in 1:ncond ) {
            denom <- denom + sum((W(RT[[i]][CR[[i]]])/G[[i]](RT[[i]][CR[[i]]]))^2,na.rm=TRUE)
        }
        denom <- sqrt(denom)
    }

    if (!OR) numer <- -1*numer
    return(list(statistic=numer/denom, p.val=c(pnorm(numer/denom),1-pnorm(numer/denom))))
}
