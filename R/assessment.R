assessment <- function(RT, CR, OR=c(TRUE, FALSE), correct=c(TRUE, FALSE), fast=c(TRUE, FALSE), detection=TRUE) {
  slow <- !fast
  incorrect <- !correct

  times <- sort(unique(unlist(RT)))
  p.correct <- unlist(lapply(CR, mean))
  
  if (OR) {
    if (detection) {
      if (correct & fast) {
        ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==1 ])
        ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==1 ])
        ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==1 ])

        channel1C  <- ecdf.channel1(times) * p.correct[2] * (1-p.correct[3])
        channel2C  <- ecdf.channel2(times) * p.correct[3] * (1-p.correct[2])
        channel12C <- ecdf.channel1(times) * p.correct[2] * (1-ecdf.channel2(times)) * p.correct[3]
        channel21C <- (1-ecdf.channel1(times)) * p.correct[2] * ecdf.channel2(times) * p.correct[3]
        channelCC <- ecdf.channel1(times) * p.correct[2] * ecdf.channel2(times) * p.correct[3]

        numer <- log(channel1C + channel2C + channel12C + channel21C + channelCC)
        denom <- log(ecdf.redundant(times) * p.correct[1])

        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Detection OR: Correct and Fast"
      }
      if (correct & slow) {
        ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==1 ])
        ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==1 ])
        ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==1 ])
        
        channel1C  <- (1-ecdf.channel1(times)) * p.correct[2] * (1-p.correct[3])
        channel2C  <- (1-ecdf.channel2(times)) * p.correct[3] * (1-p.correct[2])
        channelCC <-  (1-ecdf.channel1(times)) * p.correct[2] * (1-ecdf.channel2(times)) * p.correct[3]
        
        numer <- log(channel1C + channel2C + channelCC)
        denom <- log( (1-ecdf.redundant(times)) * p.correct[1] )

        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Detection OR: Correct and Slow"
      }
      if (incorrect & fast ) {
        ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==0 ])
        ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==0 ])
        ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==0 ])

        numer <- log( ecdf.channel1(times) * (1-p.correct[2])) + log( ecdf.channel2(times) * (1-p.correct[3]))
        denom <- log( ecdf.redundant(times) * (1-p.correct[1]) )

        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Detection OR: Incorrect and Fast"
      }
      if (incorrect & slow) {
        ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==0 ])
        ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==0 ])
        ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==0 ])

        channel1I  <- (1-ecdf.channel1(times)) * (1-p.correct[2]) * (1-p.correct[3])
        channel2I  <- (1-ecdf.channel2(times)) * (1-p.correct[2]) * (1-p.correct[3])
        channelII <-  (1-ecdf.channel1(times)) * (1-ecdf.channel2(times)) * (1-p.correct[2]) * (1-p.correct[3])

        numer <- log( channel1I + channel2I - channelII )
        denom <- log( (1-ecdf.redundant(times)) * (1-p.correct[1]) )

        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Detection OR: Incorrect and Slow"
      }
    } else {
      if (correct & fast) {
        ecdf.redundant <- ecdf(RT[[1]][CR[[1]]==1])
      
        G <- vector("list", length(RT)-1)
        for (i in 2:length(RT)) {
          g <- rep(0, length(times))
          for ( tval in RT[[i]][CR[[i]]==1] ) {
            idx <- which(times==tval)
            g[idx] <- sum(RT[[i]] > tval)
          }
          g <- g/(sum(CR[[i]])*length(RT[[-1*i+5]]))
          G[[i-1]] <- cumsum(g)
        }
      
        numer <- rep(0,length(times))
        for ( i in 2:length(RT) ) {
          numer <- numer + p.correct[i]*G[[i-1]]
        }
        numer <- log(numer)
        denom <- log(ecdf.redundant(times)*p.correct[1])
      
        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Discrimination OR: Correct and Fast"
      }
  
      if (correct & slow) {
      
        ecdf.redundant <- ecdf(RT[[1]][CR[[1]]==1])
      
        G <- vector("list", length(RT)-1)
        for (i in 2:length(RT)) {
          g <- rep(0, length(times))
          for ( tval in RT[[i]][CR[[i]]==1] ) {
            idx <- which(times==tval)
            g[idx] <- sum(RT[[i]] > tval)
          }
          g <- g/(sum(CR[[i]])*length(RT[[-1*i+5]]))
          G[[i-1]] <- rev(cumsum(rev(g)))
        }
      
        numer <- rep(0,length(times))
        for ( i in 2:length(RT) ) {
          numer <- numer + p.correct[i]*G[[i-1]]
        }
        numer <- log(numer)
        denom <- log( (1-ecdf.redundant(times))*p.correct[1] )
      
        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Discrimination OR: Correct and Slow"
      }
      
      if (incorrect & fast ) {
        p.incorrect <- 1-unlist(lapply(CR, mean))
      
        ecdf.redundant <- ecdf(RT[[1]][CR[[1]]==0])
      
        G <- vector("list", length(RT)-1)
        for (i in 2:length(RT)) {
          g <- rep(0, length(times))
          for ( tval in RT[[i]][CR[[i]]==0] ) {
            idx <- which(times==tval)
            g[idx] <- sum(RT[[i]] > tval)
          }
          g <- g/(sum(1-CR[[i]])*length(RT[[-1*i+5]]))
          G[[i-1]] <- cumsum(g)
        }
      
        numer <- rep(0,length(times))
        for ( i in 2:length(RT) ) {
          numer <- numer + p.incorrect[i]*G[[i-1]]
        }
        numer <- log(numer)
        denom <- log( ecdf.redundant(times)*p.incorrect[1] )
      
        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Discrimination OR: Incorrect and Fast"
      }
      
      if (incorrect & slow) {
        p.incorrect <- 1-unlist(lapply(CR, mean))
      
        ecdf.redundant <- ecdf(RT[[1]][CR[[1]]==0])
      
        G <- vector("list", length(RT)-1)
        for (i in 2:length(RT)) {
          g <- rep(0, length(times))
          for ( tval in RT[[i]][CR[[i]]==0] ) {
            idx <- which(times==tval)
            g[idx] <- sum(RT[[i]] > tval)
          }
          g <- g/(sum(1-CR[[i]])*length(RT[[-1*i+5]]))
          G[[i-1]] <- rev(cumsum(rev(g)))
        }
      
        numer <- rep(0,length(times))
        for ( i in 2:length(RT) ) {
          numer <- numer + p.incorrect[i]*G[[i-1]]
        }
        numer <- log(numer)
        denom <- log( (1-ecdf.redundant(times))*p.incorrect[1] )
      
        At <- numer/denom
        At[!is.finite(At)] <- NA
        A <- stepfun(times, c(NA,At))
        attributes(A)$call <- "Discrimination OR: Incorrect and Slow"
      }
    }
  } else {
    if (correct & fast) {
      ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==1 ])
      ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==1 ])
      ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==1 ])

      numer <- log( ecdf.channel1(times) * p.correct[2] ) + log( ecdf.channel2(times) * p.correct[3] ) 
      denom <- log( ecdf.redundant(times) * p.correct[1] ) 

      At <- numer/denom
      At[!is.finite(At)] <- NA
      A <- stepfun(times, c(NA,At))
      attributes(A)$call <- "AND: Correct and Fast"
    }
    if (correct & slow) {
      ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==1 ])
      ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==1 ])
      ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==1 ])

      channel1C  <- (1-ecdf.channel1(times)) * p.correct[2] * p.correct[3]
      channel2C  <- (1-ecdf.channel2(times)) * p.correct[3] * p.correct[2]
      channelCC <-  (1-ecdf.channel1(times)) * p.correct[2] * (1-ecdf.channel2(times)) * p.correct[3]

      numer <- log ( channel1C + channel2C - channelCC ) 
      denom <- log ( (1-ecdf.redundant(times)) * p.correct[1] ) 

      At <- numer/denom
      At[!is.finite(At)] <- NA
      A <- stepfun(times, c(NA,At))
      attributes(A)$call <- "AND: Correct and Slow"
    }
    if (incorrect & fast ) {
      ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==0 ])
      ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==0 ])
      ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==0 ])

      channel1I  <- ecdf.channel1(times) * (1-p.correct[2]) * p.correct[3]
      channel2I  <- ecdf.channel2(times) * (1-p.correct[3]) * p.correct[2]
      channel12I <- ecdf.channel1(times) * (1-p.correct[2]) * (1-ecdf.channel2(times)) * (1-p.correct[3])
      channel21I <- (1-ecdf.channel1(times)) * (1-p.correct[2]) * ecdf.channel2(times) * (1-p.correct[3])
      channelII <- ecdf.channel1(times) * (1-p.correct[2]) * ecdf.channel2(times) * (1-p.correct[3])

      numer <- log( channel1I + channel2I + channel12I + channel21I + channelII)
      denom <- log( ecdf.redundant(times) * (1-p.correct[1]) )

      At <- numer/denom
      At[!is.finite(At)] <- NA
      A <- stepfun(times, c(NA,At))
      attributes(A)$call <- "AND: Incorrect and Fast"

    }
    if (incorrect & slow) {
      ecdf.redundant <- ecdf(RT[[1]][ CR[[1]]==0 ])
      ecdf.channel1  <- ecdf(RT[[2]][ CR[[2]]==0 ])
      ecdf.channel2  <- ecdf(RT[[3]][ CR[[3]]==0 ])

      channel1I  <- (1-ecdf.channel1(times)) * (1-p.correct[2]) * p.correct[3]
      channel2I  <- (1-ecdf.channel2(times)) * (1-p.correct[3]) * p.correct[2]
      channelII <-  (1-ecdf.channel1(times)) * (1-p.correct[2]) * (1-ecdf.channel2(times)) * (1-p.correct[3])

      numer <- log( channel1I + channel2I + channelII )
      denom <- log( (1-ecdf.redundant(times)) * (1-p.correct[1]) )
      
      At <- numer/denom
      At[!is.finite(At)] <- NA
      A <- stepfun(times, c(NA,At))
      attributes(A)$call <- "AND: Incorrect and Slow"

    }
  }
  return(A)
}
