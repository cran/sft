\name{capacityGroup}
\alias{capacityGroup}
\title{
Capacity Analysis
}
\description{
Performs workload capacity analysis on each participant and each condition.  Plots each capacity coefficient individually and returns the results of the nonparametric null-hypothesis test for unlimited capacity independent parallel performance.
}
\usage{
capacityGroup(inData, acc.cutoff=.9, ratio=TRUE, plotCt=TRUE, ...)
}
\arguments{
  \item{inData}{Data collected from a Double Factorial Paradigm experiment in standard form.}
  \item{acc.cutoff}{Minimum accuracy for each stimulus category used in calculating the capacity coefficient.}
  \item{ratio}{Indicates whether to return the standard ratio capacity coefficient or, if FALSE, the difference form.}
  \item{plotCt}{Indicates whether or not to generate plots of the capacity coefficients.}
  \item{...}{Arguments to be passed to plot function.}
}
\details{
For the details of the OR and AND capacity coefficients, see \code{\link{capacity.or}} and \code{\link{capacity.and}}.  If accuracy in any of the stimulus categories used to calculate a capacity coefficient falls below the cutoff, NA is retuned for that value in both the statistic and the Ct matrix.
}
\value{
  \item{statistic}{Data frame with columns indicating participant, condition, capacity statistic for OR processing and associated p values and the capacity statistic for AND processing and associated p values.}
  \item{Ct.or}{Matrix with each row giving the values of the of the estimated OR capacity coefficient for one participant in one condition for values of times.  The rows match the ordering of statistic.}
  \item{Var.or}{Matrix with each row giving the values of the of the variance of the estimated OR capacity coefficient for one participant in one condition for values of times.}
  \item{Ct.and}{Matrix with each row giving the values of the of the estimated AND capacity coefficient for one participant in one condition for values of times.  The rows match the ordering of statistic.}
  \item{Var.and}{Matrix with each row giving the values of the of the variance of the estimated AND capacity coefficient for one participant in one condition for values of times.}
  \item{times}{Times at which the matrix capacity coefficients are calculated in Ct.or and Ct.and.}
}
\references{
Townsend, J.T. & Wenger, M.J. (2004). A theory of interactive parallel processing: New capacity measures and predictions for a response time inequality series. \emph{Psychological Review, 111}, 1003--1035.

Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2012). Statistical Measures for Workload Capacity Analysis.  \emph{Journal of Mathematical Psychology, XX}, XXX-XXX.
}
\author{
Joe Houpt <jhoupt@indiana.edu>
}
\seealso{
\code{\link{capacity.and}}
\code{\link{capacity.or}}
}
\examples{
data(dots)
capacityGroup(dots)
}
\keyword{ sft }
