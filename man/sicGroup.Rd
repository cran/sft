\name{sicGroup}
\alias{sicGroup}
\title{
SIC Analysis for a Group
}
\description{
Calculates the SIC for each individual in each condition of a DFP experiment.  The function will plot each individuals SIC and return the results of the test for stochastic dominance and the statistical test of SIC form.
}
\usage{
sicGroup(inData, sictest = "ks", domtest = "ks", plotSIC=TRUE, ...)
}
\arguments{
  \item{inData}{Data collected from a Double Factorial Paradigm experiment in standard form.}
  \item{sictest}{Which type of hypothesis test to use for SIC form. \"ks\" is the only test currently implemented.}
  \item{domtest}{Which type of hypothesis test to use for testing stochastic dominance relations, either as series of KS tests (\"ks\") or the dominance test based on Dirichlet process priors (\"dp\").  DP not yet implemented.}
  \item{plotSIC}{Indicates whether or not to generate plots of the survivor interaction contrasts.}
  \item{...}{Arguments to be passed to plot function.}
}
\details{
See the help page for the \code{\link{sic}} function for details of the survivor interaction contrast.
}
\value{
  \item{statistic}{Data frame with columns indicating participant, condition, scaling factor for the test statistic, value of the positive one sided KS test statistic and associated p-value, value of the negative one-sided KS test statistic and associated p-value, p-values for the ART tests and a value indicating which model is predicted to produce the indicated pattern of statistics.}
  \item{SIC}{Matrix with each row giving the values of the of the estimated SIC for one participant in one condition for values of times.  The rows match the ordering of statistic.}
  \item{Dominance}{List with each element giving the results of the KS test for survivor orderings implied by selective influence.}
  \item{times}{Times at which the matrix capacity coefficients are calculated in Ct.or and Ct.and.}
}
\references{
Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2010). The statistical properties of the survivor interaction contrast. \emph{Journal of Mathematical Psychology, 54}, 446-453.

Heathcote, A., Brown, S.D., Wagenmakers, E-J. & Eidels, A. (2010) Distribution-free tests of stochastic dominance for small samples. \emph{Journal of Mathematical Psychology, 54}, 454-463. 
}
\author{
Joe Houpt <jhoupt@indiana.edu>
}
\seealso{
\code{\link{sic}}
\code{\link{capacityGroup}}
}
\examples{
data(dots)
sicGroup(dots)
}
\keyword{ ~sft }
