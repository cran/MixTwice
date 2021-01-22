\name{mixtwice}
\alias{mixtwice}
\title{Large-scale hypothesis testing by variance mixing}
\description{
MixTwice deploys large-scale hypothesis testing in the case when testing units provide estimated effects and
estimated standard errors.  It produces empirical Bayesian local false discovery and sign rates for tests of zero effect.
}
\usage{
mixtwice(thetaHat, s2, Btheta = 15, Bsigma2 = 10, df, prop = 1)
}
\arguments{
  \item{thetaHat}{Estimated effect sizes (vector over testing units))}
  \item{s2}{Estimated squared standard errors of thetaHat (vector over testing units)}
  \item{Btheta}{Grid size parameter for effect distribution}
  \item{Bsigma2}{Grid size parameter for variance distribution}
  \item{df}{Degrees of freedom in chisquare associated with  estimated standard error}
  \item{prop}{Proportion of units randomly selected to fit the distribution, with default, \code{prop = 1} (use all units to fit the distribution).}
}
\details{
\code{mixtwice} takes estimated effects and standard errors over a set of testing units. To compute local error-rate statistics, it finds nonparametric MLEs of the underlying distributions. It is similar to "ashr", except that mixtwice allows both a mixing distribution of underlying effects theta as well as a mixing distribution over underlying variance parameters.  Furthermore, it treats the effect mixing distribution nonparametrically, but enforces the shape constraint that this distribution is unimodal with mode at theta=0.  (We do not assume symmetry). The distribution of variance parameters is also treated nonparametrically, but with no shape constraints. The observations are assumed to be emitted from a normal distribution (on estimated effects) and an independent Chi-square distribution (on estimated squared standard errors). 
}
\value{
\item{grid.theta}{Support of the estimated mixing distribution on effects}
\item{grid.sigma2}{Support of the estimated mixing distribution of variances}
\item{mix.theta}{Estimated distribution of effect size, on grid.theta}
\item{mix.sigma2}{Estimated distribution of variance, on grid.sigma2}
\item{lfdr}{Local false discovery rate for each testing unit}
\item{lfsr}{Local false sign rate for each testing unit}
}
\references{
Zheng et al. \emph{MixTwice: Large scale hypothesis testing for peptide 
arrays by variance mixing}. Technical Report, October 2020.
}
\author{
Zihao Zheng,  Michael A.Newton
}
\note{
cite the biorxiv/arxiv paper
}

\seealso{
See Also as \code{?peptide_data}
}
\examples{

### for a single group testing problem for zero effect

### basic setting, take p = 100 as an toy example

pi = 0.5 ## true value of null proportion

p = 100; n = 10 ## number of testing units and sample size of each unit

p1=(1-pi)*p; p2=pi*p ## number of non-null and null

mu=c(rnorm(round(p1), mean=0, sd=1), rep(0, round(p2)))

sd=rep(1, p)

x=NULL

for (i in 1:(p1+p2)) {
  
  xx=rnorm(n, mu[i], sd[i])
  
  x=rbind(x,xx)
}


thetaHat = rowMeans(x) ## effect size of each testing unit

s2 = apply(x, 1, sd)^2/n ## estimated variance of effect size

mm1=mixtwice(thetaHat=thetaHat, s2=s2, Btheta = 15, Bsigma2 = 10, df=n-1)

## summarize and visualize the result

# estimated mixing distribution and true mixing distribution
plot(mm1$grid.theta, cumsum(mm1$mix.theta), type = "s", 
     xlab = "grid.theta", ylab = "ecdf of theta", lwd = 2)
lines(ecdf(mu),cex = 0.1, lwd = 0.5, lty = 2, col = "red")

plot(mm1$grid.sigma2, cumsum(mm1$mix.sigma2), type = "s",
     xlab = "grid.sigma2", ylab = "ecdf of sigma2", lwd = 2)

# effect size and estimated local false discovery rate

plot(thetaHat, mm1$lfdr, pch = ".", cex = 2, xlab = "estimated effect size", ylab = "lfdr")

# true positive and false positive rate, under the level of 0.05

mean(mm1$lfsr[1:50]<=0.05) # true positive rate
mean(mm1$lfsr[51:100]<=0.05) # false positive rate

# null proportion estimation

max(mm1$mix.theta)

}
\keyword{High density peptide array}
\keyword{Large-scale hypothesis testing}
\keyword{Mixing distribution}
\keyword{Local false discovery rate}