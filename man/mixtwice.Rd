\name{mixtwice}
\alias{mixtwice}
\title{Large-scale hypothesis testing by variance mixing}
\description{
MixTwice deploys large-scale hypothesis testing in the case when testing units provide estimated effects and
estimated standard errors.  It produces empirical Bayesian local false discovery and sign rates for tests of zero effect.
}
\usage{
mixtwice(thetaHat, s2, Btheta = 15, Bsigma2 = 10, df, 
method = c("EM-pava", "AugLag"), maxit = 100, prop = 1)
}
\arguments{
  \item{thetaHat}{Estimated effect sizes (vector over testing units))}
  \item{s2}{Estimated squared standard errors of thetaHat (vector over testing units)}
  \item{Btheta}{Grid size parameter for effect distribution}
  \item{Bsigma2}{Grid size parameter for variance distribution}
  \item{df}{Degrees of freedom in chisquare associated with  estimated standard error}
  \item{method}{Method used for solving the non-parametric MLE optimization. \code{method = "EM-pava"} solves the optimization problem using EM approach with pool adjacent violator algorithm (pava) and \code{method = "EM-pava"} solves the optimization problem directly using Augmented Lagrangian approach.}
  \item{maxit}{Number of iterations in EM if \code{method = "EM-pava"}, with default, \code{maxit = 100}}
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
arrays by variance mixing}. Bioinformatics, 2021.
}
\author{
Zihao Zheng,  Michael A.Newton
}
\note{
See Zheng et al. 2021 for further detail
}
\seealso{
See Also as \code{?peptide_data}
}
\examples{
set.seed(1)

l = 1000 ## number of testing units

neach = 20 ## number of subjects in each group

pi0 = 0.8 ## null proportion

signal1 = rep(0, l)

signal2 = signal1

signal2[1:round((1-pi0)*l)] = rnorm(round((1-pi0)*l), mean = 0, sd = 3)

## I will generate the sigma^2 parameter

sigma2 = rep(1,l)

sigma = sqrt(sigma2)

## Then I can generate data

data1 <- data2 <- matrix(NA, nrow = l, ncol = neach)

for (i in 1:l) {
  
  data1[i,] = rnorm(neach, mean = rep(signal1[i], each = neach), sd = sigma[i])
  data2[i,] = rnorm(neach, mean = rep(signal2[i], each = neach), sd = sigma[i])
  
}

thetaHat = rowMeans(data2) - rowMeans(data1)

sd1 = apply(data1, 1, sd)
sd2 = apply(data2, 1, sd)

s2 = sd1^2/neach + sd2^2/neach

fit.EM = mixtwice(thetaHat, s2, Btheta = 15, Bsigma2 = 10, df = 2*neach - 2,
                   method = "EM-pava", maxit = 100, prop = 1)

## you can try to visualize the result

plot(fit.EM$grid.theta, cumsum(fit.EM$mix.theta), type = "s",
     xlab = "grid.theta", ylab = "ecdf of theta", lwd = 2)

lines(ecdf(signal2 - signal1),cex = 0.1, lwd = 0.5, lty = 2, col = "red")

legend("topleft", legend = c("fit.mix", "true.mix"), col = c("black", "red"), 
       lwd = 1, pch = 19)

## calculate false discovery rate and true positive under 0.05

oo = order(fit.EM$lfdr)

# number of discovery

x1 = sum(cumsum(fit.EM$lfdr[oo])/c(1:l) <= 0.05)

# number of true discovery

x2 = sum(cumsum(fit.EM$lfdr[oo])/c(1:l) <= 0.05 & (signal2 != 0)[oo])

# number of real positive

x3 = sum(signal2 != 0) 

# number of false discovery

x4 = sum(cumsum(fit.EM$lfdr[oo])/c(1:l) <= 0.05 & (signal2 == 0)[oo]) 

x4/x1 ## false discovery rate FDR

x2/x3 ## true positive rate

## null proportion estimation

max(fit.EM$mix.theta)

## you can also try using another method (Augmented Lagrangian), the result would be similar

# fit.AugLag = mixtwice2(thetaHat, s2, Btheta = 15, Bsigma2 = 10, df = 2*neach - 2,
#                        method = "AugLag", prop = 1)
}
\keyword{High density peptide array}
\keyword{Large-scale hypothesis testing}
\keyword{Mixing distribution}
\keyword{Local false discovery rate}
