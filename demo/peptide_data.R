readline("Load the peptdie data")

data(peptide_data)

readline("Process the data to get some useful statistics")

## effect size for peptide

thetaHat = apply(peptide_data[1:8], 1, mean) - apply(peptide_data[9:16], 1, mean)

## squared standard error for peptide

sd1 = apply(peptide_data[,1:8], 1, sd)
sd0 = apply(peptide_data[,9:16], 1, sd)

s2 = sd1^2/8 + sd0^2/8

## z-score for peptide

tstat = thetaHat/sqrt(s2)

z = qnorm(pt(tstat, df = 14))

## p-value for peptide

p = 2*(1-pnorm(abs(z)))

readline("Apply MixTwice tool using thetaHat and s2 as input")

mm = mixtwice(thetaHat = thetaHat, s2 = s2, 
              Btheta = 15, Bsigma2 = 10, df = 14, prop = 0.01)

plot(mm$grid.theta, cumsum(mm$mix.theta), type = "s", 
     xlab = "grid.theta", ylab = "ecdf of theta", lwd = 2)
plot(mm$grid.sigma2, cumsum(mm$mix.sigma2), type = "s",
     xlab = "grid.sigma2", ylab = "ecdf of sigma2", lwd = 2)

summary(mm$lfsr);sum(mm$lfsr<=0.05)

readline("Other methods for large scale testing can also be applied")

readline("BH adjustment")

p.adj = p.adjust(p, method = "BH")
sum(p.adj<=0.05)

readline("Storey's q-value")

lfdr.qvalue = qvalue::qvalue(p)$lfdr
sum(lfdr.qvalue<=0.05)

readline("Efron's locFDR")

lfdr.Efron = fdrtool::fdrtool(p, statistic = "pvalue", plot = F)$lfdr
sum(lfdr.Efron<=0.05)

readline("Adaptive shrinkage (ASH)")

lfsr.ash = get_lfsr(ashr::ash(thetaHat, sqrt(s2), df = 14))
sum(lfsr.ash<=0.05)
