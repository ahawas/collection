require(plyr)
alpha = -0.0259
sig = 0.0407
lambda = 1.512
mu = 0.0792
si = 0.1699
#ns = lapply(rpois(100,lambda), function(x) {sum(rnorm(x,alpha, sig)) } )
N = 20000000
ns.samples = rep(NA, N)
for( i in 1:N){
  j = rpois ( 1, lambda)
  yj = sum(rnorm(j, alpha, sig))
  ns.samples[i] = yj+rnorm(1, mu, si)
}
pn_mu = mean(ns.samples)
pn_sig2 = var(ns.samples)
pn_sig = var(ns.samples)^0.5
pn_skew = sum( (ns.samples - pn_mu)^3 ) /((N-1) * (pn_sig)^3)
pn_kurt = sum( (ns.samples - pn_mu)^4 ) /((N-1) * (pn_sig)^4)
print(cat(pn_mu,pn_sig2, pn_skew, pn_kurt))
theo_mean = alpha*lambda + mu
theo_var = lambda*(alpha^2 + sig^2) + si^2
theo_skew = 0 + (lambda*(alpha^3) + 3*lambda*alpha*(sig^2))/(theo_var^1.5)
theo_kurt = (lambda*(alpha^4) + 6*lambda*(alpha^2)*(sig^2) + 3*lambda*(sig^4) )/(theo_var^2) + 3
print(cat(theo_mean,theo_var, theo_skew, theo_kurt))