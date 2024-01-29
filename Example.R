source("Isotonic-function.R")

## setting
piececons_func = function(x){
  if(x <= 25){0}
  else if((25 < x) && (x <= 80)){2.5}
  else {3}
}

n = p = 100
theta_pcons = sapply(1:100, piececons_func)

## data generation
set.seed(1)
X = diag(n)
sigma = 0.25
y_pcons = rnorm(n, theta_pcons, sigma)

## estimation
result_pcons_HS = BIR(y_pcons, X, r=1/2,a_lam = 1/2, mc=2500, burn=500)
etahat_pcons_HS = apply(result_pcons_HS$eta, 2, mean)
D = lower.tri(diag(p), diag = T)
thetahat_pcons_HS = D %*% etahat_pcons_HS
low_bound_pcons_HS = c()
upper_bound_pcons_HS = c()
for(i in 1:100){
  low_bound_pcons_HS[i] = quantile((result_pcons_HS$eta %*% t(D))[, i], probs = 0.025)
  upper_bound_pcons_HS[i] = quantile((result_pcons_HS$eta %*% t(D))[, i], probs = 0.975)
}

## plot
plot(1:p, y_pcons, pch=1, lwd=2, cex.lab = 1.5, col="gray", 
     xlab="x", ylab="y",xaxt="n", ylim = c(0, 4), cex.axis =1.5, 
     type="n", main="Piecewise Constant", yaxt="n", cex.main =2)
axis(side=1, at=c(0, 20, 40, 60, 80, 100), labels=T)
axis(side=2, at=0:4, labels=T)
xx = c(1, 1:100, 100:1)
yy = c(low_bound_pcons_HS[1], upper_bound_pcons_HS, rev(low_bound_pcons_HS))
polygon(xx, yy, col="gray75", border = NA)
points(y_pcons, pch=1, lwd=1, col="gray50", cex=0.5)
points(thetahat_pcons_HS, col="grey20", type="l", lwd=2)
points(theta_pcons, col="red", type="l", lty=2, lwd=2)

