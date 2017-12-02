
powa <- function(n1=300, n2=300, p1=0.45, p2=0.55, cs=12.5, icc=0.02, cv=0.56, alpha=0.05) {
  deff <- 1 + ((cv^2 + 1) * cs - 1) * icc
  z <- qnorm(1 - alpha/2)
  q1 <- 1 - p1
  q2 <- 1 - p2
  pm <- (n1 * p1 + n2 * p2)/(n1 + n2)
  ds <- z * sqrt(deff*(1/n1 + 1/n2) * pm * (1 - pm))
  ex <- abs(p1 - p2)
  sd <- sqrt(deff*(p1 * q1/n1 + p2 * q2/n2))
  1 - pnorm((ds - ex)/sd) + pnorm((-ds - ex)/sd)
}

# run our power stuff. We want to generate for n=200,250,300,350 across
# a range of P and output the power
pow_range <- expand.grid(n=seq(200,400,by=50), p1=c(0.05,0.1,0.2,0.3,0.5), delta_p=seq(0.05,0.3,by=0.01))
pow_range$power=powa(n2=pow_range$n, p1=pow_range$p1, p2=pow_range$p1 + pow_range$delta_p)

ggplot(pow_range) + geom_line(aes(x=delta_p,y=power,col=factor(n))) + facet_wrap(~p1)

# now a computation to figure out minimal delta_p detected given power=0.8
powa_dp <- function(delta_p, n1=300, n2=300, p1=0.45, power=0.8, cs=12.5, icc=0.02, cv=0.56, alpha=0.05) {
  powa(n1, n2, p1, p1+delta_p, cs, icc, cv, alpha) - power
}

delta_pw <- function(n1=300, n2=300, p1=0.45, power=0.8, cs=12.5, icc=0.02, cv=0.56, alpha=0.05) {
  uniroot(powa_dp, interval=c(1e-3, 0.5), n1=n1, n2=n2, p1=p1, power=power, cs=cs, icc=icc, cv=cv, alpha=alpha)$root
}

pow_range <- expand.grid(n=seq(200,400,by=50), p1=seq(0.01,0.4,by=0.01), power=c(0.6,0.7,0.8,0.9))
delta_p <- numeric(nrow(pow_range))
for (i in 1:nrow(pow_range)) {
  pr <- pow_range[i,]
  delta_p[i] <- delta_pw(n1=pr$n, p1=pr$p1, power=pr$power, icc=0.28)
}
pow_range$delta_p = delta_p

ggplot(pow_range) + geom_line(aes(x=p1,y=delta_p,col=factor(power))) + facet_wrap(~n)


if (pooled) {
  p <- (p1 + p2)/2
  sdd <- sqrt(p * (1 - p) * 2 * DEFF/(m * n))
} else {
  sdd <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) * DEFF/(m * n))
}

zcrit <- qnorm(alpha/2, lower.tail = FALSE)

pnorm(abs(p1 - p2)/sdd - zcrit, lower.tail = TRUE)
