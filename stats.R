library(data.table)
library(lme4)
library(rstan)
library(purrr)
rstan_options(auto_write = TRUE)
options('mc.cores' = 4)

fitMods <- function(N, J, sy, sa, b) {
    set.seed(1)
    ma <- 0
    a <- rnorm(J, ma, sa)
    x <- rnorm(N)
    j <- c(1:J, sample(1:J, N-J, replace = T, runif(J))) # group membership ensure at least one of each
    y <- a[j] + b*x + rnorm(N, sd = sy)

    ## estimate each 'a' as a single pooling
    m0 <- lm(y ~ x)

    ## estimate each 'a' by itself
    m1 <- lm(y ~ factor(j) + x - 1)

    ## estimate two level model with lme4
    m2 <- lmer(y ~ x + (1 | j))

    ## estimate two level model with stan
    spec <- "
data {
    int N;
    int J;
    real y[N];
    real x[N];
    int grp[N];
}
parameters {
    real a[J];
    real b;
    real ma;
    real<lower=0> sa;
    real<lower=0> sy;
}
model {
    for (i in 1:N) {
        y[i] ~ normal(a[grp[i]] + b * x[i], sy);
    }

    for (j in 1:J) {
        a[j] ~ normal(ma, sa);
    }
}
"
    ## m3 <- stan(model_code = spec, data = list(N=N, J=J, y=y, x=x, grp=j), refresh = 0)

    res <- data.table(prm = c(paste0('a', 1:J), 'b', 'sy', 'sa', 'ma'),
                      obs = c(table(j), rep(N, 4)),
                      tru = c(a, b, sy, sa, ma),
                      poo = c(rep(coef(m0)[1], J), coef(m0)[1], summary(m1)$sigma, NA, NA),
                      unp = c(coef(m1), summary(m1)$sigma, NA, NA),
                      lme = c(coef(m2)$j[,1], fixef(m2)[2], summary(m2)$sigma, summary(m2)$varcor$j, fixef(m2)[1]))
                      ## stn = get_posterior_mean(m3)[, 'mean-all chains'][c(1:(J+1), J+4, J+3, J+2)])
    ## list(res=res, perf=cor(res[1:J, .(tru, unp, lme, stn)])[1,])
    ## list(res=res, perf=cor(res[1:J, .(tru, unp, lme)])[1,])
    tot.err <- sum(res$tru[1:J]^2)
    list(res=res,
         perf=c(poo=sum((res$tru[1:J] - res$poo[1:J])^2), unp=sum((res$tru[1:J] - res$unp[1:J])^2), lme=sum((res$tru[1:J] - res$lme[1:J])^2))/tot.err)
}

xl <- lapply(seq(0.01, 2, by=.05), function(v) fitMods(1000, 10, 1, v, 1.7)) # lme does better with the smaller values
xl <- lapply(c(2, 5, 10, 30, 100, 200, 500), function(v) fitMods(1000, v, 2, 1, 1.7)) # lme does better with more groups
do.call(rbind, map(xl, 'perf'))

N <- 1e4
J <- 100
B <- 50
W.x <- matrix(rnorm(10*B), ncol=B)
S.x <- t(W.x) %*% W.x + diag(runif(B, 0.1, 1))
X <- MASS::mvrnorm(N, mu = rep(0, B), Sigma = S.x, empirical = T)
j <- sort(c(1:J, sample(1:J, N-J, replace = T, runif(J)))) # group membership ensure at least one of each
mu.b <- rnorm(B)
W.b <- matrix(rnorm(10*B), ncol=B)
S.b <- t(W.b) %*% W.b + diag(runif(B, 0.1, 1))
C.b <- diag(1/sqrt(diag(S.b))) %*% S.b %*% diag(1/sqrt(diag(S.b)))
b <- MASS::mvrnorm(max(J,B), mu = mu.b, Sigma = C.b, empirical = T) # j-th row are the true betas for the j-th group
y <- do.call(rbind, lapply(1:J, function(grp) X[grp == j,] %*% b[grp,]))
jf <- factor(j, ordered = F)
dat <- cbind(y, jf, as.data.frame(X))
## h(model.matrix(y ~ X:jf - 1, data = dat))
m1 <- lm(y ~ X:jf - 1, data = dat)
m2 <- lmer(y ~ 0+X + (0+X | jf))
m1.cc <- matrix(nrow = J, ncol = B)
for (i in 1:J) {
    for (k in 1:B) {
        m1.cc[i,k] <- coef(m1)[paste0('X',k,':jf',i)]
    }
}

b[1:J,]
m1.cc
coef(m2)$jf
summary(as.numeric(m1.cc - b[1:J,]))
summary(as.numeric(as.matrix(coef(m2)$jf) - b[1:J,]))

srrs2 <- fread("~/Downloads/srrs2.dat")
mn <- srrs2[state == 'MN']
mn[, y := log(ifelse(activity == 0, .1, activity))]
mn[, x := floor]
mn[, cnty := match(county, unique(county))]
mn[, c("county", "cnty") := .(cnty, county)]

cty <- fread("~/Downloads/cty.dat")
srrs2.fips <- mn[, stfips*1000 + cntyfips]
usa.fips <- cty[, stfips*1000 + ctfips]
usa.rows <- match(unique(srrs2.fips), usa.fips)
uranium <- cty[usa.rows,"Uppm"]
mn[, u := log(uranium[county])]

lm.pooled <- mn[, lm(y ~ x + u)]
## lm.unpooled <- mn[, lm(y ~ x + factor(county) - 1)]
er.fit <- mn[, lmer(y ~ x + u + (1 | county))]
ml.fit <- stan(file='~/development/stats/stats.stan', data = list(N=nrow(mn), J=mn[, uniqueN(county)], y=mn$y, x=mn$x, u=mn$u, county=mn$county))
x <- data.table(name = rownames(get_posterior_mean(ml.fit))[1:86],
                nobs = c(mn[, .N, by = county]$N, nrow(mn)),
                mu.poo = c(coef(lm.pooled)[1] + coef(lm.pooled)[3]*log(uranium$Uppm), coef(lm.pooled)[2]),
                mu.lme = c(fixef(er.fit)[1] + fixef(er.fit)[3]*log(uranium$Uppm) + ranef(er.fit)$county[,1], fixef(er.fit)[2]),
                mu.mlm = get_posterior_mean(ml.fit)[1:86, 5])
## mu.unp = coef(lm.unpooled)[c(2:86, 1)],
## se.unp = summary(lm.unpooled)$coefficients[c(2:86, 1), 'Std. Error'],
## se.poo = c(rep(summary(lm.pooled)$coefficients[, 'Std. Error'][1], 85), summary(lm.pooled)$coefficients[, 'Std. Error'][2]),
## se.lme = c(sqrt(drop(attr(ranef(er.fit)$county, "postVar"))), sqrt(diag(vcov(er.fit)))[2]),
## se.mlm = summary(ml.fit)$summary[1:86, 'sd']


## fake data checking
N <- nrow(mn)
J <- mn[, uniqueN(county)]
b.true <- -0.7
g0.true <- 1.5
g1.true <- 0.7
sigma.y.true <- 0.9
sigma.a.true <- 0.09
a.true <- rnorm(J, g0.true + g1.true*log(uranium$Uppm), sigma.a.true)
y.fake <- rnorm(N, a.true[mn$county] + b.true*mn$x, sigma.y.true)
lm.pooled <- mn[, lm(y.fake ~ x + u)]
er.fit <- mn[, lmer(y.fake ~ x + u + (1 | county))]
ml.fit <- stan(file='~/development/stats/stats.stan', data = list(N=N, J=J, y=y.fake, x=mn$x, u=log(uranium$Uppm), county=mn$county))
x <- data.table(name = rownames(get_posterior_mean(ml.fit))[1:86],
                nobs = c(mn[, .N, by = county]$N, nrow(mn)),
                mu.tru = c(a.true, b.true),
                mu.poo = c(coef(lm.pooled)[1] + coef(lm.pooled)[3]*log(uranium$Uppm), coef(lm.pooled)[2]),
                mu.lme = c(fixef(er.fit)[1] + fixef(er.fit)[3]*log(uranium$Uppm) + ranef(er.fit)$county[,1], fixef(er.fit)[2]),
                mu.mlm = get_posterior_mean(ml.fit)[1:86, 5])
x[1:85, .(mu.tru, mu.poo, mu.lme, mu.lme)] |> cor()
