y <- as.matrix(read.table('https://raw.github.com/wiki/stan-dev/rstan/rats.txt', header = TRUE))
x <- c(8, 15, 22, 29, 36)
xbar <- mean(x)
N <- nrow(y)
T <- ncol(y)
rats_fit <- stan(file='~/development/stats/stats.stan', data = list(N=N, T=T, y=y, x=x, xbar=xbar))

library(rstan)
rstan_options(auto_write = TRUE)
options('mc.cores' = 4)

library(shinystan)
options('browser' = 'firefox')
launch_shinystan(rats_fit)

library(lme4)
srrs2 <- fread("~/Downloads/srrs2.dat")
mn <- srrs2[state == 'MN']
mn[, y := log(ifelse(activity == 0, .1, activity))]
mn[, x := floor]
mn[, cnty := match(county, unique(county))]
mn[, c("county", "cnty") := .(cnty, county)]

lm.pooled <- mn[, lm(y ~ x)]
lm.unpooled <- mn[, lm(y ~ x + factor(county) - 1)]
er.fit <- mn[, lmer(y ~ x + (1 | county))]
ml.fit <- stan(file='~/development/stats/stats.stan', data = list(N=nrow(mn), J=mn[, uniqueN(county)], y=mn$y, x=mn$x, county=mn$county))
x <- data.table(name = names(m[1:86, 5]),
                nobs = c(mn[, .N, by = county]$N, nrow(mn)),
                mu.poo = c(rep(coef(lm.pooled)[1], 85), coef(lm.pooled)[2]),
                mu.unp = coef(lm.unpooled)[c(2:86, 1)],
                mu.lme = c(coef(er.fit)$county[, 1], fixef(er.fit)[2]),
                mu.mlm = get_posterior_mean(ml.fit)[1:86, 5],
                se.poo = c(rep(summary(lm.pooled)$coefficients[, 'Std. Error'][1], 85), summary(lm.pooled)$coefficients[, 'Std. Error'][2]),
                se.unp = summary(lm.unpooled)$coefficients[c(2:86, 1), 'Std. Error'],
                se.lme = c(sqrt(drop(attr(ranef(er.fit)$county, "postVar"))), sqrt(diag(vcov(er.fit)))[2]),
                se.mlm = summary(ml.fit)$summary[1:86, 'sd'])
