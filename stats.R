library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
options('mc.cores' = 4)
## library(shinystan)
## options('browser' = 'firefox')
library(lme4)
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
