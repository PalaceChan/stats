# 1 var case
N <- 300
x1 <- data.table(x=rnorm(N), u=rnorm(N))
x2 <- data.table(x=rnorm(N), u=rnorm(N))
x1[, `:=`(y1=2*x + rnorm(N, sd=0.25), y2=3*x + rnorm(N, sd=1), y3=5*x + rnorm(N, sd=3))]
x2[, `:=`(y1=2*x + rnorm(N, sd=0.25), y2=3*x + rnorm(N, sd=1), y3=5*x + rnorm(N, sd=3))]
models <- list(x1[, lm(y1 ~ x - 1)], x1[, lm(y2 ~ x - 1)])
x2[, `:=`(m1.hat=predict(models[[1]], newdata = x2), m2.hat=predict(models[[2]], newdata = x2))]
cor(x2)

# 1 var and intercept case
N <- 300
x1 <- data.table(x=rnorm(N), u=rnorm(N))
x2 <- data.table(x=rnorm(N), u=rnorm(N))
x1[, `:=`(y1=2 + 2*x + rnorm(N, sd=0.25), y2=3 + 3*x + rnorm(N, sd=1), y3=5 + 5*x + rnorm(N, sd=3))]
x2[, `:=`(y1=2 + 2*x + rnorm(N, sd=0.25), y2=3 + 3*x + rnorm(N, sd=1), y3=5 + 5*x + rnorm(N, sd=3))]
models <- list(x1[, lm(y1 ~ x)], x1[, lm(y2 ~ x)])
x2[, `:=`(m1.hat=predict(models[[1]], newdata = x2), m2.hat=predict(models[[2]], newdata = x2))]
cor(x2)

# 2 slopes case
x1[, `:=`(y1=2*x + 7*u + rnorm(N, sd=0.25), y2=3*x + 4*u + rnorm(N, sd=1), y3=5*x + 1*u + rnorm(N, sd=3))]
x2[, `:=`(y1=2*x + 7*u + rnorm(N, sd=0.25), y2=3*x + 4*u + rnorm(N, sd=1), y3=5*x + 1*u + rnorm(N, sd=3))]
models <- list(x1[, lm(y1 ~ x + u - 1)], x1[, lm(y2 ~ x + u - 1)])
x2[, `:=`(m1.hat=predict(models[[1]], newdata = x2), m2.hat=predict(models[[2]], newdata = x2))]
cor(x2)

# 2 var orthogonal to each other
N <- 4000
M1 <- pracma::randortho(N)
M2 <- pracma::randortho(N)
x1 <- data.table(u1=M1[,1], u2=M1[,2])
x2 <- data.table(u1=M2[,1], u2=M2[,2])
x1[, `:=`(y1=2*u1 + 7*u2 + rnorm(N, sd=0.05), y2=3*u1 + 4*u2 + rnorm(N, sd=1), y3=5*u1 + 1*u2 + rnorm(N, sd=3))]
x2[, `:=`(y1=2*u1 + 7*u2 + rnorm(N, sd=0.05), y2=3*u1 + 4*u2 + rnorm(N, sd=1), y3=5*u1 + 1*u2 + rnorm(N, sd=3))]
models <- list(x1[, lm(y1 ~ u1 + u2 - 1)], x1[, lm(y2 ~ u1 + u2 - 1)])
x1[, `:=`(m1.hat=predict(models[[1]]), m2.hat=predict(models[[2]]))]
x2[, `:=`(m1.hat=predict(models[[1]], newdata = x2), m2.hat=predict(models[[2]], newdata = x2))]
cor(x2)
