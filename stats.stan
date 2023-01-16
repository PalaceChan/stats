data {
    int N;
    int J;

    real y[N];
    real x[N];
    real u[J];
    int county[N];
}

parameters {
    real a[J];
    real b;

    real g0;
    real g1;
    real<lower=0> sigma_a;
    real<lower=0> sigma_y;
}

model {
    for (i in 1:N) {
        y[i] ~ normal(a[county[i]] + b * x[i], sigma_y);
    }

    for (j in 1:J) {
        a[j] ~ normal(g0 + g1 * u[j], sigma_a);
    }
}
