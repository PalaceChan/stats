data {
    int N;
    int J;

    real y[N];
    real x[N];
    int county[N];
}

parameters {
    real a[J];
    real b;

    real mu_a;
    real<lower=0> sigma_a;
    real<lower=0> sigma_y;
}

model {
    for (i in 1:N) {
        y[i] ~ normal(a[county[i]] + b * x[i], sigma_y);
    }

    for (j in 1:J) {
        a[j] ~ normal(mu_a, sigma_a);
    }

    mu_a ~ normal(0, 100);
    sigma_a ~ uniform(0, 100);
}
