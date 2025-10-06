functions
{   
    vector sir (real t,
                vector y,
                real beta,
                real gamma,
                int N)
    // No need to pack parameters into theta
    // No need to pack data into x_r & x_i
    {   
        vector[3] dydt;
        dydt[1] = - beta * (y[2] / N) * y[1];
        dydt[2] = beta * (y[2] / N) * y[1] - gamma * y[2];
        dydt[3] = gamma * y[2];
        // Where y[1] = S, y[2] = I & y[3] = R
        return dydt;
    }
}

data {
    vector[3] y0; // Initial conditions
    real t0; // Initial time point
    int<lower = 1> n_days; // No. of days
    array[n_days] real ts; // Array of time points
    array[n_days] int cases; // Array of cases
    int N; // Population size
}

parameters 
{
    real<lower = 0> beta; // Effective contact rate
    real<lower = 0> gamma; // Rate of recovery
    real<lower = 0> phi_inv;
}

transformed parameters 
{
    real phi = 1. / phi_inv; // Dispersion factor
    // To create an array over n_days, each entry a vector of y
    array[n_days] vector[3] y = ode_rk45(sir, y0, t0, ts, beta, gamma, N);
}

model
{
    // Priors
    beta ~ normal(2, 1);
    gamma ~ normal(0.4, 0.5);
    phi_inv ~ exponential(5);
    
    // Likelihood or sampling distribution
    cases ~ neg_binomial_2(y[:, 2], phi);
    // [:, 2] selects for all rows of column 2 i.e. I_{ODE}(t)
    // This is used to "generate" the overdispersed observed cases I_{obs}(t) 
}

generated quantities
{
    real R0 = beta / gamma;
    real recovery_time = 1 / gamma;
    array[n_days] real pred_cases;
    for (t in 1:n_days) {
        pred_cases[t] = neg_binomial_2_rng(fmax(y[t, 2], 1e-5), phi);
    // fmax(..., 1e-5) ensures positive counts
    }
}
