functions 
{   
    vector sir (real t,
                vector y,
                real beta,
                real gamma,
                int N) 
    {   
        vector[3] dydt;
        dydt[1] = - beta * (y[2] / N) * y[1];
        dydt[2] = beta * (y[2] / N) * y[1] - gamma * y[2];
        dydt[3] = gamma * y[2];
        return dydt;
    }
}

data {
    vector[3] y0;
    real t0;
    int<lower = 1> n_days;
    array[n_days] real ts;
    array[n_days] int cases;
    int N;
    int compute_likelihood;
}

parameters 
{
    real<lower = 0> beta;
    real<lower = 0> gamma;
    real<lower = 0> phi_inv;
}

transformed parameters 
{
    real phi = 1. / phi_inv;
    array[n_days] vector[3] y = ode_rk45(sir, y0, t0, ts, beta, gamma, N);
}

model
{
    // Priors
    beta ~ normal(2, 1);
    gamma ~ normal(0.4, 0.5);
    phi_inv ~ exponential(5);
    
    // To toggle likelihood on/off
    // Without sampling distribution, parameters are not fitted to data
    if (compute_likelihood == 1)
        cases ~ neg_binomial_2(y[:, 2], phi);
}

generated quantities
{
    real R0 = beta / gamma;
    real recovery_time = 1 / gamma;
    array[n_days] real pred_cases;
    for (t in 1:n_days) {
        pred_cases[t] = neg_binomial_2_rng(fmax(y[t, 2], 1e-5), phi);
    }
}
