functions 
{   
    vector sir (real t,
                vector y,
                real beta,
                real delta,
                real gamma,
                real e0,
                real i0,
                int N) 
    {   
        // To define "true" initial conditions within ODE function
        vector[4] init;
        init[1] = N - e0 - i0;
        init[2] = e0;
        init[3] = i0;
        init[4] = 0;

        // To shift initial y[t] by init[t]
        vector[4] y_ip = y + init;
         
        // To pass y_ip[t] (not y[t]) to derivatives
        vector[4] dydt;
        dydt[1] = - beta * (y_ip[2] / N) * y_ip[1];
        dydt[2] = beta * (y_ip[2] / N) * y_ip[1] - delta * y_ip[2];
        dydt[3] = delta * y_ip[2] - gamma * y_ip[3];
        dydt[4] = gamma * y_ip[3];
        return dydt;
    }
}

data {
    vector[4] y0;
    real t0;
    int<lower = 1> n_days;
    array[n_days] real ts;
    array[n_days] int cases;
    int N;
}

parameters 
{
    real<lower = 0> beta;
    real<lower = 0> delta;
    real<lower = 0> gamma;
    real<lower = 0> phi_inv;
    real<lower = 0, upper = 1> p_reported;
    real<lower = 0> e0; // Initial no. of exposed people
    real<lower = 0> i0; // Initial no. of infected people
}

transformed parameters 
{
    real phi = 1. / phi_inv;

    // To initialize an empty vector
    vector[4] y0_ip = rep_vector(0.0, 4);
    // To pass y0_ip (not y0) to the ODE solver
    array[n_days] vector[4] y = ode_rk45(sir, y0_ip, t0, ts, e0, i0, beta, delta, gamma, N);

    // To initialize an empty array for incidence
    array[n_days - 1] real incidence;
    // To calculate incidence
    for (t in 1:(n_days - 1)) {
        incidence[t] = (y[t, 1] + y[t, 2] - y[(t + 1), 1] - y[(t + 1), 2]) * p_reported;
    }
}

model
{
    beta ~ normal(2, 1);
    delta ~ normal(0.4, 0.5);
    gamma ~ normal(0.4, 0.5);
    phi_inv ~ exponential(5);
    cases[1:(n_days - 1)] ~ neg_binomial_2(incidence, phi);
    p_reported ~ beta(1, 2);
    e0 ~ normal(0, 10);
    i0 ~ normal(0, 10);
}

generated quantities
{
    real R0 = beta / gamma;
    real incubation_time = 1/delta;
    real recovery_time = 1 / gamma;
    array[n_days - 1] real pred_cases;
    pred_cases = neg_binomial_2_rng(incidence, phi);
}
