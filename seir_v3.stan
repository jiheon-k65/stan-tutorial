functions 
{   
    // Function that switches eta value based on t
    real switch_eta (real t,
                    real t1,
                    real eta,
                    real nu,
                    real xi) 
    {
        return (eta + (1 - eta) / (1 + exp(xi * (t - t1 - nu))));
    }

    vector sir (real t,
                vector y,
                real beta,
                real delta,
                real gamma,
                real eta,
                real nu,
                real xi,
                real e0,
                real i0,
                real t_switch,
                int N) 
    {   
        // To account for decreased beta due to control measures
        real forcing_function = switch_eta(t, t_switch, eta, nu, xi);
        real beta_eff = beta * forcing_function; 
        
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
    vector[4] y0; // Initial conditions
    real t0; // Initial time point
    int<lower = 1> n_days; // No. of days
    array[n_days] real ts; // Array of times
    array[n_days] int cases; // Array of observed cases
    int N; // Population size
    real t_switch; // Start time of control measures
    int t_survey_start; // Start time of serological survey
    int t_survey_end; // End time of serological survey
    int n_survey_inf; // No. of people surveyed with past infection
    int n_survey_test; // Total no. of people surveyed
}

parameters 
{
    real<lower = 0> beta; // Effective contact rate
    real<lower = 0> delta; // Rate of being infectious
    real<lower = 0> gamma; // Rate of recovery
    real<lower = 0> phi_inv;
    real<lower = 0, upper = 1> eta; // Decrease in transmission factor
    real<lower = 0> nu; // Delay from introduction until measure is 50% effective
    real<lower = 0, upper = 1> xi_raw; 
    real<lower = 0> e0; // Initial no. of exposed people
    real<lower = 0> i0; // Initial no. of infected people
    real<lower = 0, upper = 1> p_reported; // Prob. of reporting
}

transformed parameters 
{
    real phi = 1. / phi_inv;
    real xi = xi_raw + 0.5;
    
    // To initialize an empty vector
    vector[4] y0_ip = rep_vector(0.0, 4);
    // To pass y0_ip (not y0) to the ODE solver
    array[n_days] vector[4] y = ode_rk45(sir, y0_ip, t0, ts, e0, i0, beta, delta, gamma, eta, nu, xi, t_switch, N);

    // To calculate prob. of past infection among those surveyed
    real<lower = 0, upper = 1> p_survey_inf;
    p_survey_inf = mean(to_vector(y[t_survey_start:t_survey_end, 4])) / N; // Mean of R components
    
    // To define an array for incidence
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
    eta ~ beta(2.5, 4);
    nu ~ exponential(1. / 5);
    xi_raw ~ beta(1, 1);
    p_reported ~ beta(1, 2);
    e0 ~ normal(0, 10);
    i0 ~ normal(0, 10);
    cases[1:(n_days - 1)] ~ neg_binomial_2(incidence, phi);
    n_survey_inf ~ binomial(n_survey_test, p_survey_inf);
}

generated quantities
{
    real R0 = beta / gamma;
    real incubation_time = 1 / delta;
    real recovery_time = 1 / gamma;

    // To calculate pred_cases
    array[n_days - 1] real pred_cases = neg_binomial_2_rng(incidence, phi);
    
    // To calculate effective reproduction no.
    array [n_days] real Reff;
    for (t in 1:n_days) {
        Reff[t] = switch_eta(t, t_switch, eta, nu, xi) * beta / gamma;
    }
}
