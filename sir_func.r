# Function to compile & fit Stan model
fit_model <- function(
    stan_file,
    data,
    seed = 0,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    max_treedepth = NULL,
    adapt_delta = NULL,
    show_message = FALSE) {
    stan_path <- paste0(stan_file, ".stan")
    stan_mod <- cmdstan_model(stan_path)
    fit_mod <- stan_mod$sample(
        data = data,
        seed = seed,
        chains = chains,
        parallel_chains = parallel_chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        max_treedepth = max_treedepth,
        adapt_delta = adapt_delta,
        show_message = show_message
    )
    return(fit_mod)
}

# Function to remove leading dots from colnames
remove_dots <- function(dt) {
    dt_adj <- setnames(
        dt,
        names(dt),
        str_replace(names(dt), "^\\.", "")
    )
    return(dt_adj)
}

# Function to calculate time point
calc_time <- function(df, event_date) {
    time <- nrow(df[df$date < event_date, ]) + 1
    return(time)
}

# Function to summarize pred_cases
summ_pred <- function(dt) {
    cases_trim <- cases[-length(cases)] # To remove last index
    dt_summary <- dt[, .(
        mean = mean(pred_cases),
        lower = quantile(pred_cases, probs = 0.025),
        upper = quantile(pred_cases, probs = 0.975)
    ),
    by = .(t)
    ]
    dt_summary[, cases := cases_trim] # To add cases
    return(dt_summary)
}

# Function to summarize pred_cases by chain
summ_pred_4chains <- function(dt) {
    cases_trim <- cases[-length(cases)] # To remove last index
    dt_summary <- dt[, .(
        mean = mean(pred_cases),
        lower = quantile(pred_cases, probs = 0.025),
        upper = quantile(pred_cases, probs = 0.975)
    ),
    by = .(chain, t)
    ]
    dt_summary[, cases := cases_trim, by = chain] # To add cases by chain
    return(dt_summary)
}

# Function to plot pred_cases
plot_pred <- function(dt) {
    plot <- ggplot(data = dt, aes(x = t)) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% CrI"), alpha = 0.2) +
        geom_line(aes(y = mean, color = "Mean"), linewidth = 1) +
        geom_point(aes(y = cases, color = "Data"), size = 1) +
        labs(x = "Days", y = "No. of newly reported cases", fill = NULL, color = NULL) +
        scale_fill_manual(values = c("95% CrI" = "darkorange")) +
        scale_color_manual(values = c("Mean" = "darkorange", "Data" = "black")) +
        theme_bw(base_size = 14) +
        theme(legend.position = "bottom")
    return(plot)
}

# Function to plot pred_cases by chain
plot_pred_4chains <- function(dt) {
    plot <- ggplot(data = dt, aes(x = t)) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% CrI"), alpha = 0.2) +
        geom_line(aes(y = mean, color = "Mean"), linewidth = 1) +
        geom_point(aes(y = cases, color = "Data"), size = 1) +
        labs(x = "Days", y = "No. of newly reported cases", fill = NULL, color = NULL) +
        facet_wrap(~chain) +
        scale_fill_manual(values = c("95% CrI" = "darkorange")) +
        scale_color_manual(values = c("Mean" = "darkorange", "Data" = "black")) +
        theme_bw(base_size = 14) +
        theme(legend.position = "bottom")
    return(plot)
}

# Function to plot trace plots for model parameters
plot_trace <- function(fitted_model, model_parm) {
    samples <- fitted_model$draws(variables = model_parm)
    color_scheme_set("brewer-Spectral")
    plot <- mcmc_trace(samples) + theme_bw(base_size = 14) + theme(legend.position = "bottom")
    return(plot)
}

# Function to plot pair plots for model parameters
plot_pairs <- function(fitted_model, model_parm) {
    # To draw posterior samples from fitted model
    samples <- fitted_model$draws(variables = model_parm, format = "matrix")
    # To create dataframe for divergent points
    diverg_mat <- fitted_model$sampler_diagnostics()[, , "divergent__"]
    n_iter <- nrow(diverg_mat)
    n_chain <- ncol(diverg_mat)
    diverg_df <- data.frame(
        Chain = rep(1:n_chain, each = n_iter),
        Iteration = rep(1:n_iter, times = n_chain),
        Parameter = "divergent__",
        Value = c(diverg_mat)
    )
    # To plot pair plot
    color_scheme_set("brightblue")
    plot <- mcmc_pairs(
        samples,
        diag_fun = "dens",
        diag_args = list(alpha = 0.6),
        off_diag_fun = "scatter",
        off_diag_args = list(size = 0.6, alpha = 0.6),
        np = if (any(diverg_df$Value == 1)) diverg_df else NULL,
        np_style = pairs_style_np(div_color = "brown1", div_shape = 3, div_size = 2),
        max_treedepth = 10
    )
    return(plot)
}

# Function to summarize Reff
summ_Reff <- function(dt) {
    dt_summary <- dt[, .(
        mean = mean(Reff),
        lower = quantile(Reff, 0.025),
        upper = quantile(Reff, 0.975)
    ),
    by = .(t)
    ]
    return(dt_summary)
}

# Function to plot Reff
plot_Reff <- function(dt) {
    plot <- ggplot(data = dt, aes(x = t)) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% CrI"), alpha = 0.2) +
        geom_line(aes(y = mean, color = "Mean"), linewidth = 1) +
        geom_vline(aes(xintercept = t_switch), linetype = 2, linewidth = 1) +
        labs(x = "Days", y = "Effective reproduction number", fill = NULL, color = NULL) +
        scale_fill_manual(values = c("95% CrI" = "darkorange")) +
        scale_color_manual(values = c("Mean" = "darkorange", "Data" = "black")) +
        theme_bw(base_size = 14) +
        theme(legend.position = "bottom")
    return(plot)
}

# Function to plot posterior distributions by chain
plot_post <- function(dt) {
    plot <- ggplot(data = dt) +
        geom_density(aes(x = value, color = factor(chain), fill = factor(chain)), alpha = 0.4) +
        labs(x = "Value", y = "Density", color = "Chain", fill = "Chain") +
        facet_wrap(~variable, scales = "free") +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        theme_bw(base_size = 14) +
        theme(legend.position = "bottom")
    return(plot)
}
