#
# Utility functions for fitting and summarising logistic growth models
# 


# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=
# ---- DATA PROCESSING ---------------------------------------------------------------------------
# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=

#' Remove downstream observations once growth rate starts decreasing.
#' 
#' @description 
#' Each time series will be cut at the first point which is 
#' $\geq$ `percent_cutoff` percent lower than the previous observation, 
#' with this point and all subsequent observations discarded.
#'  
#' Setting `percent_cutoff = 0` removes *all* observations which are less 
#' than the previous point in the same timeseries. Setting 
#' `percent_cutoff = Inf` removes no points.
#' 
#' All corresponding points from other replicates/viruses with the same 
#' treatment will also be removed. This ensures downstream fits for a 
#' given treatment always use the same number of timepoints for all
#' replicates.
#'
#' @param data A data frame which includes a `time` column.
#' @param percent_cutoff Cutoff (in percent).
#' @param timeseries_identifier Name of a column uniquely identifying 
#' individual timeseries.
#' @param counts_column Name of the column containing counts.
#' @param treatment_column Name of the column specifing unique treatments.
#' 
#' @return A data frame.
remove_decreases <- function(data, percent_cutoff, timeseries_identifier = "timeseries_id", 
                             counts_column = "infected_count", treatment_column = "treatment") {
  data %>% 
    # Remove decreases
    group_by(.data[[timeseries_identifier]]) %>% 
    arrange(.data$time) %>% 
    mutate(change_percent = .data[[counts_column]] - lag(.data[[counts_column]]),
           change_percent = .data$change_percent / lag(.data[[counts_column]]) * 100,
           early_fluctuations = .data$time <= 5 & lag(.data[[counts_column]]) < 10,  # At very low counts the change_percent method is too sensitive
           keep = is.na(.data$change_percent) | .data$change_percent >= -percent_cutoff | .data$early_fluctuations) %>% 
    mutate(index = 1:n(),
           keep = cumsum(.data$keep) == .data$index) %>%   # Only true if this and all previous observations were kept
    filter(.data$keep) %>% 
    
    # Ensure all replicates/viruses at a given concentration use the same number of points:
    mutate(Npoints = n()) %>% 
    group_by(.data[[treatment_column]]) %>% 
    filter(.data$index <= min(.data$Npoints)) %>% 
    ungroup() %>% 
    
    # Remove columns added here
    select(-one_of("change_percent", "keep", "index", "Npoints"))
}


#' Prepare data for fitting
#' 
#' @description 
#' The first observation in each timeseries is used as the starting population, `P0`.
#' Also removes downstream observations representing negative growth (see 
#' `remove_decreases()`.
#'
#' @param data A data frame which includes a `time` column. Time should be measured in days.
#' @param percent_cutoff Cutoff for negative growth (passed to `remove_decreases()`)
#' @param timeseries_identifier Name of a column uniquely identifying 
#' individual timeseries.
#' @param counts_column Name of the column containing counts.
#' @param treatment_column Name of the column specifing unique treatments.
#' 
#' @return A data frame.
prepare_data <- function(data, percent_cutoff, timeseries_identifier = "timeseries_id", 
                         counts_column = "infected_count", treatment_column = "treatment") {
  processed_data <- data %>% 
    remove_decreases(percent_cutoff = percent_cutoff, 
                     timeseries_identifier = timeseries_identifier, 
                     counts_column = counts_column, 
                     treatment_column = treatment_column)
  
  starting_pop <- processed_data %>% 
    filter(.data$time == 1) %>% 
    mutate(P0 = .data[[counts_column]]) %>% 
    select(all_of(c(timeseries_identifier, "P0")))
  
  stopifnot(length(unique(starting_pop[[timeseries_identifier]])) == nrow(starting_pop)) # One per timeseries
  
  processed_data %>% 
    mutate(time = .data$time - 1) %>% 
    filter(.data$time != 0) %>% 
    left_join(starting_pop, by = timeseries_identifier) %>% 
    ungroup()
}


# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=
# ---- MODEL FITTING ------------------------------------------------------------------------------
# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=

#' Fit a logistic growth curve model
#'
#' @param model_data A prepared data frame (see `prepare_data()`). Expected 
#' columns are: `time`, `treatment`, `virus`, `replicate` and `infected_count`
#' @param include_interaction Should an treatment x virus interaction be
#' included?
#' 
#' @details 
#' `model_data` should contain a column named "infected_count", containing the
#' reponse variable (population sizes) over time. Replicates are assumed to 
#' uniquely identify *experimental* replicates.
#' 
#' @return A list containing the model data and the fitted model
fit_logistic <- function(model_data, include_interaction = TRUE, starting_k = 3000) {
  if (include_interaction) {
    fit <- nlme(infected_count ~ (K*P0*exp(r*time)) / (K + P0*(exp(r*time) - 1)),
                fixed = list(r ~ treatment + virus + treatment:virus,
                             K ~ treatment),
                random = list(timeseries_id = K ~ 1),
                data = model_data,
                start = c(r = c(2.5, 0, 0, 0),  # Intercept + 3 effects 
                          K = c(starting_k, 0)))      # Intercept + 1 effect
  } else {
    fit <- nlme(infected_count ~ (K*P0*exp(r*time)) / (K + P0*(exp(r*time) - 1)),
                fixed = list(r ~ treatment + virus,
                             K ~ treatment),
                random = list(timeseries_id = K ~ 1),
                data = model_data,
                start = c(r = c(2.5, 0, 0),
                          K = c(starting_k, 0)))
  }
  
  list(data = model_data,
       fit = fit)
}



# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=
# ---- BOOTSTRAPPING ------------------------------------------------------------------------------
# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=

#' Generate a single hierarchical bootstrap resample of a dataset, based on a fit from `fit_logistic()`
#' 
#' @description 
#' This function is needed to allow calculation of confidence intervals, since parametric bootstrap 
#' methods are not available for non-linear (nlme) models.
#' Based on https://stats.stackexchange.com/a/231320 & https://stats.stackexchange.com/a/232949.
#'
#' @param data The data frame used to fit the model.
#' @param fit The fitted model.
#' 
#' @return A resampled data frame.
resample_data <- function(data, fit) {
  # Sample groups (random effects) (with replacement)
  # - Note: currently this allows just one group
  resampled_replicate <- sample(unique(data$replicate), size = length(unique(data$replicate)), replace = TRUE)
  
  resampled_groups <- data.frame(replicate = resampled_replicate) # use expand.grid if adding more groups here
  
  
  # Re-sample residuals (within all groups, i.e. the lowest level)
  internal_data <- data.frame(data, 
                              inner_predictions = predict(fit),  # Predict at the inner-most level
                              residuals = resid(fit))
  
  resampled_data <- data.frame()
  
  for (g in 1:nrow(resampled_groups)) {
    current_replicate <- resampled_groups$replicate[g]
    
    group_data <- internal_data %>% 
      filter(.data$replicate == current_replicate)
    
    group_data$residuals <- sample(group_data$residuals, size = nrow(group_data), replace = TRUE)
    group_data$infected_count <- group_data$inner_predictions + group_data$residuals
    
    resampled_data <- rbind(resampled_data, group_data)
  }
  
  resampled_data %>% 
    select(-one_of("inner_predictions", "residuals"))
}


#' Generate `n` hierarchical bootstrap fits.
#' 
#' @description 
#' This function both resamples the data (see `resample_data()`) 
#' and re-fits the model on these data using `fit_logistic()`. Note 
#' that not all resampled datasets may be fittable (with the default
#' initial params) - this function will retry only once per iteration.
#' As a result, the number of returned results may be lower than `n`.
#' 
#' @param fit_obj A result from calling `fit_logistic()` on the original data.
#' @param n Number of repeats required
#' 
#' @return A list of lists containing matching pairs of resampled data and 
#' model fits.
get_boot_fits <- function(fit_obj, n = 100) {
  res <- replicate(n = n, simplify = FALSE, expr = {
    try({
      boot_data <- resample_data(fit_obj$data, fit_obj$fit)
      fit_logistic(boot_data) # returns both data and fit
    })
  })
  
  # Check for failures
  failed <- sapply(res, function(x) class(x) == "try-error")
  
  if (!any(failed)) {
    return(res)
  }
  
  # Try once to get all the missing replicates:
  res2 <- replicate(n = sum(failed), simplify = FALSE, expr = {
    try({
      boot_data <- resample_data(fit_obj$data, fit_obj$fit)
      fit_logistic(boot_data) # returns both data and fit
    })
  })
  
  failed2 <- sapply(res2, function(x) class(x) == "try-error")
  
  # Return
  if (any(failed2)) 
    warning("Some fits failed after 2 tries - result will contain fewer than the requested number of bootstraps")
  
  c(res[!failed], res2[!failed2])
}


# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=
# ---- Predictions --------------------------------------------------------------------------------
# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=
#' Extract predictions from the result of `fit_logistic()`
#' 
#' @param fit_obj Fitted model object which includes data.
#' @param level Grouping level at which predictions are sought (see details)
#' @param treatment_column Name of the column specifing unique treatments.
#' 
#' @details
#' By default, predictions are returned at the population level (`level = 0`).
#' Setting `level = 1` would give observation-level predictions. See 
#' `nlme::predict.nlme()` for details.
get_logistic_predictions <- function(fit_obj, level = 0, 
                                     treatment_column = "treatment") {
  if (level == 1) {
    fit_obj$data %>% 
      mutate(prediction = predict(fit_obj$fit, newdata = ., level = level))
    
  } else if(level == 0) {
    fit_obj$data %>% 
      group_by(.data[[treatment_column]], .data$virus, .data$time) %>% 
      summarise(P0 = mean(.data$P0), .groups = "drop") %>% 
      mutate(prediction = predict(fit_obj$fit, newdata = ., level = level))
    
  } else {
    stop("Predictions implemented for level=0 and level=1 only")
  }
}


#' Extract predictions from the result of `get_boot_fits()`
#' 
#' @param fit_obj Fitted model object which includes data.
#' @param level,treatment_column parameters passed to `get_logistic_predictions()`
get_logistic_predictions_boot <- function(boot_obj, level = 0, 
                                          treatment_column = "treatment") {
  lapply(boot_obj, get_logistic_predictions, 
         level = level, treatment_column = treatment_column) %>% 
    bind_rows()
}


#' Calculate prediction confidence intervals and add these to a prediction data frame
#' 
add_ci <- function(preds, boot_preds) {
  boot_preds %>% 
    bind_rows(.id = 'percent_cutoff') %>% 
    mutate(percent_cutoff = factor(.data$percent_cutoff, levels = names(cutoffs))) %>% 
    
    group_by(.data$percent_cutoff, .data$treatment, .data$virus, .data$time) %>% 
    summarise(lower = quantile(.data$prediction, probs = 0.025),
              upper = quantile(.data$prediction, probs = 0.975),
              P0_lower = quantile(.data$P0, probs = 0.025),
              P0_upper = quantile(.data$P0, probs = 0.975),
              .groups = "drop_last") %>% 
    ungroup() %>% 
    
    left_join(preds, by = c("percent_cutoff", "treatment", "virus", "time"))
}


# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=
## ---- COEFFICIENTS ------------------------------------------------------------------------------
# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=
#' Convert coefficient names into human-readable strings
#' 
#' @param coef_names A character vector of coefficient labels.
#' @param viruses Virus present in the fitted dataset (used to identify the baseline).
#' @param treatment_name Name to use for the "treatment" effect.
#' 
#' @return A character vector.
clean_coef_names <- function(coef_names, viruses, treatment_name) {
  # Determine baseline
  virus_additional <- coef_names %>% 
    str_extract("virus.+") %>% 
    str_remove("virus") %>% 
    na.omit() %>% 
    unique()
  
  stopifnot(length(virus_additional) == 1) # Only implemented for two viruses
  
  virus_baseline <- viruses[viruses != virus_additional]
  
  # Find matching component (baseline varies)
  component <- str_extract(coef_names, "^[[:alpha:]]+")
  
  stopifnot(all(component %in% c("r", "K")))
  r_baseline <- paste0("Baseline\n(", virus_baseline, ",\n", treatment_name, "=0)")
  k_baseline <- paste0("Baseline\n(", treatment_name, "=0)")
  
  # Process names
  clean_names <- coef_names %>% 
    str_remove("^[[:alpha:]]+") %>% 
    str_remove("^\\.") %>% 
    str_replace("\\(Intercept\\)", "Baseline") %>% 
    str_replace(":", " x ") %>% 
    str_to_sentence() %>% 
    str_replace("[T|t]reatment", treatment_name) %>% 
    str_replace("Virus.+", paste0("Virus (", virus_additional, ")")) %>% 
    str_replace("x virus.+", "x virus")
    
  clean_names <- if_else(component == "r", 
                         str_replace(clean_names, "^Baseline", r_baseline),
                         str_replace(clean_names, "^Baseline", k_baseline))
  
  clean_names
}


#' Extract coefficients and clean up names
#' 
#' @param fit_obj A fitted model from `fit_logistic()`
#' @param viruses Virus present in the fitted dataset (used to identify the baseline).
#' @param treatment_name Name to use for the "treatment" effect.
#' 
#' @return A data frame.
get_coefficients <- function(fit_obj, viruses, treatment_name = "Treatment") {
  coefficients(fit_obj$fit) %>% 
    summarise_all(.funs = mean) %>% 
    pivot_longer(everything(), names_to = "coefficient", values_to = "value") %>% 
    mutate(component = str_extract(.data$coefficient, "^[[:alpha:]]+"),
           component = case_when(.data$component == "r" ~ "Growth rate",
                                 .data$component == "K" ~ "Carrying capacity",
                                 TRUE ~ .data$component),
           effect = clean_coef_names(.data$coefficient, 
                                     viruses = viruses, 
                                     treatment_name = treatment_name),
           effect = if_else(grepl(":", .data$coefficient),
                            paste0(.data$effect, "\ninteraction"),
                            .data$effect))
}

#' Calculate coefficient 
#' 
#' @param boot_fits A list of bootstrap fits, as generated by `get_boot_fits()`.
#' @param viruses Virus present in the fitted dataset (used to identify the baseline).
#' @param treatment_name Name to use for the "treatment" effect.
#' 
#' @return A data frame.
get_coefficient_ci <- function(boot_fits, viruses, treatment_name = "Treatment") {
  lapply(boot_fits, get_coefficients, viruses = viruses, treatment_name = treatment_name) %>% 
    bind_rows() %>% 
    group_by(.data$coefficient, .data$component, .data$effect) %>% 
    summarise(lower = quantile(.data$value, probs = 0.025),
              upper = quantile(.data$value, probs = 0.975),
              .groups = "drop")
}


# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=
## ---- GROWTH RATE PREDICTIONS -------------------------------------------------------------------
# =-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=-=-=-=-==-=-=-=
# Calculating growth rates is dataset specific, so all functions below expect a dataset-specific
# function capable of returning growth rates for a single fit and treatment dose


#' Calculate a confidence interval for growth rate at a given dose
#' 
#' @param boot_fits bootstrap fits from a single cutoff, as produced by `get_boot_fits`.
#' @param treatment_dose a single dose.
#' @param growth_rate_fun a function which calculates growth rates (see details).
#' 
#' @details 
#' `growth_rate_fun` should have an argument named `treatment_dose`, and should 
#' accept one further argument - a single `fit_obj`. It should return a data frame
#' with (at least) columns "virus", "treatment" and "growth_rate". The "treatment"
#' column can be used to store the `treatment_dose` used for calculations, allowing
#' the results from calling this function with varying doses to be merged easily.
#' 
#' @return a data frame.
get_growth_rate_ci <- function(boot_fits, treatment_dose, growth_rate_fun) {
  lapply(boot_fits, growth_rate_fun, treatment_dose = treatment_dose) %>% 
    bind_rows() %>% 
    group_by(.data$virus, .data$treatment) %>% 
    summarise(lower = quantile(.data$growth_rate, probs = 0.025),
              upper = quantile(.data$growth_rate, probs = 0.975),
              .groups = "drop")
}


#' Get the predicted growth rates and matching confidence intervals for a range of doses
#' 
#' @param fit_obj a single model fit, as produced by `fit_logistic()`
#' @param boot_fits bootstrap fits from a single cutoff, as produced by `get_boot_fits`.
#' @param doses a vector of doses for which predictions are sought.
#' @param growth_rate_fun a function which calculates growth rates (see details).
#' 
#' @details 
#' `growth_rate_fun` should have an argument named `treatment_dose`, and should 
#' accept one further argument - a single `fit_obj`. It should return a data frame
#' with (at least) columns "virus", "treatment" and "growth_rate". The "treatment"
#' column can be used to store the `treatment_dose` used for calculations, allowing
#' the results from calling this function with varying doses to be merged easily.
#' 
#' @return a data frame.
combine_growth_rate_and_ci <- function(fit_obj, boot_fits, doses, growth_rate_fun) {
  main_est <- lapply(doses, growth_rate_fun, fit_obj = fit_obj) %>% 
    bind_rows()
  
  ci_est <- lapply(doses, get_growth_rate_ci, 
                   boot_fits = boot_fits, growth_rate_fun = growth_rate_fun) %>% 
    bind_rows()
  
  left_join(main_est, ci_est, by = c("virus", "treatment"))
}

