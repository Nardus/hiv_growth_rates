# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Fit logistic growth regression models to IFNa data
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

library(readr)

source("scripts/utils/logistic_growth_model_utils.R")

set.seed(20201203)

# ---- Data ---------------------------------------------------------------------------------------
data_raw <- read_csv("data/IFNa_data_2020_07_31.csv", 
                     col_types = cols(.default = col_double(),
                                      Experiment = col_factor(),
                                      Virus = col_factor(),
                                      SampleID = col_factor())) %>% 
  rename(infected_count = .data$CountInfectedCellsofLiveCells,
         virus = .data$Virus,
         treatment = .data$IFNaConc,
         replicate = .data$Replicate,
         time = .data$Timepoint) %>% 
  mutate(virus = factor(.data$virus, levels = c("CH058-TF", "CH058-6MO")),  # Set TF virus as baseline
         timeseries_id = paste(.data$virus, .data$treatment, .data$replicate, sep = "_")) # Uniquely identify each timeseries


# Since cells are eventually killed by the virus, only fitting to earlier points, while the virus
# is growing. This is determined by a tolerance cutoff, specifing how large of a decrease in the
# number of infected cells is allowed before discarding all subsequent points.
data_prepared <- prepare_data(data_raw, percent_cutoff = 30)


# ---- Model fits ---------------------------------------------------------------------------------
message("Fitting models")

# Maximum likelihood fit
model_full <- fit_logistic(data_prepared)
model_reduced <- fit_logistic(data_prepared, include_interaction = FALSE)


# Hierarchical bootstraps
n_boot <- 1000

boot_fits_full <- get_boot_fits(model_full, n = n_boot)
boot_fits_reduced <- get_boot_fits(model_reduced, n = n_boot)


# ---- Data needed for plots ----------------------------------------------------------------------
message("Extracting predictions")

# Predicted growth
preds_full <- prepare_predictions(model_full, boot_fits_full)
preds_reduced <- prepare_predictions(model_reduced, boot_fits_reduced)

# Coefficients
coefs_full <- prepare_coefs(model_full, boot_fits_full, 
                            raw_data = data_raw, 
                            treatment_name = "IFN")

coefs_reduced <- prepare_coefs(model_reduced, boot_fits_reduced, 
                               raw_data = data_raw, 
                               treatment_name = "IFN")


# Fitted growth rates
# Get growth rate at a given concentration (for a single cutoff):
virus_names <- unique(data_raw$virus)

get_growth_rate_full <- function(fit_obj, treatment_dose, viruses = virus_names) {
  get_coefficients(fit_obj, viruses = viruses) %>% 
    filter(.data$component == "Growth rate") %>% 
    group_by(.data$coefficient) %>% 
    summarise(value = mean(.data$value), .groups = "drop") %>% 
    pivot_wider(names_from = .data$coefficient, values_from = .data$value) %>% 
    transmute(`CH058-TF` = .data$`r.(Intercept)` + .data$r.treatment * treatment_dose,
              `CH058-6MO` = .data$`CH058-TF` + .data$`r.virusCH058-6MO` + +.data$`r.treatment:virusCH058-6MO` * treatment_dose) %>% 
    pivot_longer(everything(), names_to = 'virus', values_to = 'growth_rate') %>% 
    mutate(treatment = treatment_dose)
}

growth_rates_full <- combine_growth_rate_and_ci(model_full, boot_fits_full,
                                                doses = seq(0, 0.5, by = 0.1),
                                                growth_rate_fun = get_growth_rate_full)


# ---- Output -------------------------------------------------------------------------------------
dir.create("output/model_fit", recursive = TRUE)

saveRDS(data_raw, "output/model_fit/data.rds")

saveRDS(model_full, "output/model_fit/full_model.rds")
saveRDS(boot_fits_full, "output/model_fit/full_boot.rds")
saveRDS(preds_full, "output/model_fit/full_predictions.rds")
saveRDS(coefs_full, "output/model_fit/full_coefficients.rds")
saveRDS(growth_rates_full, "output/model_fit/full_growth_rates.rds")

saveRDS(model_reduced, "output/model_fit/reduced_model.rds")
saveRDS(boot_fits_reduced, "output/model_fit/reduced_boot.rds")
saveRDS(preds_reduced, "output/model_fit/reduced_predictions.rds")
saveRDS(coefs_reduced, "output/model_fit/reduced_coefficients.rds")