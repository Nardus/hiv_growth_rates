# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Fit logistic growth regression models to IFNa data
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(nlme)

source("scripts/utils/logistic_growth_model_utils.R")

set.seed(20201203)

# ---- Data ---------------------------------------------------------------------------------------
data_ifn <- read_csv("data/IFNa_data_2020_07_31.csv", 
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
data_ifn <- prepare_data(data_ifn, percent_cutoff = 30)


# ---- Model fits ---------------------------------------------------------------------------------
# Maximum likelihood fit
model_full <- fit_logistic(data_ifn)
model_reduced <- fit_logistic(data_ifn, include_interaction = FALSE)


# Hierarchical bootstraps
n_boot <- 1000

boot_fits_full <- get_boot_fits(model_full, n = n_boot)
boot_fits_reduced <- get_boot_fits(model_reduced, n = n_boot)


# ---- Output -------------------------------------------------------------------------------------
dir.create("output/model_fit", recursive = TRUE)

saveRDS(model_full, "output/model_fit/full_model.rds")
saveRDS(model_reduced, "output/model_fit/reduced_model.rds")

saveRDS(boot_fits_full, "output/model_fit/full_boot.rds")
saveRDS(boot_fits_reduced, "output/model_fit/reduced_boot.rds")
