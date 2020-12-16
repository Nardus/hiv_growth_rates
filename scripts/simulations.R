# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Simulate logistic growth under different assumptions
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

library(dplyr)
library(tidyr)

source("scripts/utils/logistic_simulation_utils.R")

set.seed(20201203)


# ---- Only growth rates differ between viruses ---------------------------------------------------
param_combinations <- expand.grid(scaling_factor = seq(0.6, 0.95, by = 0.05),  # growth rates differ by up to 30%
                                  relative_impediment = seq(0, 0.9, by = 0.1),
                                  impediment_scaling = 1)  # 1 means impediment equally strong on both viruses

growth_rate_sims <- mapply(simulate_logistic,
                           scaling_factor = param_combinations$scaling_factor,
                           relative_impediment = param_combinations$relative_impediment,
                           impediment_scaling = param_combinations$impediment_scaling,
                           
                           MoreArgs = list(growth_rate_1 = 3,
                                           carrying_capacity = 1e5,
                                           times = seq(1, 10, by = 0.05)),
                           SIMPLIFY = FALSE) %>% 
  bind_rows() %>% 
  pivot_longer(starts_with("virus_"), names_to = "virus", values_to = "infected")



# ---- Only level of inhibition differs between viruses -------------------------------------------
param_combinations <- expand.grid(scaling_factor = 1,                        # growth rates equal
                                  relative_impediment = seq(0, 0.9, by = 0.1),
                                  impediment_scaling = seq(1.1, 1.8, by = 0.1))  # Impediment up to twice as strong on virus 2

inhibition_sims <- mapply(simulate_logistic,
                          scaling_factor = param_combinations$scaling_factor,
                          relative_impediment = param_combinations$relative_impediment,
                          impediment_scaling = param_combinations$impediment_scaling,
                          
                          MoreArgs = list(growth_rate_1 = 3,
                                          carrying_capacity = 1e5,
                                          times = seq(1, 10, by = 0.05)),
                          SIMPLIFY = FALSE) %>% 
  bind_rows() %>% 
  pivot_longer(starts_with("virus_"), names_to = "virus", values_to = "infected")



# ---- Output -------------------------------------------------------------------------------------
dir.create("output/simulations", recursive = TRUE)

saveRDS(growth_rate_sims, "output/simulations/growth_rate_differs.rds")
saveRDS(inhibition_sims, "output/simulations/inhibition_level_differs.rds")
