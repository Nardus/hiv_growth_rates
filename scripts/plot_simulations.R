# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Plot simulations
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

library(cowplot)

source("scripts/utils/simulation_plot_utils.R")
source("scripts/utils/plotting_constants.R")

growth_rate_sims <- readRDS("output/simulations/growth_rate_differs.rds")
inhibition_sims <- readRDS("output/simulations/inhibition_level_differs.rds")


p1 <- plot_simulation(growth_rate_sims, "growth_rate", plot_theme = PRISM_THEME)
p2 <- plot_simulation(inhibition_sims, "impediment_scaling", plot_theme = PRISM_THEME)


combined <- plot_grid(p1, p2, nrow = 2,
                      labels = c("A", "B"))

dir.create("output/plots/", recursive = TRUE)
ggsave2("output/plots/simulations.pdf", combined, width = 7, height = 8)



