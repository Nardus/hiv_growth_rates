# Main workflow: generate simulations and fit logistic growth models

.PHONY: all
all: output/plots/simulations.pdf


# Simulations:
output/simulations/growth_rate_differs.rds:
	Rscript scripts/simulations.R


# Plots:
output/plots/simulations.pdf: output/simulations/growth_rate_differs.rds
	Rscript scripts/plot_simulations.R
