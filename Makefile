# Main workflow: generate simulations and fit logistic growth models

.PHONY: all
all: output/plots/simulations.pdf \
     output/plots/model_fits.pdf \
	 output/plots/explanations.pdf


# Simulations:
output/simulations/growth_rate_differs.rds:
	Rscript scripts/simulations.R


output/plots/simulations.pdf: output/simulations/growth_rate_differs.rds
	Rscript scripts/plot_simulations.R


# Model fits:
output/model_fit/full_model.rds: data/IFNa_data_2020_07_31.csv
	Rscript scripts/fit_growth_models.R


output/plots/model_fits.pdf: output/model_fit/full_model.rds
	Rscript scripts/plot_model_fits.R


output/plots/explanations.pdf: output/simulations/growth_rate_differs.rds data/IFNa_data_2020_07_31.csv
	Rscript scripts/plot_explanations.R
