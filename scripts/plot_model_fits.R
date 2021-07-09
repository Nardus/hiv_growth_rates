# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Plot raw data, model fits and coefficients
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

library(nlme)
library(scales)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggtext)
library(gridtext)
library(grid)

source("scripts/utils/plotting_constants.R")

data_raw <- read_csv("data/IFNa_data_2020_07_31.csv", 
                     col_types = cols(.default = col_double(),
                                      Experiment = col_factor(),
                                      Virus = col_factor(),
                                      SampleID = col_factor())) %>% 
  rename(virus = .data$Virus,
         treatment = .data$IFNaConc,
         replicate = .data$Replicate,
         time = .data$Timepoint,
         percent_infected = .data$PerInfectedofLiveCells,
         infected_count = .data$CountInfectedCellsofLiveCells) %>% 
  mutate(timeseries_id = paste(.data$virus, .data$treatment, .data$replicate, sep = "_")) # Uniquely identify each timeseries


preds_reduced <- readRDS("output/model_fit/reduced_predictions.rds")

model_full <- readRDS("output/model_fit/full_model.rds")
model_reduced <- readRDS("output/model_fit/reduced_model.rds")

growth_rates_full <- readRDS("output/model_fit/full_growth_rates.rds")

coefs_full <- readRDS("output/model_fit/full_coefficients.rds")
coefs_reduced <- readRDS("output/model_fit/reduced_coefficients.rds")


# ---- Readable virus names -----------------------------------------------------------------------
fix_virus_names <- function(virus) {
  as.character(virus) %>% 
    str_replace_all("[-_]", " ") %>% 
    str_replace("6MO$", "CC")
}

data_raw <- data_raw %>% 
  mutate(virus = fix_virus_names(.data$virus))

preds_reduced <- preds_reduced %>% 
  mutate(virus = fix_virus_names(.data$virus))

growth_rates_full <- growth_rates_full %>% 
  mutate(virus = fix_virus_names(.data$virus),
         virus = factor(.data$virus, VIRUS_PLOTTING_ORDER))


# ---- Plot data ----------------------------------------------------------------------------------
# Plotting each virus as a separate layer so points and lines of different colours don't overlap
virus_1 <- filter(data_raw, .data$virus == VIRUS_PLOTTING_ORDER[1])
virus_2 <- filter(data_raw, .data$virus == VIRUS_PLOTTING_ORDER[2])

# Labeller to add dose units to facets:
dose_labeller <- function(dose) {
  paste(dose, "pg/µ*l* IFN")
}

dose_labeller <- as_labeller(dose_labeller)


p_data <- ggplot(data_raw, aes(x = time, y = percent_infected, colour = virus)) +
  geom_point(colour = "grey70", size = 0.3) +
  
  stat_summary(fun.data = mean_se, geom = "point", size = 0.5, data = virus_2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, data = virus_2) +
  stat_summary(fun.data = mean_se, geom = "line", size = 0.6, data = virus_2) +
  
  stat_summary(fun.data = mean_se, geom = "point", size = 0.5, data = virus_1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, data = virus_1) +
  stat_summary(fun.data = mean_se, geom = "line", size = 0.6, data = virus_1) +
  
  facet_wrap(vars(treatment), nrow = 1, scales = "free_y", labeller = dose_labeller) +
  
  scale_y_log10(limits = c(0.1, 100), 
                labels = label_number(big.mark = " ", accuracy = 0.1)) +
  scale_x_continuous(breaks = seq(1, 10, by = 2)) +
  scale_colour_manual(values = VIRUS_COLOURS_DATA,
                      guide = guide_legend(reverse = TRUE)) +
  labs(x = "Day", y = " % GFP<sup>+</sup> cells", colour = NULL) +
  PRISM_THEME +
  theme(legend.position = "top",
        legend.margin = margin(t = 2.5, r = 2.5, b = -10, l = 2.5),
        axis.title.y = element_markdown(),
        strip.text = element_markdown())


# ---- Model fits ---------------------------------------------------------------------------------
plot_preds <- function(predictions, raw_data, labeller = dose_labeller) {
  ymax <- max(predictions$prediction, raw_data$infected_count)
  
  raw_data_1 <- filter(raw_data, .data$virus == paste(VIRUS_PLOTTING_ORDER[1], "(raw data)"))
  raw_data_2 <- filter(raw_data, .data$virus == paste(VIRUS_PLOTTING_ORDER[2], "(raw data)"))
  
  predictions_1 <- filter(predictions, .data$virus == VIRUS_PLOTTING_ORDER[1])
  predictions_2 <- filter(predictions, .data$virus == VIRUS_PLOTTING_ORDER[2])
  
  ggplot(predictions, aes(y = prediction, x = time, color = virus)) +
    # Raw data:
    geom_line(aes(y = infected_count, group = timeseries_id), size = 0.1, data = raw_data_2) +
    geom_point(aes(y = infected_count), shape = 20, size = 0.1, data = raw_data_2) +
    
    geom_line(aes(y = infected_count, group = timeseries_id), size = 0.1, data = raw_data_1) +
    geom_point(aes(y = infected_count), shape = 20, size = 0.1, data = raw_data_1) +
    
    # Predictions:
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8, data = predictions_2) +
    geom_line(size = 0.5, data = predictions_2) +
    geom_point(size = 0.5, data = predictions_2) +
    
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8, data = predictions_1) +
    geom_line(size = 0.5, data = predictions_1) +
    geom_point(size = 0.5, data = predictions_1) +
    
    # Other layers:
    facet_wrap(vars(treatment), nrow = 1, scales = "free_y", labeller = labeller) +
    
    scale_y_log10(limits = c(10, ymax), labels = label_number(big.mark = " ")) +
    scale_x_continuous(breaks = seq(1, 10, by = 2)) +
    scale_colour_manual(values = c(VIRUS_COLOURS_DATA, VIRUS_COLOURS_FAINT),
                        guide = FALSE) +
    
    labs(x = "Day", y = "GFP<sup>+</sup> cells") +
    PRISM_THEME +
    theme(axis.title.y = element_markdown(),
          strip.text = element_markdown(),
          panel.spacing = unit(2, "points"))
}

data_raw <- data_raw %>% 
  mutate(virus = paste(.data$virus, "(raw data)"))

p_reduced <- plot_preds(preds_reduced, data_raw)


# ---- Coefficients -------------------------------------------------------------------------------
make_coef_plot <- function(coefs) {
  ggplot(coefs, aes(x = value, y = effect, xmin = lower, xmax = upper, colour = model)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey80") +
    geom_point(position = position_dodge2(width = 0.7, preserve = "single"), size = 0.8) +
    geom_errorbarh(position = position_dodge2(width = 1, preserve = "single"), height = 0.7) +
    
    facet_wrap(vars(component), nrow = 2, scales = "free") +
    scale_colour_manual(values = c("#AA3377", "#EE6677")) +
    scale_x_continuous(labels = label_number(big.mark = " "),
                       expand = expansion(mult = 0.12)) +
    
    labs(x = "Effect size", y = "Effect", colour = "Model") +
    PRISM_THEME
}

# Plotting order
effect_order <- c("Baseline", "Difference between\nviruses", 
                  "IFN sensitivity\n(1 pg/µl)",
                  "Additional IFN-sensitivity\nof CH058 CC (1 pg/µl)",
                  "IFN toxicity (1 pg/µl)")

combined_coefs <- list("Differential sensitivity" = coefs_full, 
                       "Constant sensitivity" = coefs_reduced) %>% 
  bind_rows(.id = "model") %>% 
  mutate(model = factor(.data$model, levels = c("Differential sensitivity", "Constant sensitivity")),
         effect = case_when(.data$effect == "IFN dose" ~ "IFN sensitivity\n(1 pg/µl)",
                            .data$effect == "Additional IFN-sensitivity\n(CH058 CC)" ~ "Additional IFN-sensitivity\nof CH058 CC (1 pg/µl)",
                            .data$effect == "IFN toxicity" ~ "IFN toxicity (1 pg/µl)",
                            TRUE ~ .data$effect),
         effect = factor(.data$effect, levels = rev(effect_order)))


# Plot
p_coef_r <- combined_coefs %>% 
  filter(.data$component == "Growth rate") %>% 
  make_coef_plot() +
  guides(colour = FALSE)

p_coef_k <- combined_coefs %>% 
  filter(.data$component == "Carrying capacity") %>% 
  make_coef_plot() +
  guides(colour = FALSE)


p_all_coef <- plot_grid(p_coef_r, p_coef_k, 
                        nrow = 2, rel_heights = c(2, 1.3), 
                        align = "v", axis = "lr")


# ---- Model comparison ---------------------------------------------------------------------------
# Use this as the legend for p_all_coef
legend_plot <- make_coef_plot(combined_coefs)
legend <- get_legend(legend_plot)

# Get and format summary stats
model_comparison <- anova(model_full$fit, model_reduced$fit)

lrt <- model_comparison[2, ]
lrt$df <- abs(diff(model_comparison$df))
lrt <- sprintf("χ<sup>2</sup> (%d) = %2.3f, p = %2.3f", 
               lrt$df, lrt$L.Ratio, lrt$`p-value`)

delta_aic <- sprintf("ΔAIC  = %2.3f", diff(model_comparison$AIC))

model_table <- c(lrt, delta_aic) %>% 
  paste("<span style='font-size:6pt'>", ., "</span>")

# Compose plot
p_model_table <- ggdraw(legend, xlim = c(0, 2)) +
  draw_line(x = 0.95, y = c(0.25, 0.56)) +
  draw_line(x = c(0.93, 0.9542), y = 0.25) +
  draw_line(x = c(0.93, 0.9542), y = 0.56) +
  draw_grob(richtext_grob(model_table,
                          x = 0.975, y = c(0.49, 0.37), 
                          hjust = 0, align_widths = TRUE))


# ---- Growth rates -------------------------------------------------------------------------------
p_growthrate <- ggplot(growth_rates_full, aes(x = treatment, y = growth_rate, colour = virus)) +
  geom_line(aes(group = virus), colour = "grey80", linetype = 2, size = 0.6) +
  
  geom_point(size = 0.8) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02) +
  
  scale_x_continuous(breaks = unique(growth_rates_full$treatment)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_colour_manual(values = VIRUS_COLOURS_DATA) +
  
  labs(x = "IFN dose (pg/µ*l*)", y = "Achieved growth rate", colour = NULL) +
  PRISM_THEME +
  theme(axis.title.x = element_markdown(),
        legend.position = c(0.2, 0.2),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 18))



# ---- Combine ------------------------------------------------------------------------------------
combined_top <- plot_grid(p_data, p_reduced,
                          nrow = 2, rel_heights = c(1.2, 1),
                          align = "v", axis = "lr",
                          labels = c("A", "B"))

combined_lower_right <- plot_grid(p_model_table, p_growthrate,
                                  nrow = 2, rel_heights = c(1, 2), 
                                  labels = c("", "D"), hjust = -1)

combined_lower <- plot_grid(p_all_coef, combined_lower_right,
                            ncol = 2, rel_widths = c(1.2, 1),
                            labels = c("C", ""))

combined_all <- plot_grid(combined_top, combined_lower, ncol = 1)


# ---- Output -------------------------------------------------------------------------------------
dir.create("output/plots", recursive = TRUE)

quartz(file = "output/plots/model_fits.pdf", type = "pdf", width = 7, height = 6.2)
print(combined_all)
dev.off()
