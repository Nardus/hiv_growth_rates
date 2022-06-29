# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Plot raw data, model fits and coefficients
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

library(nlme)
library(scales)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggtext)
library(gridtext)
library(grid)

source("scripts/utils/plotting_constants.R")
source("scripts/utils/logistic_growth_model_utils.R")

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


preds_full <- readRDS("output/model_fit/full_predictions.rds")
preds_reduced <- readRDS("output/model_fit/reduced_predictions.rds")

model_full <- readRDS("output/model_fit/full_model.rds")
model_reduced <- readRDS("output/model_fit/reduced_model.rds")

growth_rates_full <- readRDS("output/model_fit/full_growth_rates.rds")

coefs_full <- readRDS("output/model_fit/full_coefficients.rds")
coefs_reduced <- readRDS("output/model_fit/reduced_coefficients.rds")

boot_full <- readRDS("output/model_fit/full_boot.rds")


# ---- Readable virus names -----------------------------------------------------------------------
fix_virus_names <- function(virus) {
  as.character(virus) %>% 
    str_replace_all("[-_]", " ") %>% 
    str_replace("6MO$", "CC")
}

data_raw <- data_raw %>% 
  mutate(virus = fix_virus_names(.data$virus))

preds_full <- preds_full %>% 
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
  labs(x = "Day", y = " % GFP<sup>+</sup> cells", colour = NULL,
       title = "In vitro") +
  PRISM_THEME +
  theme(legend.margin = margin(t = 2.5, r = 2.5, b = -10, l = 2.5),
        axis.title.y = element_markdown(),
        strip.text = element_markdown())


# ---- Model fits ---------------------------------------------------------------------------------
plot_preds <- function(predictions, raw_data, labeller = dose_labeller) {
  ymax <- max(predictions$prediction, raw_data$infected_count)
  virus_colours_data <- VIRUS_COLOURS_DATA
  names(virus_colours_data) <- paste(names(virus_colours_data), "(model)")
  legend_order <- c(names(virus_colours_data)[1], names(VIRUS_COLOURS_FAINT)[1],
                    names(virus_colours_data)[2], names(VIRUS_COLOURS_FAINT)[2])
  
  raw_data_1 <- filter(raw_data, .data$virus == paste(VIRUS_PLOTTING_ORDER[1], "(in vitro)"))
  raw_data_2 <- filter(raw_data, .data$virus == paste(VIRUS_PLOTTING_ORDER[2], "(in vitro)"))
  
  predictions_1 <- filter(predictions, .data$virus == paste(VIRUS_PLOTTING_ORDER[1], "(model)"))
  predictions_2 <- filter(predictions, .data$virus == paste(VIRUS_PLOTTING_ORDER[2], "(model)"))
  
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
    scale_colour_manual(breaks = legend_order,
                        values = c(virus_colours_data, VIRUS_COLOURS_FAINT)) +
    
    labs(x = "Day", y = "GFP<sup>+</sup> cells", colour = NULL) +
    PRISM_THEME +
    theme(axis.title.y = element_markdown(),
          strip.text = element_markdown(),
          panel.spacing = unit(2, "points"))
}

data_raw <- data_raw %>% 
  mutate(virus = paste(.data$virus, "(in vitro)"))

preds_full <- preds_full %>% 
  mutate(virus = paste(.data$virus, "(model)"))

preds_reduced <- preds_reduced %>% 
  mutate(virus = paste(.data$virus, "(model)"))


p_full <- plot_preds(preds_full, data_raw) +
  ggtitle("Differential sensitivity model fit")

p_reduced <- plot_preds(preds_reduced, data_raw) +
  ggtitle("Constant sensitivity model fit")


# ---- Coefficients -------------------------------------------------------------------------------
make_coef_plot <- function(coefs) {
  ggplot(coefs, aes(x = value, y = effect, xmin = lower, xmax = upper, colour = model)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey80") +
    geom_point(position = position_dodge2(width = 0.7, preserve = "single"), size = 0.8) +
    geom_errorbarh(position = position_dodge2(width = 1, preserve = "single"), height = 0.7) +
    
    scale_colour_manual(values = c("#AA3377", "#EE6677")) +
    scale_x_continuous(labels = label_number(big.mark = " "),
                       expand = expansion(mult = 0.12)) +
    
    labs(x = "Effect size", y = "Effect", colour = "Model") +
    PRISM_THEME
}

# Add 0 for removed effect in constant sensitivity model
constant_effect <- data.frame(coefficient = "r.treatment:virusCH058-6MO",
                              component = "Growth rate",
                              effect = "Additional IFN-sensitivity\n(CH058 CC)",
                              value = 0,
                              lower = NA,
                              upper = NA)
coefs_reduced <- bind_rows(coefs_reduced, constant_effect)

# Plotting order
effect_order <- c("Common growth rate",
                  "Baseline carrying capacity",
                  "Reduction in growth\nrate for CC", 
                  "Common IFN sensitivity\n(0.5 pg/µl)",
                  "Additional IFN sensitivity\nfor CC (0.5 pg/µl)",
                  "IFN toxicity (0.5 pg/µl)")

# Coefficients associated with dose (showing 0.5 instead 1pg/ul dose effect, matching experiment range)
dose_coef_labels <- c("r.treatment", "r.treatment:virusCH058-6MO", "K.treatment")

# Convert doses and set labels
rename_effects <- function(data) {
  data %>% 
    mutate(effect = case_when(.data$effect == "Baseline" & .data$component == "Growth rate" ~ "Common growth rate",
                              .data$effect == "Baseline" ~ "Baseline carrying capacity",
                              .data$effect == "Difference between\nviruses" ~ "Reduction in growth\nrate for CC",
                              .data$effect == "IFN dose" ~ "Common IFN sensitivity\n(0.5 pg/µl)",
                              .data$effect == "Additional IFN-sensitivity\n(CH058 CC)" ~ "Additional IFN sensitivity\nfor CC (0.5 pg/µl)",
                              .data$effect == "IFN toxicity" ~ "IFN toxicity (0.5 pg/µl)",
                              TRUE ~ .data$effect))
}

combined_coefs <- list("Differential sensitivity" = coefs_full, 
                       "Constant sensitivity" = coefs_reduced) %>% 
  bind_rows(.id = "model") %>% 
  mutate(value = if_else(.data$coefficient %in% dose_coef_labels, .data$value * 0.5, .data$value),  # Convert dose
         lower = if_else(.data$coefficient %in% dose_coef_labels, .data$lower * 0.5, .data$lower),
         upper = if_else(.data$coefficient %in% dose_coef_labels, .data$upper * 0.5, .data$upper)) %>% 
  rename_effects() %>% 
  mutate(model = factor(.data$model, levels = c("Differential sensitivity", "Constant sensitivity")),
         effect = factor(.data$effect, levels = rev(effect_order)))


# Plot
p_coef_r <- combined_coefs %>% 
  filter(.data$component == "Growth rate") %>% 
  make_coef_plot() +
  geom_text(aes(label = "(fixed)"), nudge_y = 0.2, hjust = -0.1, size = 2,
            data = filter(combined_coefs, .data$value == 0)) +
  labs(subtitle = "Fitted growth rate effects") +
  guides(colour = FALSE)

p_coef_k <- combined_coefs %>% 
  filter(.data$component == "Carrying capacity") %>% 
  make_coef_plot() +
  labs(subtitle = "Fitted carrying capacity effects") +
  guides(colour = FALSE)



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
p_model_table <- ggdraw(legend, xlim = c(0, 2.5)) +
  draw_line(x = 0.88, y = c(0.25, 0.56)) +
  draw_line(x = c(0.86, 0.8842), y = 0.23) +
  draw_line(x = c(0.86, 0.8842), y = 0.56) +
  draw_grob(richtext_grob(model_table,
                          x = 0.905, y = c(0.51, 0.35), 
                          hjust = 0, align_widths = TRUE))


# ---- Growth rates -------------------------------------------------------------------------------
growth_rates_full$plot_label <- "Fitted growth rate\n"

p_growthrate <- ggplot(growth_rates_full, aes(x = treatment, y = growth_rate, colour = virus)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey20") +
  
  geom_line(aes(group = virus), colour = "grey80", linetype = 2, size = 0.6) +
  
  geom_point(size = 0.8) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.02) +
  
  scale_x_continuous(breaks = unique(growth_rates_full$treatment)) +
  scale_y_continuous(breaks = seq(-3, 3, by = 1), limits = c(-3, 3)) +
  scale_colour_manual(values = VIRUS_COLOURS_DATA) +
  
  facet_wrap(vars(plot_label)) +
  
  labs(x = "IFN dose (pg/µl)", y = "Achieved growth rate", colour = NULL) +
  PRISM_THEME +
  theme(legend.position = c(0.8, 0.9),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 18))


# Inset: illustrate expected patterns
treatments <- unique(data_raw$treatment)
baseline <- with(coefs_full, value[coefficient == "r.(Intercept)"])
treat_effect <- with(coefs_full, value[coefficient == "r.treatment"])
virus_effect <- with(coefs_full, value[coefficient == "r.virusCH058-6MO"]) - 0.1 # larger, for clarity
extra_effect <- 2.5 # Hypothetical, for illustration

expectations <- bind_rows(
  data.frame(label = "Same\nsensitivity",
             treatment = treatments,
             TF = baseline + treat_effect*treatments,
             CC = baseline + treat_effect*treatments + virus_effect),
  data.frame(label = "CC more\nsensitive",
             treatment = treatments,
             TF = baseline + treat_effect*treatments,
             CC = baseline + treat_effect*treatments + virus_effect + -extra_effect*treatments),
  data.frame(label = "CC less\nsensitive",
             treatment = treatments,
             TF = baseline + treat_effect*treatments,
             CC = baseline + treat_effect*treatments + virus_effect + extra_effect*treatments)
) %>% 
  pivot_longer(c("TF", "CC"), names_to = "virus", values_to = "growth_rate") %>% 
  mutate(label = factor(.data$label, levels = c("Same\nsensitivity", 
                                                "CC more\nsensitive",
                                                "CC less\nsensitive")))

inset_colours <- c("TF" = unname(VIRUS_COLOURS_DATA["CH058 TF"]),
                   "CC" = unname(VIRUS_COLOURS_DATA["CH058 CC"]))

inset_plot <- ggplot(expectations, aes(x = treatment, y = growth_rate, colour = virus)) +
  geom_line() +
  facet_grid(cols = vars(label)) +
  scale_colour_manual(values = inset_colours) +
  ggtitle("Expected patterns") +
  theme_nothing() +
  theme(plot.title = element_text(size = 6, colour = "grey10", face = "bold"),
        strip.text = element_text(size = 5, colour = "grey10", face = "bold",
                                  margin = margin(t = 3, r = 3, b = 3, l = 3)),
        plot.margin = margin(t = 3.6, r = 3.6, b = 3.6, l = 3.6),
        panel.spacing = unit(0.1, "lines"),
        plot.background = element_rect(colour = "grey20", fill = NA, size = 0.6))

p_growthrate <- p_growthrate +
  annotation_custom(ggplotGrob(inset_plot), 
                    xmin = -0.01, xmax = 0.36,
                    ymin = -0.5, ymax = -3.05)


# ---- Illustrate individual effects making up growth rate ----------------------------------------
doses <- unique(data_raw$treatment) %>% 
  expand.grid(coefficient = coefs_full$coefficient,
              virus = factor(c("CH058-TF", "CH058-6MO"), levels = c("CH058-TF", "CH058-6MO")),
              dose = .) %>% 
  mutate(virus_binary = as.numeric(.data$virus) - 1)

calculate_effect <- function(fit_obj, doses_df = doses) {
  effect_order <- c("Common growth\nrate",
                    "Baseline carrying capacity",
                    "Reduction in growth\nrate for CC", 
                    "Common IFN\nsensitivity",
                    "Additional IFN\nsensitivity for CC",
                    "IFN toxicity")

  get_coefficients(fit_obj, viruses = c("CH058-6MO", "CH058-TF"), treatment_name = "IFN") %>% 
    left_join(doses, by = "coefficient") %>% 
    mutate(value = if_else(grepl("treatment", .data$coefficient), 
                           .data$value * .data$dose,
                           .data$value),
           value = if_else(grepl("virus", .data$coefficient),
                           .data$value * .data$virus_binary,
                           .data$value)) %>% 
    mutate(effect = case_when(.data$effect == "Baseline" & .data$component == "Growth rate" ~ "Common growth\nrate",
                              .data$effect == "Baseline" ~ "Baseline carrying capacity",
                              .data$effect == "Difference between\nviruses" ~ "Reduction in growth\nrate for CC",
                              .data$effect == "IFN dose" ~ "Common IFN\nsensitivity",
                              .data$effect == "Additional IFN-sensitivity\n(CH058 CC)" ~ "Additional IFN\nsensitivity for CC",
                              .data$effect == "IFN toxicity" ~ "IFN toxicity",
                              TRUE ~ .data$effect),
           effect = factor(.data$effect, levels = effect_order))
}

effects_ci <- lapply(boot_full, calculate_effect) %>% 
  bind_rows() %>% 
  group_by(.data$component, .data$coefficient, .data$effect, .data$virus, .data$dose) %>% 
  summarise(lower = quantile(.data$value, probs = 0.025),
            upper = quantile(.data$value, probs = 0.975),
            .groups = "drop")

effects_ml <- calculate_effect(model_full) %>% 
  left_join(effects_ci, by = c("coefficient", "component", "effect", "virus", "dose")) %>% 
  mutate(virus = fix_virus_names(.data$virus),
         virus = factor(.data$virus, levels = VIRUS_PLOTTING_ORDER))

growth_effects <- effects_ml %>% 
  filter(.data$component == "Growth rate")

# Make each plot separately so we can add annotations between them:
plot_effect <- function(effect_name, effect_data, y_label) {
  effect_data <- effect_data %>% 
    filter(.data$effect == effect_name)
  
  effect_data_TF <- effect_data %>% 
    filter(grepl("TF", .data$virus))
  
  effect_data_CC <- effect_data %>% 
    filter(grepl("CC", .data$virus))
  
  ggplot(effect_data, aes(x = dose, y = value, colour = virus)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey20") +
    
    geom_line(aes(group = virus), colour = "grey80", linetype = 2, size = 0.6) +
    
    geom_point(size = 0.8, data = effect_data_CC) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.04, data = effect_data_CC) +
    geom_point(size = 0.8, data = effect_data_TF) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.04, data = effect_data_TF) +
    
    facet_wrap(vars(effect), scales = "free", nrow = 1) +
    
    scale_y_continuous(breaks = seq(-3, 3, by = 1), limits = c(-3, 3)) +
    scale_colour_manual(values = VIRUS_COLOURS_DATA, guide = FALSE) +
    
    labs(x = "IFN dose (pg/µl)", y = y_label, colour = NULL) +
    PRISM_THEME
}

p_growth_effects <- lapply(sort(unique(growth_effects$effect)), plot_effect,
                           effect_data = growth_effects,
                           y_label = "Effect on growth rate")


p_growth_effects <- plot_grid(plotlist = c(p_growth_effects, list(p_growthrate)), 
                              nrow = 1, rel_widths = c(1, 1, 1, 1, 2),
                              labels = c("", "+", "+", "+", "="), hjust = -0.5, vjust = 2,
                              label_size = 10)



# ---- Combine ------------------------------------------------------------------------------------
combined_top <- plot_grid(p_data, p_full, p_reduced,
                          nrow = 3, rel_heights = c(1.2, 1, 1),
                          align = "v", axis = "lr",
                          labels = c("A", "B", "C"), label_size = 10)

combined_coef <- plot_grid(p_model_table, p_coef_k,
                           nrow = 2, rel_heights = c(1.3, 2))

combined_middle <- plot_grid(p_coef_r, combined_coef, 
                             nrow = 1, rel_heights = c(2, 1.3))

combined_all <- plot_grid(combined_top, combined_middle, p_growth_effects, 
                          ncol = 1, rel_heights = c(3, 1.2, 1.5),
                          labels = c("", "D", "E"), label_size = 10)


# ---- Output -------------------------------------------------------------------------------------
dir.create("output/plots", recursive = TRUE)

cairo_pdf(file = "output/plots/model_fits2.pdf", width = 7, height = 8)
print(combined_all)
dev.off()
