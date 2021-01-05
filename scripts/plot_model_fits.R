# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Plot raw data and model fits
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

library(scales)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggtext)

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


preds_full <- readRDS("output/model_fit/full_predictions.rds")
preds_reduced <- readRDS("output/model_fit/reduced_predictions.rds")


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


# ---- Plot data ----------------------------------------------------------------------------------
# Plotting each virus as a separate layer so points and lines of differnt colours don't overlap
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
  scale_colour_manual(values = VIRUS_COLOURS_DATA) +
  labs(x = "Day", y = " % GFP<sup>+</sup> cells", colour = NULL) +
  PRISM_THEME +
  theme(legend.position = "top",
        legend.margin = margin(t = 2.5, r = 2.5, b = -10, l = 2.5),
        axis.title.y = element_markdown(),
        strip.text = element_markdown())


# ---- Model fits ---------------------------------------------------------------------------------
plot_preds <- function(predictions, raw_data, labeller = dose_labeller) {
  ymax <- max(predictions$prediction, raw_data$infected_count)
  
  ggplot(predictions, aes(y = prediction, x = time, color = virus)) +
    # Raw data:
    geom_line(aes(y = infected_count, group = timeseries_id), size = 0.1, data = raw_data) +
    geom_point(aes(y = infected_count), shape = 20, size = 0.1, data = raw_data) +
    
    # Predictions:
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.8) +
    geom_line(size = 0.5) +
    geom_point(size = 0.5) +
    
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

p_full <- plot_preds(preds_full, data_raw)
p_reduced <- plot_preds(preds_reduced, data_raw)


# ---- Combine ------------------------------------------------------------------------------------
combined <- plot_grid(p_data, p_full, p_reduced,
                      nrow = 3, rel_heights = c(1.2, 1, 1),
                      align = "v", axis = "lr",
                      labels = LETTERS[1:3])


# ---- Output -------------------------------------------------------------------------------------
dir.create("output/plots", recursive = TRUE)
ggsave("output/plots/model_fits.pdf", combined, width = 7, height = 5)
