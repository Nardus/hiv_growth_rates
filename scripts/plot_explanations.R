# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Plot explanation for the effect of growth rate difference
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

library(cowplot)
library(ggtext)
library(scales)
library(dplyr)
library(tidyr)
library(readr)

source("scripts/utils/simulation_plot_utils.R")
source("scripts/utils/plotting_constants.R")

growth_rate_sims <- readRDS("output/simulations/growth_rate_differs.rds")

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


# ---- Readable virus names -----------------------------------------------------------------------
fix_virus_names <- function(virus) {
  as.character(virus) %>% 
    str_replace_all("[-_]", " ") %>% 
    str_replace("6MO$", "CC")
}

data_raw <- data_raw %>% 
  mutate(virus = fix_virus_names(.data$virus))


# ---- Simulations --------------------------------------------------------------------------------
sims <- growth_rate_sims %>% 
  arrange(rev(.data$impediment_scaling), rev(.data$virus), .data$scaling_factor) %>% 
  mutate(time = time + 1,
         virus = str_replace(.data$virus, "_", " "),
         virus = str_to_sentence(.data$virus),
         relative_impediment = sprintf("%3.0f%%", .data$relative_impediment * 100),
         scaling_factor = factor(.data$scaling_factor, 
                                 levels = rev(unique(.data$scaling_factor))),
         impediment_scaling = factor(.data$impediment_scaling, 
                                     levels = rev(unique(.data$impediment_scaling)))) %>% 
  
  filter(.data$relative_impediment %in% c("  0%", " 70%")) %>% 
  filter(.data$scaling_factor == "0.8")


sim_plot <- ggplot(sims, aes(x = time, y = infected/10000*100, shape = virus, colour = virus)) +
  geom_line(size = 0.6) +
  
  facet_grid(cols = vars(relative_impediment)) + 
  
  scale_y_continuous(limits = c(0, 100), labels = label_number(big.mark = " ", accuracy = 0.1)) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  scale_colour_manual(values = VIRUS_COLOURS_SIMULATED, guide = FALSE) +
  
  labs(x = "Day", y = "% infected", shape = NULL, colour = NULL) +
  
  PRISM_THEME +
  theme(plot.margin = margin(t = 20.5, r = 20.5, b = 5.5, l = 5.5))


# ---- Observed data ----------------------------------------------------------------------------------
data_raw <- data_raw %>% 
  filter(treatment %in% c(0, 0.5))

# Plotting each virus as a separate layer so points and lines of different colours don't overlap
virus_1 <- filter(data_raw, .data$virus == VIRUS_PLOTTING_ORDER[1])
virus_2 <- filter(data_raw, .data$virus == VIRUS_PLOTTING_ORDER[2])

# Labeller to add dose units to facets:
dose_labeller <- function(dose) {
  paste(dose, "pg/µ*l* IFN")
}

dose_labeller <- as_labeller(dose_labeller)


data_plot <- ggplot(data_raw, aes(x = time, y = percent_infected, colour = virus)) +
  geom_point(colour = "grey70", size = 0.3) +
  
  stat_summary(fun.data = mean_se, geom = "point", size = 0.5, data = virus_2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, data = virus_2) +
  stat_summary(fun.data = mean_se, geom = "line", size = 0.6, data = virus_2) +
  
  stat_summary(fun.data = mean_se, geom = "point", size = 0.5, data = virus_1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.8, data = virus_1) +
  stat_summary(fun.data = mean_se, geom = "line", size = 0.6, data = virus_1) +
  
  facet_wrap(vars(treatment), nrow = 1, labeller = dose_labeller) +
  
  scale_y_continuous(limits = c(0, 100), labels = label_number(big.mark = " ", accuracy = 0.1)) +
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  scale_colour_manual(values = VIRUS_COLOURS_DATA, guide = FALSE) +
  labs(x = "Day", y = " % GFP<sup>+</sup> cells", colour = NULL) +
  PRISM_THEME +
  theme(legend.position = "top",
        legend.margin = margin(t = 2.5, r = 2.5, b = -10, l = 2.5),
        axis.title.y = element_markdown(),
        strip.text = element_markdown())


# ---- Combine ------------------------------------------------------------------------------------
combined <- plot_grid(sim_plot, data_plot, ncol = 1,
                      align = "h", axis = "lr",
                      labels = c("A", "B"), label_size = 10)

dir.create("output/plots/", recursive = TRUE)
ggsave2("output/plots/explanations.pdf", combined, width = 6, height = 6)
