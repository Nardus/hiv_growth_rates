# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Plot fitted growth rates and coefficients
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

library(nlme)
library(scales)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggtext)
library(gridtext)
library(grid)

source("scripts/utils/plotting_constants.R")

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

growth_rates_full <- growth_rates_full %>% 
  mutate(virus = fix_virus_names(.data$virus),
         virus = factor(.data$virus, VIRUS_PLOTTING_ORDER))



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
  
  labs(x = "IFN dose (pg/µ*l*)", y = "Growth rate", colour = NULL) +
  PRISM_THEME +
  theme(axis.title.x = element_markdown(),
        legend.position = c(0.2, 0.2),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 18))



# ---- Combine ------------------------------------------------------------------------------------
combined_right <- plot_grid(p_model_table, p_growthrate,
                            nrow = 2, rel_heights = c(1, 2), 
                            labels = c("", "B"), hjust = -1)

combined_all <- plot_grid(p_all_coef, combined_right,
                          ncol = 2, rel_widths = c(1.2, 1),
                          labels = c("A", ""))



# ---- Output -------------------------------------------------------------------------------------
dir.create("output/plots", recursive = TRUE)

quartz(file = "output/plots/coefficients.pdf", type = "pdf", width = 7, height = 3)
print(combined_all)
dev.off()
