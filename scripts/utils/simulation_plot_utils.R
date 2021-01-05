# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Plotting functions for simulated data
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

require(dplyr)
require(tidyr)
require(stringr)
require(ggplot2)
require(scales)
require(cowplot)


#' Plot simulated population sizes over time
#' 
#' @param sims a data frame of population sizes over time at different parameter combinations
#' @param varying_param the key parameter that varies in these simulations (plotted along) 
#'        vertical axis
plot_simulation <- function(sims, varying_param = c("growth_rate", "impediment_scaling"), 
                            plot_theme = theme_bw()) {
  
  # Switch depending on parameter displayed:
  varying_param = match.arg(varying_param)
  
  if (varying_param == "growth_rate") {
    if (n_distinct(sims$impediment_scaling) != 1)
      stop("Multiple values for impediment_scaling detected while varying_param = 'growth_rate'")
    
    facetting <- facet_grid(rows = vars(scaling_factor), cols = vars(relative_impediment))
    horizontal_label <- expression(bold(textstyle("Growth rate inhibition"))%->%~"")
    vertical_label <- expression(bold(textstyle("Growth rate difference"))~(r[2]/r[1])%->%~"")
    
  } else {
    if (n_distinct(sims$growth_rate_1) != 1)
      stop("Multiple values for growth_rate detected while varying_param = 'impediment_scaling'")
    
    facetting <- facet_grid(rows = vars(impediment_scaling), cols = vars(relative_impediment))
    horizontal_label <- expression(bold(textstyle("Baseline growth rate inhibition"))%->%~"")
    vertical_label <- expression(bold(textstyle("Difference in inhibition strength"))%->%~"")
      
  }
  
  # Fix plotting order and labels:
  sims <- sims %>% 
    arrange(rev(.data$impediment_scaling), rev(.data$virus), .data$scaling_factor) %>% 
    mutate(virus = str_replace(.data$virus, "_", " "),
           virus = str_to_sentence(.data$virus),
           relative_impediment = sprintf("%3.0f%%", .data$relative_impediment * 100),
           scaling_factor = factor(.data$scaling_factor, 
                                   levels = rev(unique(.data$scaling_factor))),
           impediment_scaling = factor(.data$impediment_scaling, 
                                       levels = rev(unique(.data$impediment_scaling))))
  
  # Plot
  p <- ggplot(sims, aes(x = time, y = infected, shape = virus, colour = virus)) +
    geom_line(size = 0.6) +
    
    facetting +
    geom_rect(xmin = -Inf, ymin = -Inf, xmax = Inf, ymax = Inf,
              colour = "grey80", fill = NA) + 
    
    scale_y_log10(labels = label_number(big.mark = " "),
                  expand = expansion(mult = c(0.15, 0.1))) +
    scale_x_continuous(breaks = seq(min(sims$time), max(sims$time), by = 4),
                       labels = seq(min(sims$time), max(sims$time), by = 4)) +
    scale_colour_manual(values = VIRUS_COLOURS_SIMULATED, guide = FALSE) +
    
    labs(x = "Day", y = "Infected cells", shape = NULL, colour = NULL) +
    
    plot_theme +
    theme(plot.margin = margin(t = 20.5, r = 20.5, b = 5.5, l = 5.5))
  
  
  ggdraw(p) + 
    draw_label(horizontal_label, x = 0.55, y = 0.96) +
    draw_label(vertical_label, x = 0.98, y = 0.5, angle = -90)
}
