# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Constants used across all plots to ensure a uniform look
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

require(ggplot2)
require(colorspace)

# Theme designed to match plots from Graphpad Prism:
PRISM_THEME <- theme_classic(base_size = 9, 
                             base_family = "sans", 
                             base_line_size = 0.6,
                             base_rect_size = 0.6) +
  theme(strip.background = element_blank(),
        plot.title = element_text(colour = "grey10", face = "bold"),
        plot.subtitle = element_text(colour = "grey10", face = "bold"),
        axis.text = element_text(size = 8, colour = "grey10", face = "bold"),
        axis.title = element_text(colour = "grey10", face = "bold"),
        strip.text = element_text(size = 8, face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.line = element_line(colour = "grey10"),
        legend.key.size = unit(1, "lines"),
        legend.key.height = unit(0.8, "lines"))


# Colour schemes:
VIRUS_COLOURS_SIMULATED <- c("Virus 1" = "#3366CC", "Virus 2" = "#FF9933")
VIRUS_COLOURS_DATA <- c("CH058 TF" = "#3366CC", "CH058 CC" = "#FF9933")

VIRUS_PLOTTING_ORDER <- c("CH058 TF", "CH058 CC")

# Fainter version of these colours
VIRUS_COLOURS_FAINT <- lighten(VIRUS_COLOURS_DATA, amount = 0.6)
names(VIRUS_COLOURS_FAINT) <- paste(VIRUS_PLOTTING_ORDER, "(in vitro)")