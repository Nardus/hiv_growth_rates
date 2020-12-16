# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Constants used across all plots to ensure a uniform look
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Theme designed to match plots from Graphpad Prism:
PRISM_THEME <- theme_classic(base_size = 12, 
                             base_family = "sans", 
                             base_line_size = 0.6,
                             base_rect_size = 0.6) +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 8, face = "bold"),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))


# Colour schemes:
VIRUS_COLOURS_SIMULATED <- c("Virus 1" = "#0077bb", "Virus 2" = "#ee7733")
VIRUS_COLOURS_DATA <- c("CH058 TF" = "#0077bb", "CH058 CC" = "#ee7733")