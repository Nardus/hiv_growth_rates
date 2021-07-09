# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Utility functions for simulating logistic growth
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

require(dplyr)


#' Logistic growth function
#' 
#' @param times vector of times for which population sizes should be returned
#' @param growth_rate the growth rate
#' @param carrying_capacity maximum carrying capacity
#' @param starting_pop population size at time = 0
#' 
#' @return a vector of population sizes
logistic_fun <- function(times, growth_rate, carrying_capacity, starting_pop) {
  (carrying_capacity * starting_pop * exp(growth_rate*times)) / (carrying_capacity + starting_pop * (exp(growth_rate*times) - 1))
}


#' Simulate logistic growth of two viruses under influence of something impeding growth
#' 
#' @param times vector of times for which population sizes should be returned
#' @param growth_rate_1 the growth rate of virus 1
#' @param scaling_factor scaling factor for virus 2's growth rate relative to that of virus 1: `r2 = r1 + s`
#' @param relative_impediment size of growth rate impediment (as proportion of `growth_rate_1`), e.g.
#'                            0.4 means decrease `growth_rate_1` by 40%
#' @param impediment_scaling scaling factor for virus 2's impediment. For example, `impediment_scaling = 2`
#'                           would mean virus 2 experiences a growth impediment twice as large as that of 
#'                           virus 1
#' @param carrying_capacity carrying capacity for a single virus (virus 1 and virus 2 are assumed to be 
#'                          growing in separate environments, they are not competing)
#' @param starting_pop number of initially infected cells
#' 
#' @return a data frame of population sizes
#' 
simulate_logistic <- function(times, growth_rate_1, scaling_factor, relative_impediment, impediment_scaling, 
                              carrying_capacity, starting_pop = 100) {
  # Vectorised over times only:
  stopifnot(length(growth_rate_1) == 1)
  stopifnot(length(scaling_factor) == 1)
  stopifnot(length(relative_impediment) == 1)
  stopifnot(length(starting_pop) == 1)
  
  actual_impediment = relative_impediment * growth_rate_1 # Make impediment size relative
  
  r1 = growth_rate_1 - actual_impediment
  r2 = growth_rate_1*scaling_factor - actual_impediment*impediment_scaling
  
  r1 = max(0, r1)
  r2 = max(0, r2)
  
  pop_1 <- logistic_fun(times, growth_rate = r1, 
                        carrying_capacity = carrying_capacity, starting_pop = starting_pop)
  
  pop_2 <- logistic_fun(times, growth_rate = r2, 
                        carrying_capacity = carrying_capacity, starting_pop = starting_pop)
  
  # Ensure no negative pop sizes:
  pop_1 <- if_else(pop_1 < 0, 0, pop_1)
  pop_2 <- if_else(pop_2 < 0, 0, pop_2)
  
  data.frame(time = c(0, times),
             growth_rate_1 = growth_rate_1,
             scaling_factor = scaling_factor,
             relative_impediment = relative_impediment,
             impediment_scaling = impediment_scaling,
             carrying_capacity = carrying_capacity,
             virus_1 = c(starting_pop, pop_1),
             virus_2 = c(starting_pop, pop_2))
}
