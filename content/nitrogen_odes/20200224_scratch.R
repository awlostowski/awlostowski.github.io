remove(list = ls()) # clear all workspace variables
cat("\014")         # clear command line

library(rstudioapi)
library(deSolve)
library(ggplot2)
library(dplyr)
library(cowplot); theme_set(theme_cowplot())

# Set working directory to source file location
source_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(source_path))

# Scratch code - February 24, 2020
# A system of ODE simulating inorganic nitrogen cycling in an antarctic river.


######### LOAD DATA ######### 
# load the obsrevations - shout out Tyler K.
dat <- read.csv("relict_channel_data.csv")

######### INITIALIZE ######### 
# initialize the model domain, channel distance in meters
dist = seq(1530,4000,by = 1)

# initial conditions, N-species compositions at the upper boundary condition. 
init <- c(Ng = 10,  # concentration of nigrogen in glacier ice (ug/L)
          Nb = 1)   # concentration of black-mat derived N at x = 0

# end members: delN-15 concentrations in black mats and glacier ice
n15g = -10; 
n15b = 0;   

# ??? Why are these end members different - refresh on Tyler's paper. 


######### BUILD THE MODEL ######### 
# The model equations embedded in a function
model <- function(distance, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    
    # system of equations
    # dNg <- -gamma * (Ng)           # Change is glacier-derived nitrogen with distance
    # dNb <- phi - gamma * (Nb)    # change in black mat derived nitrogen with distance
    
    dNg <- -gamma * (Ng^2/(Ng + Nb))           # Change is glacier-derived nitrogen with distance
    dNb <- phi - (gamma * (Nb^2/(Nb + Ng)))    # change in black mat derived nitrogen with distance
    
    list(c(dNg, dNb))
  })
}

######### DIAGNOSTICS ######### 
###############################################################
# # fiddle stix
# parameters <- c(0.02,0.01)
# names(parameters) <- c("phi", "gamma")
# 
# test <- ode(y = init, times = dist, func = model, parms = parameters)
# test <- data.frame(test)
# 
# test <- test %>%
#   mutate(n15 = ((n15g * Ng) + (n15b * Nb))/(Ng + Nb))
# 
# p <- ggplot(test, aes(x = time, y = n15)) +
#   geom_line(color = "blue")
#   # geom_line(aes(x = time, y = Nb), color = "black")
# 
# print(p)
#################################################################
 

######### OPTIMIZATION #########    
# a function for the objective function - the optimization routine will minimize this. 
RSS <- function(parameters) {
  
  names(parameters) <- c("phi", "gamma")
  
  # run the ode system
  out <- ode(y = init, # initial conditions
             times = dat$distance * 1000, # distance (converted from km to m)
             func = model, # the system of odes
             parms = parameters) # parameter values

  fit <- data.frame(out)
  
  fit <- fit %>%
    mutate(n15 = ((n15g * Ng) + (n15b * Nb))/(Ng + Nb))
  
  # return: RMSE with respect to delN-15
  sum((dat$d15N - fit$n15)^2)
  
}


# optimize with some sensible conditions
Opt <- optim(c(0.2, 1),          # initial guess
             RSS,                  # objective function
             method = "L-BFGS-B",  # optimization routine
             lower = c(0, 0),      # lower bound of search range
             upper = c(1, 1))      # upper bound of search range

Opt_par <- setNames(Opt$par, c("phi", "gamma")) # rename the optimal parameter set


######### VISUALIZATION ######### 
# run the model with the optimal parameters
fit <- data.frame(ode(y = init, times = dist, func = model, parms = Opt_par))

fit <- fit %>%
  mutate(n15 = ((n15g * Ng) + (n15b * Nb))/(Ng + Nb))

p <- ggplot(fit, aes(time, n15)) +
  geom_line() +
  geom_point(data = dat, aes(x = distance*1000, y = d15N)) +
  labs(x = "Distance (meters)")

print(p)

######## SENSITIVITY ANALYSIS ######## 
# construct a monte-carlo sensitivity test
# itterate through random parameter conbinations - retain objective function output
# plot objective function by parameter value, for each parameter. 

