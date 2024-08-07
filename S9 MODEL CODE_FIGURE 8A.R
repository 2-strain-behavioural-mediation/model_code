
rm(list=ls())
# setwd()
library(deSolve)
library(ggplot2)
library(reshape2)
library(viridis)
library(viridisLite)
library(dplyr)
library(tidyverse)
library(tidyr)
library(scales)
library(splines)
library(operators)
library(magrittr)
library(readr)

SIR_2strains_1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N1 <- (S1 + I1 + R1)
    N2 <- (S2 + I2 + R2)
    # USING BASE KAPPA VALUE IF KAPPA NOT ADJUSTED FOR INFECTION NUMBERS ABOVE TARGET
    current_kappa <- ifelse(policy_kappa_threshold == 0, kappa, policy_kappa_threshold)
    # USING THIS KAPPA VALUE IN OVERALL FUNCTION ONLY ON OR AFTER START OF BEHAVIOURAL MEDIATION
    effective_kappa <- ifelse(t >= start_kappa, max(current_kappa, kappa), 0)
    
    lam1 <- beta_1 * I1 / ((1 + effective_kappa * (I1 + I2)) * N1)
    lam2 <- beta_2 * I2 / ((1 + effective_kappa * (I1 + I2)) * N2)
    
    dS1 <- mu*N1 + omega1*R1 - lam1*S1 - mu*S1
    dS2 <- mu*N2 + omega2*R2 - lam2*S2 - mu*S2
    dI1 <- lam1*S1 - nu1*I1 - mu*I1 - alpha*I1
    dR1 <- nu1*I1 - mu*R1 - omega1*R1
    dD1 <- alpha*I1
    dI2 <- lam2*S2 - nu2*I2 - mu*I2 - alpha*I2
    dR2 <- nu2*I2 - mu*R2 - omega2*R2
    dD2 <- alpha*I2
    
    dN1 <- dS1 + dI1 + dR1
    dN2 <- dS2 + dI2 + dR2
    
    dpol_inc_1 <- dI1
    dpol_inc_2 <- dI2
    dpol_inc_tot <- dI1 + dI2 - 
      ((dI1 * I2 + I1 * dI2) * (N1 * N2) - (I1 * I2) * (dN1 * N2 + N1 * dN2)) / (N1 * N2)^2
    # POLICY KAPPA THRESHOLD STEPWISE FUNCTION: INCREASING KAPPA WHEN INC_TOT ABOVE TARGET INC_HAT
    dpolicy_kappa_threshold <- ifelse(pol_inc_tot > pol_inc_hat, 
                                      (1 - pol_inc_hat / pol_inc_tot) - policy_kappa_threshold, 
                                      -policy_kappa_threshold)
    
    list(c(dS1, dS2, dI1, dR1, dD1, dI2, dR2, dD2, 
           dpol_inc_1, dpol_inc_2, dpol_inc_tot, dpolicy_kappa_threshold))
  })
}

# INITIAL CONDITIONS AND PARAMETERS
initN1 <- 999
initN2 <- 999
initI1 <- 1
initR1 <- 0
initD1 <- 0
initI2 <- 0
initR2 <- 0
initD2 <- 0
initS1 <- initN1 - initI1 - initR1
initS2 <- initN2 - initI2 - initR2
istate <- c(S1 = initS1, S2 = initS2, I1 = initI1, R1 = initR1, D1 = initD1,
            I2 = initI2, R2 = initR2, D2 = initD2,
            pol_inc_1 = 0, pol_inc_2 = 0, pol_inc_tot = 0, policy_kappa_threshold = 0)

parameters_1 <- c(
  mu = (1/(50*52*7)), # background death rate
  beta_1 = 0.5, # underlying rate of transmission for variant 1 = p * c
  beta_2 = 0.8, # underlying rate of transmission for variant 2 = p * c
  alpha = (167.7/(100000*52*7)), # disease induced death rate gov.uk = 167.7 per 100000 annualised
  nu1 = 1/10, # rate of recovery = 1/(avg duration of infection)
  nu2 = 1/10, # rate of recovery = 1/(avg duration of infection)
  omega1 = 0, # loss of immunity rate = 1/(avg duration of immunity) [assume 0]
  omega2 = 0, # loss of immunity rate = 1/(avg duration of immunity) [assume 0]
  kappa = 0.1, # reduction in transmission due to behaviour change
  start_kappa = 10, # START OF BEHAVIOURAL MEDIATION
  start_I2 = 5, # START OF VARIANT 2 TRANSMISSION
  pol_inc_hat = 20 # target transmission rate R(t) for population
)

# SETTING SIMULATION TIME
time_start <- 0
time_stop <- 500
deltat <- 0.1
tps <- seq(time_start, time_stop, by = deltat)

# INITIALISING I2 = 1 WHEN T = START_I2
I2_init <- data.frame(var = "I2", time = parameters_1["start_I2"], value = 1, method = "replace")

# RUNNING ODE AND STORING IN OUT_1 DATAFRAME
out_1 <- as.data.frame(deSolve::ode(y = istate, times = tps, func = SIR_2strains_1,
                                    parms = parameters_1, events = list(data = I2_init),
                                    rtol = 1e-12, hmax = 1/120))

out_1 <- out_1 %>%
  mutate(pol_inc_hat = parameters_1["pol_inc_hat"])

summary(out_1)


# SAVING TO CSV
write_csv(out_1, "SIR_2strains_threshold_240730.csv")
print("Data has been saved to SIR_2strains_threshold_240730.csv")

# PLOTTING INCIDENCE: POL_INC_1, POL_INC_2, POL_INC_HAT, POL_INC_TOT, AND POLICY_KAPPA_THRESHOLD V TIME
# RESHAPING DATA FROM WIDE TO LONG FORMAT
# CONVERTING DATA FRAME TO TIBBLE FOR PLOTTING
out_long <- out_1 %>%
  dplyr::select(time, pol_inc_1, pol_inc_2, pol_inc_tot, pol_inc_hat, policy_kappa_threshold) %>%
  pivot_longer(cols = -time, names_to = "variable", values_to = "value")

# PLOTTING POL_INC_TOT, POL_INC_HAT ON LEFT AXIS AND POLICY_KAPPA_THRESHOLD ON RIGHT AXIS
ggplot(out_long, aes(x = time, y = value, color = variable)) +
  geom_line(data = . %>% filter(variable != "policy_kappa_threshold")) +
  scale_color_manual(values = c("pol_inc_1" = "blue", 
                                "pol_inc_2" = "red", 
                                "pol_inc_tot" = "black", 
                                "pol_inc_hat" = "orange", 
                                "policy_kappa_threshold" = "green"),
                     name = "Legend") +
  scale_y_continuous(name = "Infected Individuals",
                     sec.axis = sec_axis(~ . / 100, name = "Policy Kappa Threshold")) +
  labs(title = "Policy Kappa Targeting Threshold Infection Numbers",
       x = "Time") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),  
    axis.text = element_text(size = 16),  
    plot.title = element_text(size = 20), 
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  geom_line(data = filter(out_long, variable == "policy_kappa_threshold"), 
            aes(y = value * 100), linetype = "dashed", linewidth = 1.5)

# SAVING THE PLOT
ggsave("Policy_threshold_plot.png", width = 20, height = 6)
