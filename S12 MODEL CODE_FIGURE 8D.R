
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

# USING SIR_2STRAINS_1 FUNCTION TO MODEL POLICY KAPPA TOOL TARGETING REPRODUCTION NUMBERS
SIR_2strains_1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N1 <- S1 + I1 + R1
    N2 <- S2 + I2 + R2
    # DEFINING RT FOR VARIANT 1 AND VARIANT 2
    rt_1_base <- (beta_1 * S1) / (N1 * (nu1 + mu + alpha))
    rt_2_base <- (beta_2 * S2) / (N2 * (nu2 + mu + alpha))
    # USING PLACEHOLDER TO SIMPLIFY
    w1 <- I1 / N1
    w2 <- I2 / N2
    # SOLVING FOR RT FOR TOTAL POPULATION
    rt_tot_base <- (rt_1_base * w1 + rt_2_base * w2 - w1 * w2 * (rt_1_base + rt_2_base)) / 
      (w1 + w2 - w1 * w2)
    # POLICY KAPPA RT STEPWISE FUNCTION: INCREASING KAPPA WHEN RT ABOVE TARGET RT_HAT
    if (rt_tot_base > rt_hat) {
      kappa_solution <- (rt_tot_base - rt_hat) / ((I1 + I2) * (rt_tot_base - rt_hat + rt_hat/rt_tot_base))
    } else {
      kappa_solution <- 0
    }
    # USING THIS KAPPA VALUE IN OVERALL FUNCTION ONLY ON OR AFTER START OF BEHAVIOURAL MEDIATION
    effective_kappa <- ifelse(t >= start_kappa, max(kappa_solution, kappa), 0)
    # ADDING BACK THIS TIMING ADJUSTED KAPPA VALUE TO FINAL RT EQUATIONS FOR EACH POPULATION
    rt_1 <- rt_1_base / (1 + effective_kappa * (I1 + I2))
    rt_2 <- rt_2_base / (1 + effective_kappa * (I1 + I2))
    rt_tot <- (rt_1 * w1 + rt_2 * w2 - w1 * w2 * (rt_1 + rt_2)) / (w1 + w2 - w1 * w2)
    
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
    
    list(c(dS1, dS2, dI1, dR1, dD1, dI2, dR2, dD2),
         c(rt_1 = rt_1, rt_2 = rt_2, rt_tot = rt_tot, kappa_solution = kappa_solution))
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
            I2 = initI2, R2 = initR2, D2 = initD2)

parameters_1 <- c(
  mu = 1/(50*52*7), # background death rate
  beta_1 = 0.5, # underlying rate of transmission for variant 1 = p * c
  beta_2 = 0.8, # underlying rate of transmission for variant 2 = p * c
  alpha = (167.7/(100000*52*7)), # disease induced death rate gov.uk = 167.7 per 100000 annualised
  nu1 = 1/10, # rate of recovery = 1/(avg duration of infection)
  nu2 = 1/10, # rate of recovery = 1/(avg duration of infection)
  omega1 = 0, # loss of immunity rate = 1/(avg duration of immunity) [assume 0]
  omega2 = 0, # loss of immunity rate = 1/(avg duration of immunity) [assume 0]
  kappa = 0.1, # reduction in transmission due to behaviour change
  start_kappa = 30, # START OF BEHAVIOURAL MEDIATION
  start_I2 = 20, # START OF VARIANT 2 TRANSMISSION
  rt_hat = 1 # target transmission rate R(t) for population
)

# SETTINGS FOR SIMULATION TIME
time_start <- 0
time_stop <- 500
deltat <- 0.1
tps <- seq(time_start, time_stop, by = deltat)

# INITIALISING I2=1 WHEN T=START_KAPPA
I2_init <- data.frame(var = "I2", time = parameters_1["start_I2"], value = 1, method = "replace")

# RUNNING ODE
out_1 <- as.data.frame(ode(y = istate, times = tps, func = SIR_2strains_1,
                           parms = parameters_1, events = list(data = I2_init),
                           method = "lsoda", atol = 1e-12, rtol = 1e-12))

out_1$rt_hat <- parameters_1["rt_hat"]

# SAVING TO CSV
write_csv(out_1, "SIR_2strains_threshold_240731.csv")
print("Data has been saved to SIR_2strains_threshold_240731.csv")

# PRINTING SUMMARY
summary(out_1)

# RESHAPING DATA FROM WIDE TO LONG FORMAT
# CONVERTING DATA FRAME TO TIBBLE FOR PLOTTING
out_long <- out_1 %>%
  dplyr::select(time, rt_1, rt_2, rt_tot, rt_hat, kappa_solution) %>%
  pivot_longer(cols = -time, names_to = "variable", values_to = "value")

# RENAMING KAPPA_SOLUTION TO POLICY_KAPPA_RT FOR LEGEND
out_long$variable <- factor(out_long$variable, 
                            levels = c("rt_1", "rt_2", "rt_tot", "rt_hat", "kappa_solution"),
                            labels = c("rt_1", "rt_2", "rt_tot", "rt_hat", "policy_kappa_rt"))

# PLOTTING RT_1, RT_2, RT_TOT, RT_HAT ON LEFT AXIS AND KAPPA_SOLUTION ON RIGHT AXIS
ggplot(out_long, aes(x = time, y = value, color = variable)) +
  geom_line(data = . %>% filter(variable != "policy_kappa_rt" & variable != "rt_tot")) +
  geom_line(data = . %>% filter(variable == "rt_tot"), linewidth = 1.5) +
  scale_color_manual(values = c("rt_1" = "blue", "rt_2" = "red", "rt_tot" = "black", 
                                "rt_hat" = "orange", "policy_kappa_rt" = "green"),
                     name = "Legend") +
  scale_y_continuous(name = "Reproduction Number R(t)",
                     sec.axis = sec_axis(~ . / max(out_long$value[out_long$variable != "policy_kappa_rt"]), 
                                         name = "Policy Kappa R(t)",
                                         breaks = seq(0, 1, 0.2))) +
  labs(x = "Time", title = "Policy Kappa Targeting Reproduction Numbers") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  geom_line(data = filter(out_long, variable == "policy_kappa_rt"), 
            aes(y = value * max(out_long$value[out_long$variable != "policy_kappa_rt"])), 
            linetype = "dashed", linewidth = 1.5)

# SAVING THE PLOT
ggsave("policy_kappa_rt_2.png", width = 20, height = 6, dpi = 300)
print("Plot has been saved as policy_kappa_rt_2.png")