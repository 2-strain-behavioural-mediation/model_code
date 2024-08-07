# BEHAVIOURAL MEDIATION SCENARIO ANALYSIS: RANKING DIFFERENT BEHAVIOURAL MECHANISMS
rm(list=ls())
# setwd()
library(dplyr)
library(deSolve)
library(ggplot2)
library(tidyr)

# ELO RATING FUNCTION
# K FACTOR SET EQUAL TO 32 AS PER CHESS GAME CONVENTION
elo_update <- function(ER1, ER2, ES1, EK = 32) {
  if (is.na(ER1) || is.na(ER2) || is.na(ES1)) {
    return(c(ER1, ER2))
  }
  EP1 <- 1 / (1 + 10^((ER2 - ER1) / 400))
  ER1_new <- ER1 + EK * (ES1 - EP1)
  ER2_new <- ER2 + EK * ((1 - ES1) - (1 - EP1))
  return(c(ER1_new, ER2_new))
}

# BEHAVIOURAL MODEL 1: TOTAL NUMBER INFECTED
SIR_2strains_E1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N1 <- S1 + I1 + R1
    N2 <- S2 + I2 + R2
    
    pre_peak <- 0
    post_peak <- 1
    
    lam1 <- ifelse(t < start_kappa, beta_1 * I1 / ((1 + pre_peak * kappa * (I1 + I2)) * N1), 
                   beta_1 * I1 / ((1 + post_peak * kappa * (I1 + I2)) * N1))
    lam2 <- ifelse(t < start_kappa, beta_2 * I2 / ((1 + pre_peak * kappa * (I1 + I2)) * N2), 
                   beta_2 * I2 / ((1 + post_peak * kappa * (I1 + I2)) * N2))
    
    dS1 <- mu*N1 + omega1*R1 - lam1*S1 - mu*S1
    dS2 <- mu*N2 + omega2*R2 - lam2*S2 - mu*S2
    dI1 <- lam1*S1 - nu1*I1 - mu*I1 - alpha*I1
    dR1 <- nu1*I1 - mu*R1 - omega1*R1
    dD1 <- alpha*I1
    dI2 <- lam2*S2 - nu2*I2 - mu*I2 - alpha*I2
    dR2 <- nu2*I2 - mu*R2 - omega2*R2
    dD2 <- alpha*I2
    
    list(c(dS1, dS2, dI1, dR1, dD1, dI2, dR2, dD2, growth_rate))
  })
}

# BEHAVIOURAL MODEL 2: TOTAL DISEASE INDUCED DEATHS
SIR_2strains_E2 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N1 <- S1 + I1 + R1
    N2 <- S2 + I2 + R2
    
    pre_peak <- 0
    post_peak <- 1
    
    lam1 <- ifelse(t < start_kappa, beta_1 * I1 / ((1 + pre_peak * alpha * kappa * (I1 + I2)) * N1), beta_1 * I1 / ((1 + post_peak * alpha * kappa * (I1 + I2)) * N1))
    lam2 <- ifelse(t < start_kappa, beta_2 * I2 / ((1 + pre_peak * alpha * kappa * (I1 + I2)) * N2), beta_2 * I2 / ((1 + post_peak * alpha * kappa * (I1 + I2)) * N2))
    
    dS1 <- mu*N1 + omega1*R1 - lam1*S1 - mu*S1
    dS2 <- mu*N2 + omega2*R2 - lam2*S2 - mu*S2
    dI1 <- lam1*S1 - nu1*I1 - mu*I1 - alpha*I1
    dR1 <- nu1*I1 - mu*R1 - omega1*R1
    dD1 <- alpha*I1
    dI2 <- lam2*S2 - nu2*I2 - mu*I2 - alpha*I2
    dR2 <- nu2*I2 - mu*R2 - omega2*R2
    dD2 <- alpha*I2
    
    list(c(dS1, dS2, dI1, dR1, dD1, dI2, dR2, dD2, growth_rate))
  })
}

# BEHAVIOURAL MODEL 3: PER CAPITA GROWTH
SIR_2strains_E3 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N1 <- S1 + I1 + R1
    N2 <- S2 + I2 + R2
    
    pre_peak <- 0
    post_peak <- 1
    
    lam1 <- ifelse(t < start_kappa, beta_1 * I1 / ((1 + pre_peak * kappa * growth_rate) * N1), beta_1 * I1 / ((1 + post_peak * kappa * growth_rate) * N1))
    lam2 <- ifelse(t < start_kappa, beta_2 * I2 / ((1 + pre_peak * kappa * growth_rate) * N2), beta_2 * I2 / ((1 + post_peak * kappa * growth_rate) * N2))
    
    dS1 <- mu*N1 + omega1*R1 - lam1*S1 - mu*S1
    dS2 <- mu*N2 + omega2*R2 - lam2*S2 - mu*S2
    dI1 <- lam1*S1 - nu1*I1 - mu*I1 - alpha*I1
    dR1 <- nu1*I1 - mu*R1 - omega1*R1
    dD1 <- alpha*I1
    dI2 <- lam2*S2 - nu2*I2 - mu*I2 - alpha*I2
    dR2 <- nu2*I2 - mu*R2 - omega2*R2
    dD2 <- alpha*I2
    
    growth_rate <- (dI1 + dI2)/(I1 + I2)
    
    list(c(dS1, dS2, dI1, dR1, dD1, dI2, dR2, dD2, growth_rate))
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
            I2 = initI2, R2 = initR2, D2 = initD2, growth_rate = 0)

parameters_E <- c(
  mu = (1/(50*52*7)),
  beta_1 = 0.5,
  beta_2 = 0.8,
  alpha = (168/100000),
  nu1 = 1/10,
  nu2 = 1/10,
  omega1 = 0,
  omega2 = 0,
  kappa = 0.1,
  start_kappa = 16.78,
  start_I2 = 10.0
)

# SIMULATION TIME SETTINGS
time_start <- 0
time_stop <- 500
deltat <- 0.1
tps <- seq(time_start, time_stop, by = deltat)


# INITIALIZING I2 = 1 WHEN T = START_I2
I2_init <- data.frame(var = "I2", time = parameters_E["start_I2"], value = 1, method = "replace")

# RUNNING SIMULATIONS FOR EACH BEHAVIOURAL MODEL
out_E1 <- as.data.frame(ode(y = istate, times = tps, func = SIR_2strains_E1, 
                            parms = parameters_E, events = list(data = I2_init),
                            rtol=1e-6, atol=1e-6, hmax=1/120))

out_E2 <- as.data.frame(ode(y = istate, times = tps, func = SIR_2strains_E2, 
                            parms = parameters_E, events = list(data = I2_init),
                            rtol=1e-6, atol=1e-6, hmax=1/120))

out_E3 <- as.data.frame(ode(y = istate, times = tps, func = SIR_2strains_E3, 
                            parms = parameters_E, events = list(data = I2_init),
                            rtol=1e-6, atol=1e-6, hmax=1/120))

print(str(out_E3))
print(head(out_E3))
print(tail(out_E3))

# CALCULATING RT_1 AND RT_2 VALUES FOR OUT_E1
# ESTIMATING RT_1
out_E1 <- out_E1 %>%
  mutate(Rt_1 = (parameters_E["beta_1"] * out_E1[[2]])/
           ((out_E1[[2]] + out_E1[[4]] + out_E1[[5]])*
              (1 + parameters_E["kappa"] * (out_E1[[4]] + out_E1[[7]]))*
              (parameters_E["nu1"] + parameters_E["mu"] + parameters_E["alpha"])))
# ESTIMATING RT_2
out_E1 <- out_E1 %>%
  mutate(Rt_2 = (parameters_E["beta_2"]*out_E1[["S2"]])/
           ((out_E1[["S2"]] + out_E1[["I2"]] + out_E1[["R2"]])*
              (1 + parameters_E["kappa"] * (out_E1[["I1"]] + out_E1[["I2"]]))*
              (parameters_E["nu2"] + parameters_E["mu"] + parameters_E["alpha"])))

# CALCULATING RT_1 AND RT_2 VALUES FOR OUT_E2
# ESTIMATING RT_1
out_E2 <- out_E2 %>%
  mutate(Rt_1 = (parameters_E["beta_1"] * out_E2[[2]])/
           ((out_E2[[2]] + out_E2[[4]] + out_E2[[5]])*
        (1 + parameters_E["alpha"] * parameters_E["kappa"] * (out_E2[[4]] + out_E2[[7]]))*
              (parameters_E["nu1"] + parameters_E["mu"] + parameters_E["alpha"])))
# ESTIMATING RT_2
out_E2 <- out_E2 %>%
  mutate(Rt_2 = (parameters_E["beta_2"]*out_E2[["S2"]])/
           ((out_E2[["S2"]] + out_E2[["I2"]] + out_E2[["R2"]])*
        (1 + parameters_E["alpha"] * parameters_E["kappa"] * (out_E2[["I1"]] + out_E2[["I2"]]))*
              (parameters_E["nu2"] + parameters_E["mu"] + parameters_E["alpha"])))

# CALCULATING RT_1 AND RT_2 VALUES FOR OUT_E3
# ESTIMATING RT_1
out_E3 <- out_E3 %>%
  mutate(Rt_1 = (parameters_E["beta_1"] * out_E3[[2]])/
           ((out_E3[[2]] + out_E3[[4]] + out_E3[[5]])*
              (1 + parameters_E["kappa"] * (out_E3[["growth_rate"]])) *
              (parameters_E["nu1"] + parameters_E["mu"] + parameters_E["alpha"])))
# ESTIMATING RT_2
out_E3 <- out_E3 %>%
  mutate(Rt_2 = (parameters_E["beta_2"]*out_E3[["S2"]])/
           ((out_E3[["S2"]] + out_E3[["I2"]] + out_E3[["R2"]])*
              (1 + parameters_E["kappa"] * (out_E3[["growth_rate"]])) *
              (parameters_E["nu2"] + parameters_E["mu"] + parameters_E["alpha"])))

# RELATIVE FITNESS RT_1/RT_2 AND ELO RATINGS CALCULATIONS
# STARTING AT 1500 AS IS CONVENTION
calculate_elo <- function(E_df) {
  elo_E1 <- 1500
  elo_E2 <- 1500
  elo_E3 <- 1500
  
  elo_ratings <- data.frame(time = E_df$time, E1 = NA, E2 = NA, E3 = NA)
  elo_ratings$E1[1] <- 1500
  elo_ratings$E2[1] <- 1500
  elo_ratings$E3[1] <- 1500
  
  # DEFINING RELATIVE FITNESS AS RT_1/RT_2 COMPARED WITH POTENTIAL FITNESS BETA_1/BETA_2
  for (i in 2:nrow(E_df)) {
    rf_E1 <- out_E1$Rt_1[i] / out_E1$Rt_2[i]
    rf_E2 <- out_E2$Rt_1[i] / out_E2$Rt_2[i]
    rf_E3 <- out_E3$Rt_1[i] / out_E3$Rt_2[i]
    
  # PAIRWISE COMPARISON OF RELATIVE FITNESS ACROSS BEHAVIOURAL MODELS  
    elo_E1 <- elo_update(elo_E1, elo_E2, ifelse(rf_E1 < rf_E2, 1, ifelse(rf_E1 > rf_E2, 0, 0.5)))[1]
    elo_E1 <- elo_update(elo_E1, elo_E3, ifelse(rf_E1 < rf_E3, 1, ifelse(rf_E1 > rf_E3, 0, 0.5)))[1]
    
    elo_E2 <- elo_update(elo_E2, elo_E1, ifelse(rf_E2 < rf_E1, 1, ifelse(rf_E2 > rf_E1, 0, 0.5)))[1]
    elo_E2 <- elo_update(elo_E2, elo_E3, ifelse(rf_E2 < rf_E3, 1, ifelse(rf_E2 > rf_E3, 0, 0.5)))[1]
    
    elo_E3 <- elo_update(elo_E3, elo_E1, ifelse(rf_E3 < rf_E1, 1, ifelse(rf_E3 > rf_E1, 0, 0.5)))[1]
    elo_E3 <- elo_update(elo_E3, elo_E2, ifelse(rf_E3 < rf_E2, 1, ifelse(rf_E3 > rf_E2, 0, 0.5)))[1]
    
    elo_ratings$E1[i] <- elo_E1
    elo_ratings$E2[i] <- elo_E2
    elo_ratings$E3[i] <- elo_E3
    
    # ADDING DEBUG PRINT STATEMENTS
    if (i %% 1000 == 0) {
      print(paste("elo_E1:", elo_E1, "elo_E2:", elo_E2, "elo_E3:", elo_E3))
    }
    
  }
  
  return(elo_ratings)
}

elo_ratings <- calculate_elo(out_E1)

# USING PRINT STATEMENTS FOR TROUBLESHOOTING
print(summary(elo_ratings))
print(tail(elo_ratings))

# PLOTTING ELO RATINGS
ggplot(elo_ratings, aes(x = time)) +
  geom_line(aes(y = E1, color = "E1: Total number infected")) +
  geom_line(aes(y = E2, color = "E2: Total disease induced deaths")) +
  geom_line(aes(y = E3, color = "E3: Per capita growth")) +
  labs(title = "Variant 2 Relative Fitness Scores Across Behavioural Models",
       subtitle = "Higher Elo scores indicate greater relative fitness of variant 2",
       x = "Time", y = "Variant 2 Relative Fitness Score") +
  scale_color_manual(values = c("E1: Total number infected" = "red",
                                "E2: Total disease induced deaths" = "blue",
                                "E3: Per capita growth" = "green")) +
  theme_minimal() +
  theme(legend.title = element_blank())

ggsave("variant_2_advantage_scores.png", width = 10, height = 6)
