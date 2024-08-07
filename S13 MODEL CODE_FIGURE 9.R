
# QUANTIFYING POSITIVE FEEDBACK AND MODEL UNCERTAINTY
rm(list=ls())
# setwd()
library(dplyr)
library(deSolve)
library(ggplot2)
library(tidyr)
library(patchwork)

# POSITIVE FEEDBACK MEASUREMENT
# BEHAVIOURAL MODEL 1: TOTAL NUMBER INFECTED
SIR_2strains_feedback <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    N1 <- (S1 + I1 + R1)
    N2 <- (S2 + I2 + R2)
    
    pre_peak <- 0
    post_peak <- 1
    
    lam1 <- ifelse(t < start_kappa, 
                   beta_1 * I1 / ((1 + pre_peak * kappa * (I1 + I2)) * N1), 
                   beta_1 * I1 / ((1 + post_peak * kappa * (I1 + I2)) * N1))
    
    lam2 <- ifelse(t < start_kappa, 
                   beta_2 * I2 / ((1 + pre_peak * kappa * (I1 + I2)) * N2), 
                   beta_2 * I2 / ((1 + post_peak * kappa * (I1 + I2)) * N2))
    
    dS1 <- mu*N1 + omega1*R1 - lam1*S1 - mu*S1
    dS2 <- mu*N2 + omega2*R2 - lam2*S2 - mu*S2
    dI1 <- lam1*S1 - nu1*I1 - mu*I1 - alpha*I1
    dR1 <- nu1*I1 - mu*R1 - omega1*R1
    dD1 <- alpha*I1
    
    dI2 <- lam2*S2 - nu2*I2 - mu*I2 - alpha*I2
    dR2 <- nu2*I2 - mu*R2 - omega2*R2
    dD2 <- alpha*I2
    
    list(c(dS1, dS2, dI1, dR1, dD1, dI2, dR2, dD2))
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

parameters_f <- c(
  mu = (1/(50*52*7)), # background death rate
  beta_1 = 0.5, # underlying rate of transmission for variant 1 = p * c
  beta_2 = 0.8, # underlying rate of transmission for variant 2 = p * c
  alpha = (167.7/(100000*52*7)), # disease induced death rate gov.uk = 167.7 per 100000 annualised
  nu1 = 1/10, # rate of recovery = 1/(avg duration of infection)
  nu2 = 1/10, # rate of recovery = 1/(avg duration of infection)
  omega1 = 0, # loss of immunity rate = 1/(avg duration of immunity) [assume 0]
  omega2 = 0, # loss of immunity rate = 1/(avg duration of immunity) [assume 0]
  kappa = 0.1, # reduction in transmission due to behaviour change
  start_kappa = 16.78, # START OF BEHAVIOURAL MEDIATION
  start_I2 = 10 # START OF VARIANT 2 TRANSMISSION
)

# SIMULATION TIME SETTINGS
time_start <- 0
time_stop <- 500
deltat <- 0.1
tps <- seq(time_start, time_stop, by = deltat)

# INITIALIZING I2 = 1 WHEN T = START_I2
I2_init_f <- data.frame(var = "I2", time = parameters_f["start_I2"], value = 1, method = "replace")

# RUNNING SIMULATIONS FOR EACH BEHAVIOURAL MODEL
out_f <- as.data.frame(ode(y = istate, times = tps, func = SIR_2strains_feedback, 
                            parms = parameters_f, events = list(data = I2_init_f),
                            rtol=1e-6, atol=1e-6, hmax=1/120))

# CALCULATING CUMULATIVE BEHAVIOURAL EFFECT (CBE)
out_f$BE <- 1 / (1 + parameters_f["kappa"] * (out_f$I1 + out_f$I2))

# CALCULATING VARIANT INTERACTION INDEX (VII)
out_f$dI1 <- c(0, diff(out_f$I1) / diff(out_f$time))
out_f$dI2 <- c(0, diff(out_f$I2) / diff(out_f$time))
out_f$VI <- abs(out_f$dI1 * out_f$dI2)

# CREATING NEW COLUMN FOR POSITIVE_FEEDBACK AS THE PRODUCT OF BE AND VI
out_f$positive_feedback <- out_f$BE * out_f$VI

summary(out_f)


# RESHAPING DATA FOR GGPLOT
plot_data <- out_f %>%
  dplyr::select(time, I1, I2, positive_feedback) %>%
  gather(key = "variable", value = "value", -time)

# RENAMING VARIABLES
plot_data$variable <- factor(plot_data$variable, 
                             levels = c("I1", "I2", "positive_feedback"),
                             labels = c("Variant 1", "Variant 2", "Positive Feedback"))

# PLOTTING IN GGPLOT
ggplot(plot_data, aes(x = time, y = value, color = variable)) +
  geom_line(data = subset(plot_data, variable != "Positive Feedback"), linewidth = 1) +
  geom_line(data = subset(plot_data, variable == "Positive Feedback"), linewidth = 2) +
  scale_color_manual(values = c("Variant 1" = "blue", "Variant 2" = "red", "Positive Feedback" = "green")) +
  labs(title = "Variant 1 and Variant 2 Prevalence and Positive Feedback Effect",
       x = "Time",
       y = "Total Infections (Prevalence)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size = 16), 
    axis.text = element_text(size = 16),  
    plot.title = element_text(size = 20)  
  )

# SAVING PLOT
ggsave("positive_feedback.png", width = 20, height = 6)
