# APPARENT COMPETITION CHARTS USING SIR_2STRAINS_1
rm(list=ls())
# setwd()
library(deSolve)
library(numDeriv)
library(pracma)
library(Matrix)
library(ggplot2)
library(reshape2)
library(viridis)
library(viridisLite)
library(dplyr)
library(tidyverse)
library(tidyr)
library(scales)
library(rgl)
library(plot3D)
library(RColorBrewer)
library(scatterplot3d)
library(stringr)
library(data.table)
library(plotly)
library(raster)
library(purrr)
library(gridExtra)
library(patchwork)
library(splines)
library(operators)
library(magrittr)

# USING BEHAVIOURAL VARIABLE TARGETING TOTAL NUMBER INFECTED
SIR_2strains_1 <- function(t, state, parameters_1) {
  with(as.list(c(state, parameters_1)),{
    
    N1 <- (S1 + I1 + R1)
    N2 <- (S2 + I2 + R2)
    
    # USING KAPPA IN FORCE OF INFECTION LAMBDA TERMS ONLY ON OR AFTER START OF BEHAVIOURAL MEDIATION
    effective_kappa <- ifelse(t >= start_kappa, kappa, 0)
    
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

    list(c(dS1, dS2, dI1, dR1, dD1, dI2, dR2, dD2))
  })
}

# DEFINING INITIAL CONDITIONS AND PARAMETERS
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
  mu = (1/(50*52*7)), # background death rate
  beta_1 = 0.5, # underlying rate of transmission for variant 1 = p * c
  beta_2 = 0.8, # underlying rate of transmission for variant 2 = p * c
  alpha = (167.7/(100000*52*7)), # disease induced death rate 240503_gov.uk = (168/100000)
  nu1 = 1/10, # rate of recovery = 1/(avg duration of infection)
  nu2 = 1/10, # rate of recovery = 1/(avg duration of infection)
  omega1 = 0, # loss of immunity rate = 1/(avg duration of immunity) [assume 0]
  omega2 = 0, # loss of immunity rate = 1/(avg duration of immunity) [assume 0]
  kappa = 0.1, # reduction in transmission due to behaviour change
  start_kappa = 30, # START OF BEHAVIOURAL MEDIATION
  start_I2 = 21.2 # START OF VARIANT 2 TRANSMISSION
)

# SIMULATION TIME SETTINGS
time_start <- 0
time_stop <- 500
deltat<- 0.1
tps <- seq(time_start , time_stop , by = deltat)

# INITIALIZING I2 = 1 WHEN T = START_I2
I2_init <- data.frame(var = "I2", time = parameters_1["start_I2"], value = 1, method = "replace")

# RUNNING SYSTEM OF EQUATIONS
out_1 <- as.data.frame(deSolve::ode(y = istate, times = tps, func = SIR_2strains_1, 
                                    parms = parameters_1, events = list(data = I2_init),
                                    rtol=1e-12, hmax=1/120))

# FIND PEAK OF VARIANT 1
peak_time <- tps[which.max(out_1$I1)]
# PRINTING TIMING OF PEAK PREVALENCE OF VARIANT 1
peak_time
# CALCULATING MAX EXTINCTION RATE OF VARIANT 1
out_1 <- out_1 %>%
  mutate(logged_values_1 = log(out_1[["I1"]])) %>%
  mutate(logged_diff_1 = c(0, diff(logged_values_1))) %>%
  mutate(growth_rate_1 = 10 * (exp(logged_diff_1) - 1))
# PRINTING OUT MAXIMUM EXTINCTION RATE (MINIMUM GROWTH RATE) OF VARIANT 1
min(out_1[["growth_rate_1"]])
# PRINTING OUT PREAK PREVALENCE OF VARIANT 1
max(out_1[,4])
# PRINTING OUT PREAK PREVALENCE OF VARIANT 2
max(out_1[,7])

# RESHAPING FROM WIDE TO LONG FORMAT
out_long <- out_1 %>%
  dplyr::select(time, I1, I2) %>%
  pivot_longer(cols = -time, names_to = "variable", values_to = "value")

# USING GGPLOT TO CHART
ggplot(out_long, aes(x = time, y = value, color = variable)) +
  geom_line(linewidth = 1.2) +  # Increased line thickness
  geom_vline(aes(xintercept = parameters_1["start_kappa"], color = "Start of behavioural mediation"),
             linetype = "dashed", linewidth = 1.2) +  # Increased vertical line thickness
  scale_color_manual(values = c("I1" = "black", "I2" = "blue", "Start of behavioural mediation" = "orange"),
                     name = "Legend", labels = c("Variant 1", "Variant 2", "Start of behavioural mediation")) +
  scale_y_continuous(name = "Total Infections (Prevalence)", limits = c(0, 625)) +
  labs(x = "Time", title = "Variant 2 Introduced at Peak of Variant 1 with Behavioural Mediation Starting on Day 30") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "solid", "dashed"), linewidth = 1.2)))

# SAVING THE PLOT
ggsave("variant_prevalence_plot_3.png", width = 20, height = 6, dpi = 300)
print("Plot has been saved as variant_prevalence_plot_3.png")