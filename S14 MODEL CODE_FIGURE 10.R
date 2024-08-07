
# ASSUMING NO CROSS IMMUNITY
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

# BEHAVIOURAL MODEL TARGETING TOTAL NUMBER INFECTED
SIR_2strains_1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    
    N1 <- (S1 + I1 + R1)
    N2 <- (S2 + I2 + R2)
    
    # USING KAPPA ONLY AT OR AFTER START_KAPPA
    effective_kappa_1 <- ifelse(t >= start_kappa, kappa, 0)
    
    lam1 <- beta_1 * I1 / ((1 + effective_kappa_1 * (I1 + I2)) * N1)
    lam2 <- beta_2 * I2 / ((1 + effective_kappa_1 * (I1 + I2)) * N2)
      
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
  alpha = (167.7/(100000*52*7)),  # disease induced death rate gov.uk = (167.7/100000)
  nu1 = 1/10, # rate of recovery = 1/(avg duration of infection)
  nu2 = 1/10, # rate of recovery = 1/(avg duration of infection)
  omega1 = 0, # loss of immunity rate = 1/(avg duration of immunity) [assume 0]
  omega2 = 0, # loss of immunity rate = 1/(avg duration of immunity) [assume 0]
  kappa = 0.1, # reduction in transmission due to behaviour change
  start_kappa = 16.78, # START OF BEHAVIOURAL MEDIATION
  start_I2 = 10 # START OF VARIANT 2 TRANSMISSION
)

# SIMULATION
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


# TYPE REPRODUCTION NUMBER
# MATRIX F=NEW INFECTIONS/TRANSMISSION, MATRIX V=INTER-COMPARTMENT TRANSITION
# SOURCES: DIECKMANN ET AL (2012); BJOERNSTAD, O.N (2018)
# ALL NEW INFECTIONS
F1 = quote((beta_1 * S1 * I1) / ((S1 + I1 + R1) * (1 + kappa * (I1 + I2))))
F2 = quote((beta_2 * S2 * I2) / ((S2 + I2 + R2) * (1 + kappa * (I1 + I2))))
# ALL LOSSES (TRANSMISSION + TRANSITION)
Vm1 = quote((nu1 + mu + alpha) * I1)
Vm2 = quote((nu2 + mu + alpha) * I2)
# ALL TRANSMISSION GAINS
Vp1 = 0
Vp2 = 0
# SUBTRACTING VP FROM VM
V1 = substitute(a - b, list(a = Vm1, b = Vp1))
V2 = substitute(a - b, list(a = Vm2, b = Vp2))
# CREATING PARTIAL DERIVATIVES TO CALCULATE JACOBIANS
f11 = D(F1, "I1"); f12 = D(F1, "I2")
f21 = D(F2, "I1"); f22 = D(F2, "I2")
v11 = D(V1, "I1"); v12 = D(V1, "I2")
v21 = D(V2, "I1"); v22 = D(V2, "I2")
# INCORPORATING PARAMETER VALUES
# paras=list(S1=1, S2=1, I1=0, R1=0, I2=0, R2=0, parameters_1)
# CREATING MATRICES FOR EIGENVECTOR OPERATIONS
f=with(as.list(c(istate, parameters_1)),
       matrix(c(eval(f11),eval(f12),eval(f21),
                eval(f22)), nrow=2, byrow=TRUE))
v=with(as.list(c(istate, parameters_1)),
       matrix(c(eval(v11),eval(v12),eval(v21),
                eval(v22)), nrow=2, byrow=TRUE))
# CALCULATING DOMINANT EIGENVALUE # if just R0 - max(eigen(f %*% solve(v))$values)
K_mat <- f %*% solve(v)
K_mat
# USING NEXT GENERATION MATRIX K TO FIND TYPE REPRODUCTION NUMBER
# K_MAT IS NEXT GENERATION MATRIX CALCULATED ABOVE
# I_MAT IS IDENTITY MATRIX IN R^N (WHERE R^N = N-FOLD CARTESTIAN PRODUCT)
# P_MAT IS N X N PROJECTION MATRIX OF SUBSPACE COVERED BY FIRST L UNIT VECTORS IN R^N 
# E_MAT IS N X L MATRIX OF FIRST L UNIT VECTORS ALIGNED IN R^N
# M_MAT IS L X L MATRIX DEFINED BELOW
# T_R_N IS SPECTRAL RADIUS OR LARGEST ABSOLUTE VALUES (COMPLEX MODULUS) OF EIGENVALUES OF M_MAT
l <- 2  # NUMBER OF TYPES TO TARGET WITH CONTROL (1<=L<=N)
n <- 2  # NUMBER OF TOTAL TYPES
I_mat <- diag(n) # CREATING N X N IDENTITY MATRIX
U_mat <- matrix(c(rep(1,l), rep(0,n-l)), nrow=n, ncol=l) # BUILDING PROJECTION MATRIX
P_mat <- U_mat %*% t(U_mat) # PROJECTION MATRIX
E_mat <- matrix(c(rep(1,l), rep(0,n*(l-1))), nrow=n, ncol=l)
# CALCULATING THE L X L MATRIX M_MAT AND THE TYPE REPRODUCTION NUMBER T (T_R_N = T)
M_mat <- t(E_mat) %*% P_mat %*% K_mat %*% solve(I_mat-(I_mat-P_mat) %*% K_mat) %*% E_mat
T_R_N <- max(Re(eigen(M_mat)$values))
T_R_N

# CREATING WHOLE POPULATION NUMBER
out_1 <- out_1 %>%
  mutate(S_prop = ((out_1[["S1"]]+out_1[["S2"]])/2)/
           (((out_1[[2]] + out_1[[4]] + out_1[[5]])+
               (out_1[[3]] + out_1[[7]] + out_1[[8]]))/2))
# MULTIPLYING TYPE REPRODUCTION NUMBE BY S(T) (WHERE R(T) = R0 * PROP. SUSCEPTIBLE)
out_1 <- out_1 %>%
  mutate(Typ_rep_num = T_R_N * S_prop)

out_1
summary(out_1)

# PLOTTING VARIANT 1, VARIANT 2, AND TYP_REP_NUM V TIME
# RESHAPING DATA FROM WIDE TO LONG FORMAT
# CONVERTING DATA FRAME TO TIBBLE FOR PLOTTING
typ_rep_num_max <- max(out_1$Typ_rep_num)
out_long <- out_1 %>%
  dplyr::select(time, I1, I2, Typ_rep_num) %>%
  pivot_longer(cols = -time, names_to = "variable", values_to = "value")

# RENAMING VARIABLES
out_long$variable <- factor(out_long$variable, 
                            levels = c("I1", "I2", "Typ_rep_num"),
                            labels = c("Variant 1", "Variant 2", "Type Reproduction Number"))

# PLOTTING I1 AND I2 ON LEFT AXIS AND TYP_REP_NUM ON RIGHT AXIS
ggplot(out_long, aes(x = time, y = value, color = variable)) +
  geom_line(data = . %>% filter(variable != "Type Reproduction Number"), linewidth = 1.2) +
  scale_color_manual(values = c("Variant 1" = "blue", "Variant 2" = "red", "Type Reproduction Number" = "green"),
                     name = "Legend") +
  scale_y_continuous(name = "Variant 1, Variant 2 Prevalence",
                     sec.axis = sec_axis(~ . / 100, name = "Type Reproduction Number", breaks = seq(0, 5, by = 1))) +
  labs(x = "Time", title = "Variant Prevalence and Type Reproduction Number") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20)
  ) +
  geom_line(data = filter(out_long, variable == "Type Reproduction Number"), 
            aes(y = value * 100), linetype = "dashed", linewidth = 1.2)
