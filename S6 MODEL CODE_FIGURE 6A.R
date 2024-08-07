
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
  kappa = 1.0, # reduction in transmission due to behaviour change
  start_kappa = 10, # START OF BEHAVIOURAL MEDIATION
  start_I2 = 10 # START OF VARIANT 2 TRANSMISSION
)

# Simulation
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

# finding time of peak of first variant: run after out_1 execution
peak_time <- tps[which.max(out_1$I1)]
# peak_time <- tps[which.max(out_1$incidence_1)]
peak_time

# plotting I1, I2
plot(out_1[,'I1']~tps, ylim = c(0, 600), xlab='Time')
lines (out_1[,'I2']~tps, col = "magenta", ylim = c(0, 600))
legend('topright', legend=c('Variant 1', 'Variant 2'), col=c('black','magenta'), lty=1)

out_1
# CALCULATING GROWTH RATE I1
out_1 <- out_1 %>%
  mutate(logged_values_1 = log(out_1[["I1"]]))
# CALCULATING DIFFERENCE BETWEEN LOGGED VALUES
out_1 <- out_1 %>%
  mutate(logged_diff_1 = c(0, diff(logged_values_1)))
# EXPONENTIATING LOGGED DIFFERENCES, SUBTRACING 1, X 10 FOR DAILY RATE

# CHECKING MIN GROWTH RATE 1 METRICS
out_1 <- out_1 %>%
  mutate(growth_rate_1 = 10 * (exp(logged_diff_1) - 1))
min(out_1[["growth_rate_1"]])

# New code to calculate and print min_FES_out
FES_1 <- out_1[[5]]  # D1 column
FES_2 <- out_1[[8]]  # D2 column
max_FES_1 <- max(FES_1, na.rm = TRUE)
max_FES_2 <- max(FES_2, na.rm = TRUE)
min_FES_out <- max_FES_1 / max_FES_2

# Print the result
cat("Minimum FES1/FES2 ratio (min_FES_out):", min_FES_out, "\n")

out_1 <- out_1 %>%
  mutate(FES_1 = out_1[[5]],
         FES_2 = out_1[[8]])

# CHECKING FES OUTPUT METRICS
out_1 <- out_1 %>%
  mutate(FES_tot = FES_1/FES_2)
max(out_1[["FES_1"]], na.rm = TRUE) / max(out_1[["FES_2"]], na.rm = TRUE)
min_FES_out = max(out_1[["FES_1"]], na.rm = TRUE) / max(FES_2, na.rm = TRUE)

summary_1_FES <- out_1 %>%
  summarise(max_FES_1 = max(FES_1, na.rm = TRUE),
            max_FES_2 = max(FES_2, na.rm = TRUE),
            min_FES_out_23 = max(FES_1, na.rm = TRUE) / max(FES_2, na.rm = TRUE),
            min_growth_rate_1 = min(growth_rate_1, na.rm = TRUE),
            .groups = 'drop')

# Print min_FES_out_23
cat("Minimum FES1/FES2 ratio (min_FES_out_23):", summary_1_FES$min_FES_out_23, "\n")


# GLOBAL SIMULATION SCRIPT
# MULTIVARIABLE SENSITIVITY ANALYSIS BETA_2, KAPPA, START_I2, START_KAPPA, STORING RESULTS
# SIMULATING START_KAPPA V START_I2
# CALCULATING METRICS, SENSITIVITY ANALYSIS OF KAPPA AND BETA_2, STORING RESULTS
beta_2_sensitivity_vector <-seq(0,1,by=0.1)
kappa_sensitivity_vector <-seq(0,1,by=0.1)
start_I2_sensitivity_vector <-seq(0,50,by=5)
start_kappa_sensitivity_vector <-seq(0,50,by=5)
result_global_1 <-matrix(0,nrow = 1,ncol = 25)
colnames(result_global_1)<-c("time","S1","S2","I1","R1","D1", "I2","R2", "D2","FES_1","FES_2","FES_out","logged_values_1","logged_diff_1","growth_rate_1","logged_values_2","logged_diff_2","growth_rate_2","Rt_1","Rt_2","growth_1_check","beta_2","kappa","start_I2","start_kappa")
result_global_1<-as.data.frame(result_global_1)
# CREATING OUTPUT FILE DIRECTORY
dir.create("simulation_results", showWarnings = FALSE)
# CREATING A COUNTER FOR FILE NAMING
file_counter <- 1
# LOOPING THROUGH SENSITIVITY VECTORS
for (i in 1:length(beta_2_sensitivity_vector)) {
  for (j in 1:length(kappa_sensitivity_vector)) {
    # CREATING A DATA FRAME TO STORE RESULTS FROM EACH I, J COMBINATION
    result_ij <- data.frame()
    
    
    for (q in 1:length(start_I2_sensitivity_vector)) {
      for (r in 1:length(start_kappa_sensitivity_vector)) {
        
    # CREATING A COPY OF PARAMETERS FOR EACH ITERATION
    parameters_global_1 <- parameters_1
    # UPDATING PARAMETERS FOR SIMULATIONS
    parameters_global_1["beta_2"]<-beta_2_sensitivity_vector[i]
    parameters_global_1["kappa"]<-kappa_sensitivity_vector[j]
    parameters_global_1["start_I2"]<-start_I2_sensitivity_vector[q]
    parameters_global_1["start_kappa"]<-start_kappa_sensitivity_vector[r]
    
    #ADDING PRINT STATEMENT TO SHOW PROGRESS
    cat("Running: i =", i, ", j =", j, ", q =", q, ", r =", r, "\n")
    
    I2_init_global_1 <- data.frame(var = "I2", time = parameters_global_1["start_I2"], value = 1, method = "replace")
    #FINDING NUMERICAL SOLUTION IN ODE SYSTEM
    out_global_1 <- as.data.frame(deSolve::ode(y = istate, times = tps, func = SIR_2strains_1, parms = parameters_global_1, events = list(data = I2_init_global_1), rtol=1e-12,hmax=1/120))
    # ADDING COLUMNS WITH MUTATE FUNCTION
    out_global_1 <- setNames(out_global_1, c("time","S1","S2","I1","RI","D1", "I2","R2", "D2"))
    # CALCULATING FINAL EPIDEMIC SIZE FES_1, FES_2
    out_global_1 <- out_global_1 %>%
      mutate(FES_1 = out_global_1[[5]],
             FES_2 = out_global_1[[8]])
    # CALCULATING FINAL EPIDEMIC SIZE FES_OUT
    out_global_1 <- out_global_1 %>%
      mutate(FES_out = FES_1/FES_2)
    # CALCULATING GROWTH RATE I1
    out_global_1 <- out_global_1 %>%
      mutate(logged_values_1 = log(out_global_1$I1))
    # CALCULATING DIFFERENCE BETWEEN LOGGED VALUES
    out_global_1 <- out_global_1 %>%
      mutate(logged_diff_1 = c(0, diff(logged_values_1)))
    # EXPONENTIATING LOGGED DIFFERENCES, SUBTRACING 1, X 10 FOR DAILY RATE
    out_global_1 <- out_global_1 %>%
      mutate(growth_rate_1 = 10 * (exp(logged_diff_1) - 1))
    # CALCULATING GROWTH RATE I2
    out_global_1 <- out_global_1 %>%
      mutate(logged_values_2 = log(out_global_1[["I2"]]))
    # CALCULATING DIFFERENCE BETWEEN LOGGED VALUES
    out_global_1 <- out_global_1 %>%
      mutate(logged_diff_2 = c(0, diff(logged_values_2)))
    # EXPONENTIATING LOGGED DIFFERENCES, SUBTRACING 1, X 10 FOR DAILY RATE
    out_global_1 <- out_global_1 %>%
      mutate(growth_rate_2 = 10 * (exp(logged_diff_2) - 1))
    # ESTIMATING RT_1
    out_global_1 <- out_global_1 %>%
      mutate(Rt_1 = (parameters_global_1["beta_1"] * out_global_1[[2]])/
        ((out_global_1[[2]] + out_global_1[[4]] + out_global_1[[5]])*
        (1 + parameters_global_1["kappa"] * (out_global_1[[4]] + out_global_1[[7]]))*
        (parameters_global_1["nu1"] + parameters_global_1["mu"] + parameters_global_1["alpha"])))
    # ESTIMATING RT_2
    out_global_1 <- out_global_1 %>%
      mutate(Rt_2 = (parameters_global_1["beta_2"]*out_global_1[["S2"]])/
        ((out_global_1[["S2"]] + out_global_1[["I2"]] + out_global_1[["R2"]])*
        (1 + parameters_global_1["kappa"] * (out_global_1[["I1"]] + out_global_1[["I2"]]))*
        (parameters_global_1["nu2"] + parameters_global_1["mu"] + parameters_global_1["alpha"])))
    # CHECKING GROWTH RATE USING RT_1
    out_global_1 <- out_global_1 %>%
      mutate(growth_1_check = ((parameters_global_1["nu1"] + parameters_global_1["mu"] + parameters_global_1["alpha"]) * (out_global_1[["Rt_1"]]-1)))
    
    # ADDING SIMULATION COLUMNS TO EXISTING STRUCTURE
aux_mat_global_1 <-cbind(out_global_1[,c(1:21)], rep(beta_2_sensitivity_vector[i],length(tps)), rep(kappa_sensitivity_vector[j],length(tps)), rep(start_I2_sensitivity_vector[q],length(tps)), rep(start_kappa_sensitivity_vector[r],length(tps)))
colnames(aux_mat_global_1)<-c("time","S1","S2","I1","R1","D1", "I2","R2", "D2","FES_1","FES_2","FES_out","logged_values_1","logged_diff_1","growth_rate_1","logged_values_2","logged_diff_2","growth_rate_2","Rt_1","Rt_2","growth_1_check","beta_2","kappa","start_I2","start_kappa")
    # ADDING RESULTS TO RESULT_IJ
    result_ij <- rbind(result_ij, aux_mat_global_1)
      }
    }
    # SAVING RESULTS FOR THIS UNIQUE I, J COMBINATION IN SIMULATION RESULTS FOLDER
    filename <- file.path("simulation_results", sprintf("sim_%03d.csv", file_counter))
    write.csv(result_ij, file = filename, row.names = FALSE)
    cat("Completed and saved results for i =", i, ", j =", j, ", File:", filename, "\n")
    # UPDATING FILE COUNTER
    file_counter <- file_counter + 1
  }
}
cat("Completed simulation. Total files saved:", file_counter - 1, "\n")






# READING CSV FILES TO CREATE SUMMARY BOXPLOT TABLE OF MULTIVARIABLE SENSITIVITY ANALYSIS
# READING ALL CSV FILES FROM THE SIMULATION_RESULTS DIRECTORY
csv_files <- list.files(path = "simulation_results", pattern = "*.csv", full.names = TRUE)

# INITIALISING AN EMPTY DATA FRAME TO STORE SUMMARY RESULTS
summary_table_1 <- data.frame()

# INITIALISING A COUNTER TO TRACK PROGRESS
file_count <- 0
total_files <- length(csv_files)

# LOOPING THROUGH EACH CSV FILE
for (file in csv_files) {
  # USING INCREMENT COUNTER AND PRINTING PROGRESS
  file_count <- file_count + 1
  cat(sprintf("Processing file %d of %d: %s\n", file_count, total_files, basename(file)))
  
  # READING THE CSV FILE
  data_csv <- read.csv(file)
  
  # CHECKING IF REQUIRED COLUMNS EXIST
  required_columns <- c("beta_2", "kappa", "start_I2", "start_kappa", "FES_1", "FES_2", "growth_rate_1")
  missing_columns <- setdiff(required_columns, names(data_csv))
  
  if (length(missing_columns) > 0) {
    cat(sprintf("Warning: Missing columns in file %s: %s\n", basename(file), paste(missing_columns, collapse = ", ")))
    next  # Skip to the next file
  }
  
  # GROUPING BY FOUR PARAMETERS AND CALCULATING MAX FES_1, MAX FES_2, MIN_FES_OUT, MIN GROWTH_RATE_1
  summary <- data_csv %>%
    group_by(beta_2 = as.factor(beta_2), 
             kappa = as.factor(kappa), 
             start_I2 = as.factor(start_I2), 
             start_kappa = as.factor(start_kappa)) %>%
    summarise(max_FES_1 = max(FES_1, na.rm = TRUE),
              max_FES_2 = max(FES_2, na.rm = TRUE),
              min_FES_out = max(FES_1, na.rm = TRUE) / max(FES_2, na.rm = TRUE),
              min_growth_rate_1 = min(growth_rate_1, na.rm = TRUE),
              .groups = 'drop')
  
  # APPENDING TO SUMMARY TABLE
  summary_table_1 <- bind_rows(summary_table_1, summary)
}

# REMOVING ANY DUPLICATE ROWS
summary_table_1 <- distinct(summary_table_1)

# ENSURING WE HAVE ALL COMBINATIONS (11x11x11x11 = 14641 ROWS)
all_combinations <- expand.grid(
  beta_2 = as.factor(seq(0, 1, by = 0.1)),
  kappa = as.factor(seq(0, 1, by = 0.1)),
  start_I2 = as.factor(seq(0, 50, by = 5)),
  start_kappa = as.factor(seq(0, 50, by = 5))
)

# MERGING ALL COMBINATIONS WITH OUR SUMMARY RESULTS
summary_table_1 <- left_join(all_combinations, summary_table_1, 
                             by = c("beta_2", "kappa", "start_I2", "start_kappa"))

# REPLACING NA WITH APPROPRIATE VALUES FOR MAX AND MIN CALCULATIONS
summary_table_1 <- summary_table_1 %>%
  mutate(
    max_FES_1 = replace_na(max_FES_1, -Inf),
    max_FES_2 = replace_na(max_FES_2, -Inf),
    min_FES_out = replace_na(min_FES_out, Inf),
    min_growth_rate_1 = replace_na(min_growth_rate_1, Inf)
  )

# REMOVING ROWS WHERE START_KAPPA < START_I2
summary_table_1 <- summary_table_1 %>%
  mutate(
    start_kappa = as.numeric(as.character(start_kappa)),
    start_I2 = as.numeric(as.character(start_I2))
  ) %>%
  filter(start_kappa >= start_I2) %>%
  mutate(
    start_kappa = as.factor(start_kappa),
    start_I2 = as.factor(start_I2)
  )

# WRITING SUMMARY TABLE TO CSV FILE
write.csv(summary_table_1, file = "simulation_results/summary_table_1.csv", row.names = FALSE)

# PRINTING CONFIRMATION
cat("\nSummary table created and saved as 'summary_table_1.csv' in the simulation_results directory.\n")
cat("Number of rows in summary table:", nrow(summary_table_1), "\n")
print(str(summary_table_1))







# CREATING BOXPLOT SUMMARY GRID
# READ THE SUMMARY TABLE
summary_table_1 <- read.csv("simulation_results/summary_table_1.csv")

# CREATING BOXPLOT FUNCTION FOR A SPECIFIC VARIABLE
create_boxplot <- function(data, var_name, var_values, metric) {
  data %>%
    filter(!!sym(var_name) %in% var_values) %>%
    ggplot(aes(x = !!sym(metric), y = factor(!!sym(var_name), levels = rev(var_values)))) +
    geom_boxplot(outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "blue") +
    theme_minimal() +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, face = "bold", size = 14),
          axis.title.x = element_blank(),
          axis.text.y = element_text(face = "bold", size = 12),
          axis.text.x = element_text(size = 12),
          panel.grid.major.x = element_line(color = "green"),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(t = 20, r = 10, b = 10, l = 50, unit = "pt")) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ylab(var_name)
}

# DEFINING THE CORRECT VALUES FOR EACH VARIABLE
beta_2_values <- c("0", "0.2", "0.4", "0.6", "0.8", "1")
kappa_values <- c("0", "0.2", "0.4", "0.6", "0.8", "1")
start_I2_values <- c("0", "10", "20", "30", "40", "50")
start_kappa_values <- c("0", "10", "20", "30", "40", "50")

# CREATING BOXPLOTS FOR EACH VARIABLE AND METRIC
plot_beta_2_fes <- create_boxplot(summary_table_1, "beta_2", beta_2_values, "min_FES_out")
plot_beta_2_growth <- create_boxplot(summary_table_1, "beta_2", beta_2_values, "min_growth_rate_1")

plot_kappa_fes <- create_boxplot(summary_table_1, "kappa", kappa_values, "min_FES_out")
plot_kappa_growth <- create_boxplot(summary_table_1, "kappa", kappa_values, "min_growth_rate_1")

plot_start_I2_fes <- create_boxplot(summary_table_1, "start_I2", start_I2_values, "min_FES_out")
plot_start_I2_growth <- create_boxplot(summary_table_1, "start_I2", start_I2_values, "min_growth_rate_1")

plot_start_kappa_fes <- create_boxplot(summary_table_1, "start_kappa", start_kappa_values, "min_FES_out")
plot_start_kappa_growth <- create_boxplot(summary_table_1, "start_kappa", start_kappa_values, "min_growth_rate_1")

# COMBINING ALL PLOTS
combined_plot <- (plot_beta_2_fes | plot_beta_2_growth) /
  (plot_kappa_fes | plot_kappa_growth) /
  (plot_start_I2_fes | plot_start_I2_growth) /
  (plot_start_kappa_fes | plot_start_kappa_growth) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Distribution of Final Epidemic Size 1/Final Epidemic Size 2 and Minimum Growth Rate 1",
    subtitle = "For different values of beta_2, kappa, start_I2, and start_kappa",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
                  plot.subtitle = element_text(hjust = 0.5, size = 16))
  )

# ADDING LABELS TO THE COMBINED PLOT
combined_plot <- combined_plot +
  plot_annotation(
    caption = "Four selected parameters on y-axis. Two metrics on x-axis. Boxplots show distributions for a selected parameter value held constant as other parameters vary.",
    theme = theme(plot.caption = element_text(hjust = 0.5, margin = margin(t = 20), size = 14))
  ) &
  theme(plot.margin = margin(t = 40, r = 10, b = 10, l = 10, unit = "pt"))

# ADDING METRIC LABELS AT THE TOP
combined_plot <- combined_plot +
  plot_annotation(
    title = " Column 1: Final Epidemic Size Variant 1/Final Epidemic Size Variant 2     Column 2: Minimum Growth Rate Variant 1",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 20))
    )
  )

# SAVING THE COMBINED PLOT
ggsave("boxplot_summary_combined_No Cross Immunity.png", combined_plot, width = 14, height = 18, dpi = 300)

print("Combined boxplot has been saved as 'boxplot_summary_combined_No Cross Immunity.png'.")

