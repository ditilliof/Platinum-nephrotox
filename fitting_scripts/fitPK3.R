library(cmdstanr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(loo)
library(posterior)
library(bayesplot)
library(conflicted)

conflict_prefer("filter", "dplyr")

# Get the current date and time
current_time <- Sys.time()

# Format the current date and time as a time stamp
time_stamp <- format(current_time, "%Y-%m-%d_%H-%M-%S")

# Create the folder name using the time stamp
folder_name <- paste0("../Stanresults/PK3/", time_stamp)

# Create the folder
dir.create(folder_name)

# Load the data
full_S <- read.csv("../data/in vivo/full_S.csv") %>% select(-1) %>% mutate(StateVar="S") %>% filter(Time!=0)

full_Accu <- read.csv("../data/in vivo/full_Accu.csv") %>%
  select(-1) %>%
  mutate(StateVar="Accu") %>%
  filter(!(Time == 0 | Time == 2 | Time == 4))

full_PK <- rbind(semi_join(full_S,full_Accu, by = "Time"),full_Accu)
pars = c("k_PlasKid", "k_ePlas", "k_KidPlas", "k_eKid", "k_AccuKid", "k_KidAccu","scale")
ic = c("Plasma0")

# Prepare data list for Stan model
data_list <- list(
  N1 = length(unique(full_S$Time)),
  N2 = length(unique(full_Accu$Time)),
  t0 = 0,
  ts1 = unique(full_S$Time),
  ts2 = unique(full_Accu$Time),
  y1 = as.matrix(full_S %>% dplyr::select(-SD) %>% 
                  pivot_wider(names_from = "StateVar", values_from = Score) %>%
                  dplyr::select(-Time)),
  y2 = as.matrix(full_Accu %>% dplyr::select(-SD) %>% 
                   pivot_wider(names_from = "StateVar", values_from = Score) %>%
                   dplyr::select(-Time)),
  sigma1 = as.matrix(full_S %>% dplyr::select(-Score) %>% 
                       pivot_wider(names_from = "StateVar", values_from = SD) %>%
                       dplyr::select(-Time)),
  sigma2 = as.matrix(full_Accu %>% dplyr::select(-Score) %>%
                       pivot_wider(names_from = "StateVar", values_from = SD) %>%
                       dplyr::select(-Time)),
  KidneyPt0 = 0,
  AccuPt0 = 0
)

# Set options for parallel computing
options(mc.cores = parallel::detectCores())

# Model parameters
chains <- 4
iter <- 10000
warmup <- 6000
thin <- 1

# Compile the Stan model
compiled_model <- cmdstan_model("../Stan_files/PK3.stan")

# Fit the model using cmdstanr
fit_pk_3 <- compiled_model$sample(
  data = data_list,
  parallel_chains = chains,
  iter_warmup = warmup,
  iter_sampling = iter - warmup,
  output_dir = folder_name,
  refresh = 1,
  seed = 1234
)

# Print folder name to debug
print(folder_name)

# Construct and print the final file path to debug
final_file_path <- paste0(folder_name, "/fit_pk_3.rds")
print(final_file_path)

# Save the fitted model object
saveRDS(fit_pk_3, file = final_file_path)

# Generate and save trace and density plots
fit_pk_trace <- mcmc_trace(fit_pk_3$draws(variables = c(pars, ic)))
fit_pk_density <- mcmc_dens_overlay(fit_pk_3$draws(variables = c(pars, ic)))

# Save the plots correctly
ggsave(filename = paste0(folder_name, "/fit_pk_trace.png"), plot = fit_pk_trace)
ggsave(filename = paste0(folder_name, "/fit_pk_density.png"), plot = fit_pk_density)

# Get and save the summary of the model fit
fit_pk_summary <- fit_pk_3$summary()
write.csv(fit_pk_summary, paste0(folder_name,"/fit_pk_summary.csv"))

# Extract draws for further analysis
draws_pk <- fit_pk_3$draws(format = "draws_list", variables = c(pars, ic))
write.csv(draws_pk, paste0(folder_name, "/draws_pk.csv"))
saveRDS(fit_pk_3, file = paste0(folder_name, "/", "fitPK3.rds"))
