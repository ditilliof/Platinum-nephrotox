# Load the necessary libraries
library("tidyverse")
library("colorspace")
library("scales")
library("data.table")
library("readxl")
library("stringr")
library("stats")
library("parallel")
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")
library(lubridate)
library(loo)
library(conflicted)

conflict_prefer("filter", "dplyr")

# Set the path to the Stan model file
stan_model_path <- "../Stan_files/invivoqAOP.stan"

# Compile the model
compiled_model <- cmdstan_model(stan_model_path)

# Save the compiled executable path
compiled_exe_path <- compiled_model$exe_file()

# Create the folder for results
current_time <- Sys.time()

time_stamp <- format(current_time, "%Y-%m-%d_%H-%M-%S")

folder_name <- paste0("../Stanresults/invivoqAOP/", "/", time_stamp)

dir.create(folder_name)

file.copy("../Stan_files/invivoqAOP.stan", 
          file.path(folder_name, "invivoqAOP.stan"))

# Prompt the user for used data (choice between IGS and EGs)
method <- readline("Enter data choice (choose between IGS and EGs): ")

# Load the data based on user input
if (method == "IGS") {
  CISQAOP_FDT_vivo <- read_csv("../data/in vivo/CISQAOP_FDT_2024_invivodata_IGS.csv") %>%
    select(-1) %>%
    mutate(sds = ifelse(sds == 0, sds + 0.05, sds))
  
} else if (method == "EGs") {
  CISQAOP_FDT_vivo <- read_csv("../data/in vivo/CISQAOP_FDT_2024_invivodata_EGs2.csv") %>%
    select(-1) %>%
    filter(TIMEPOINT != 192) %>%
    mutate(sds = ifelse(sds == 0, sds + 0.05, sds))
}

# Prepare data for Stan model
PK_pars <- read.csv("../Stanresults/PK3/2024-07-22_15-15-34/fit_pk_summary.csv") %>%
  select(-1) %>%
  column_to_rownames(var = "variable")

pars = c("k_kidDD", "d_DD", "maxdeath", "hillDD", "k_hillnec","d_CD","k_CDINF",
         "d_INF","k_INFKF","hillKF","h1","INF0") #add p0 in case you want to estimate it

data_list <- list(
  NDD = length(unique((CISQAOP_FDT_vivo %>% filter(StateVar == "DD"))$TIMEPOINT)),
  NINF = length(unique((CISQAOP_FDT_vivo %>% filter(StateVar == "INF"))$TIMEPOINT)),
  NCD = length(unique((CISQAOP_FDT_vivo %>% filter(StateVar == "CD"))$TIMEPOINT)),
  NKF = length(unique((CISQAOP_FDT_vivo %>% filter(StateVar == "KF"))$TIMEPOINT)),
  NTOT = length(unique(CISQAOP_FDT_vivo$TIMEPOINT)),
  t0 = 0,
  tsdd = unique((CISQAOP_FDT_vivo %>% filter(StateVar == "DD"))$TIMEPOINT),
  tsinf = unique((CISQAOP_FDT_vivo %>% filter(StateVar == "INF"))$TIMEPOINT),
  tscd = unique((CISQAOP_FDT_vivo %>% filter(StateVar == "CD"))$TIMEPOINT),
  tskf = unique((CISQAOP_FDT_vivo %>% filter(StateVar == "KF"))$TIMEPOINT),
  tot = sort(unique(CISQAOP_FDT_vivo$TIMEPOINT)),
  idd = base::match(unique((CISQAOP_FDT_vivo %>% filter(StateVar == "DD"))$TIMEPOINT), sort(unique(CISQAOP_FDT_vivo$TIMEPOINT))),
  iinf = base::match(unique((CISQAOP_FDT_vivo %>% filter(StateVar == "INF"))$TIMEPOINT), sort(unique(CISQAOP_FDT_vivo$TIMEPOINT))),
  icd = base::match(unique((CISQAOP_FDT_vivo %>% filter(StateVar == "CD"))$TIMEPOINT), sort(unique(CISQAOP_FDT_vivo$TIMEPOINT))),
  ikf = base::match(unique((CISQAOP_FDT_vivo %>% filter(StateVar == "KF"))$TIMEPOINT), sort(unique(CISQAOP_FDT_vivo$TIMEPOINT))),
  ydd = as.matrix(CISQAOP_FDT_vivo %>% filter(StateVar == "DD") %>% dplyr::select(-sds) %>%
                    group_by(TIMEPOINT) %>%
                    pivot_wider(names_from = c("StateVar"), values_from = means) %>%
                    ungroup() %>% dplyr::select(-c(TIMEPOINT))),
  yinf = as.matrix(CISQAOP_FDT_vivo %>% filter(StateVar == "INF") %>% dplyr::select(-sds) %>%
                     group_by(TIMEPOINT) %>%
                     pivot_wider(names_from = c("StateVar"), values_from = means) %>%
                     ungroup() %>% dplyr::select(-c(TIMEPOINT))),
  ycd = as.matrix(CISQAOP_FDT_vivo %>% filter(StateVar == "CD") %>% dplyr::select(-sds) %>%
                    group_by(TIMEPOINT) %>%
                    pivot_wider(names_from = c("StateVar"), values_from = means) %>%
                    ungroup() %>% dplyr::select(-c(TIMEPOINT))),
  ykf = as.matrix(CISQAOP_FDT_vivo %>% filter(StateVar == "KF") %>% dplyr::select(-sds) %>%
                    group_by(TIMEPOINT) %>%
                    pivot_wider(names_from = c("StateVar"), values_from = means) %>%
                    ungroup() %>% dplyr::select(-c(TIMEPOINT))),
  sigmadd = as.matrix(CISQAOP_FDT_vivo %>% filter(StateVar == "DD") %>% dplyr::select(-means) %>%
                        group_by(TIMEPOINT) %>%
                        pivot_wider(names_from = c("StateVar"), values_from = sds) %>%
                        ungroup() %>% dplyr::select(-c(TIMEPOINT))),
  sigmainf = as.matrix(CISQAOP_FDT_vivo %>% filter(StateVar == "INF") %>% dplyr::select(-means) %>%
                         group_by(TIMEPOINT) %>%
                         pivot_wider(names_from = c("StateVar"), values_from = sds) %>%
                         ungroup() %>% dplyr::select(-c(TIMEPOINT))),
  sigmacd = as.matrix(CISQAOP_FDT_vivo %>% filter(StateVar == "CD") %>% dplyr::select(-means) %>%
                        group_by(TIMEPOINT) %>%
                        pivot_wider(names_from = c("StateVar"), values_from = sds) %>%
                        ungroup() %>% dplyr::select(-c(TIMEPOINT))),
  sigmakf = as.matrix(CISQAOP_FDT_vivo %>% filter(StateVar == "KF") %>% dplyr::select(-means) %>%
                        group_by(TIMEPOINT) %>%
                        pivot_wider(names_from = c("StateVar"), values_from = sds) %>%
                        ungroup() %>% dplyr::select(-c(TIMEPOINT))),
  Plasma0 = PK_pars[["Plasma0", 'mean']],
  Kid0 = 0,
  Accu0 = 0,
  DD0 = 0,
  CD0 = 0,
  KF0 = 0,
  p0 = 10, #comment if want to estimate initial state
  k_PlasKid = PK_pars[['k_PlasKid', 'mean']],
  k_ePlas = PK_pars[['k_ePlas', 'mean']],
  k_KidPlas = PK_pars[['k_KidPlas', 'mean']],
  k_eKid = PK_pars[['k_eKid', 'mean']],
  k_KidAccu = PK_pars[['k_KidAccu', 'mean']],
  k_AccuKid = PK_pars[['k_AccuKid', 'mean']],
  scale = PK_pars[['scale', 'mean']],
  d_p = 0.005
  #d_p = 0 #case with no decay
)

random_seed <- sample.int(.Machine$integer.max, 1)
cat("Random Seed:", random_seed, "\n")

options(mc.cores = parallel::detectCores())  # Use all available cores
chains <- 4 # Number of chains
iter <- 10000  # Number of iterations
warmup <- 6000  # Number of warmup iterations
thin <- 1  # Thinning parameter

# Load the precompiled model using the executable path
compiled_exe_path <- "../Stan_files/invivoqAOP"
compiled_model <- cmdstan_model(stan_file = stan_model_path, exe_file = compiled_exe_path)

# This is for cmdstanr
fit_cisplatin_vivo <- compiled_model$sample(data = data_list, 
                                            parallel_chains = getOption("mc.cores", 4), 
                                            chains = chains,
                                            iter_warmup = warmup, 
                                            iter_sampling = iter - warmup,
                                            output_dir = folder_name,
                                            refresh = 1,
                                            seed = random_seed
)


color_scheme_set("mix-blue-red")
fit_cisvivo_trace <- mcmc_trace(fit_cisplatin_vivo$draws(format="df",variables = pars))
tiff(paste0(folder_name,"/", "Trace-plot_", time_stamp, ".png"), units="in", 
     width=9, height=5.3, res=700, compression = 'lzw')
print(fit_cisvivo_trace)
dev.off()


fit_cisvivo_density <- mcmc_dens_overlay(fit_cisplatin_vivo$draws(format="df",variables = pars))
tiff(paste0(folder_name,"/", "Density-plot_", time_stamp, ".png"), units="in", 
     width=9, height=5.3, res=700, compression = 'lzw')
print(fit_cisvivo_density)
dev.off()


draws_cisvivo <- fit_cisplatin_vivo$draws(format = "draws_matrix",variables = pars)
fit_cisvivo_summary <- fit_cisplatin_vivo$summary()




sumParcisvivo<-data.frame(fit_cisvivo_summary)%>% textshape::column_to_rownames("variable") %>% 
  select(1,3,7)



parSetsvivo<- draws_cisvivo


filepath_draws = paste0(folder_name, "/draws_cisvivo")
filepath_sum = paste0(folder_name, "/sumPar_cisvivo")


write.csv(draws_cisvivo, paste0(filepath_draws, time_stamp, ".csv"))
write.csv(sumParcisvivo, paste0(filepath_sum, time_stamp, ".csv"))


log_lik = fit_cisplatin_vivo$draws("logLikelihood_total", format = "matrix")
loo_result <- loo(log_lik)

filepath_loo = paste0(folder_name, "/loo_result")
saveRDS(loo_result, paste0(filepath_loo, time_stamp, ".rds"))
saveRDS(fit_cisplatin_vivo, file = paste0(folder_name, "/", "fitvivo.rds"))



