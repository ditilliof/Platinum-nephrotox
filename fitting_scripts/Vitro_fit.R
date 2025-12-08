#Load the libraries

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


#options(error = traceback)

# Get the current date and time
current_time <- Sys.time()

# Path to the Stan model file
stan_model_path <- "../Stan_files/invitroqAOP.stan"

# Format the current date and time as a time stamp
time_stamp <- format(current_time, "%Y-%m-%d_%H-%M-%S")

# Prompt the user for used data (choice between IGS and EGs)
method <- readline("Enter data choice (choose between IGS and EGs): ")

# Create the folder name using the time stamp
folder_name <- paste0("../Stanresults/invitroqAOP/", "/", time_stamp)

# Create the folder
dir.create(folder_name)

file.copy("../Stan_files/invitroqAOP.stan", 
          file.path(folder_name, "invitroqAOP.stan"))
dls = c(1.0,2.5,5.0,10.0,20.0,30.0,50.0)


if (method == IGS) {
#   #Using IGS
  CISQAOP_FDT_RPTECdata <- read_csv("../data/in vitro/IGS_invitroqAOP_20240822.csv") %>%
    select(-1) %>%
    filter(DOSE %in% dls, StateVar == "DD")%>% group_by(StateVar, TIME, DOSE) %>%
    summarize(mean = mean(Score, na.rm = TRUE), sd = sd(Score, na.rm = TRUE)) %>%
    ungroup()
} else if (method == EGs) {
  # #Using EGs
  CISQAOP_FDT_RPTECdata <- read.csv("../data/in vitro/invitro160.csv") %>%
    filter(DOSE %in% dls)
}

Necro_AOP = read_csv("../data/in vitro/NecroAOPnocc_fulltime.csv") %>% select(-1) %>%
  filter(DOSE %in% dls, TIME > 20)

pars = c("k_cDD","degDD","hillDD","k_inter","maxdeath","k_hillnec","p","h")

Fina_value = 6.2 * 3600     #flow rate cisplatin from apical medium into cell in um3/h 
Finb_value = 0.64 * 3600    #flow rate cisplatin from basolateral medium into cell in um3/h
Fmin_value = 0.04 * 3600    #flow rate metabolite into cell in um3/h

data_list <- list(
  N_DD = length(unique(CISQAOP_FDT_RPTECdata$TIME)),
  N_NEC = length(unique(Necro_AOP$TIME)),
  M = length(unique(CISQAOP_FDT_RPTECdata$DOSE)),
  t0 = 0,
  ts = unique(CISQAOP_FDT_RPTECdata$TIME),
  ts1 = unique(Necro_AOP$TIME),
  tot = sort(unique(c(unique(CISQAOP_FDT_RPTECdata$TIME), unique(Necro_AOP$TIME)))),
  i_ts = match(unique(CISQAOP_FDT_RPTECdata$TIME), sort(unique(c(unique(CISQAOP_FDT_RPTECdata$TIME), unique(Necro_AOP$TIME))))),
  i_ts1 = match(unique(Necro_AOP$TIME), sort(unique(c(unique(CISQAOP_FDT_RPTECdata$TIME), unique(Necro_AOP$TIME))))),
  
  y_DD = as.matrix(CISQAOP_FDT_RPTECdata %>%
                     dplyr::select(-sd) %>%
                     pivot_wider(names_from = "DOSE", values_from = mean) %>%
                     dplyr::select(-c(TIME, StateVar))),
  y_NEC = as.matrix(Necro_AOP %>%
                      dplyr::select(-sd) %>%
                      pivot_wider(names_from = "DOSE", values_from = mean) %>%
                      dplyr::select(-c(TIME, StateVar))),
  sigma_DD = as.matrix(CISQAOP_FDT_RPTECdata %>%
                         dplyr::select(-mean) %>%
                         pivot_wider(names_from = "DOSE", values_from = sd) %>%
                         dplyr::select(-c(TIME, StateVar))),
  sigma_NEC = as.matrix(Necro_AOP %>%
                          dplyr::select(-mean) %>%
                          pivot_wider(names_from = "DOSE", values_from = sd) %>%
                          dplyr::select(-c(TIME, StateVar))),
  
  NEC0 = 0,
  Vapi = 1e12,
  Vbas = 2e12,
  Vcell = 2005,
  Ncell = 2e6,
  Fina = Fina_value,
  Fouta = Fina_value / 1.8,
  Finb = Finb_value,
  Foutb = Finb_value / 51,
  Kmet = 8.4 * 3600,
  Fmin = Fmin_value,
  Fmout = Fmin_value / 100,
  
  Qapi3 = dls[1] * 301.1 * 1e-15 * 1e12,
  Qbas3 = dls[1] * 301.1 * 1e-15 * 2e12,
  Qapi4 = dls[2] * 301.1 * 1e-15 * 1e12,
  Qbas4 = dls[2] * 301.1 * 1e-15 * 2e12,
  Qapi5 = dls[3] * 301.1 * 1e-15 * 1e12,
  Qbas5 = dls[3] * 301.1 * 1e-15 * 2e12,
  Qapi6 = dls[4] * 301.1 * 1e-15 * 1e12,
  Qbas6 = dls[4] * 301.1 * 1e-15 * 2e12,
  Qapi7 = dls[5] * 301.1 * 1e-15 * 1e12,
  Qbas7 = dls[5] * 301.1 * 1e-15 * 2e12,
  Qapi8 = dls[6] * 301.1 * 1e-15 * 1e12,
  Qbas8 = dls[6] * 301.1 * 1e-15 * 2e12,
  Qapi9 = dls[7] * 301.1 * 1e-15 * 1e12,
  Qbas9 = dls[7] * 301.1 * 1e-15 * 2e12,
  
  Qcell0 = 0,
  DD0 = 0,
  Qinter0 = 0
)



options(mc.cores = parallel::detectCores())  # Use all available cores
#options(mc.cores = 1)
chains <- 4 # Number of chains
iter <- 500 # Number of iterations
warmup <- 400 # Number of warmup iterations
thin <- 1  # Thinning parameter


compiled_model <- cmdstan_model(stan_file = stan_model_path, force_recompile = TRUE)


fit_invitro <- compiled_model$sample(
  data = data_list,
  parallel_chains = getOption("mc.cores", 4),
  chains = chains,
  iter_warmup = warmup,
  iter_sampling = iter - warmup,
  output_dir = folder_name,
  refresh = 1
)



color_scheme_set("mix-blue-red")
fit_trace <- mcmc_trace(fit_invitro$draws(format="df",variables = pars)) # 
tiff(paste0(folder_name,"/", "Trace-plot_", time_stamp, ".png"), units="in", width=9, height=5.3, res=700, compression = 'lzw')
print(fit_trace)
dev.off()


fit_density <- mcmc_dens_overlay(fit_invitro$draws(format="df", variables = pars))
tiff(paste0(folder_name,"/", "Density-plot_", time_stamp, ".png"), units="in", width=9, height=5.3, res=700, compression = 'lzw')
print(fit_density)
dev.off()


draws_cism <- fit_invitro$draws(format = "draws_matrix",variables = pars)
fit_summary <- fit_invitro$summary()


sumParcism<-data.frame(fit_summary)%>% textshape::column_to_rownames("variable") %>% 
  select(1,3,7)



parSetsm<- draws_cism


filepath_draws = paste0(folder_name, "/draws_cisDDNEC")
filepath_sum = paste0(folder_name, "/sumPar_cisDDNEC")


write.csv(draws_cism, paste0(filepath_draws, time_stamp, ".csv"))
write.csv(sumParcism, paste0(filepath_sum, time_stamp, ".csv"))
saveRDS(fit_invitro, file = paste0(folder_name, "/", "fitvitro.rds"))


loo_output3 = fit_invitro$loo(variables = "logLikelihood_total3")
loo_output4 = fit_invitro$loo(variables = "logLikelihood_total4")
loo_output5 = fit_invitro$loo(variables = "logLikelihood_total5")
loo_output6 = fit_invitro$loo(variables = "logLikelihood_total6")
loo_output7 = fit_invitro$loo(variables = "logLikelihood_total7")
loo_output8 = fit_invitro$loo(variables = "logLikelihood_total8")
loo_output9 = fit_invitro$loo(variables = "logLikelihood_total9")

saveRDS(loo_output3, paste0(folder_name, "/loo_output3_", time_stamp, ".rds"))
saveRDS(loo_output4, paste0(folder_name, "/loo_output4_", time_stamp, ".rds"))
saveRDS(loo_output5, paste0(folder_name, "/loo_output5_", time_stamp, ".rds"))
saveRDS(loo_output6, paste0(folder_name, "/loo_output6_", time_stamp, ".rds"))
saveRDS(loo_output7, paste0(folder_name, "/loo_output7_", time_stamp, ".rds"))
saveRDS(loo_output8, paste0(folder_name, "/loo_output8_", time_stamp, ".rds"))
saveRDS(loo_output9, paste0(folder_name, "/loo_output9_", time_stamp, ".rds"))




