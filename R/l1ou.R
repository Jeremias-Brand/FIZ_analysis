source('./R/functions.R')
library(ape)
library(l1ou)
library(foreach)
library(tidyverse)
library(writexl)

BS_iterations = 100
cores = 4
runmodel = TRUE
bootstrap = TRUE

if (runmodel) {
  # load the data for sperm only
  comp_data <- readRDS('output/comp_data_main.rds')
  traits <-  data.frame(comp_data[[2]][,c(2,5,4)])
  rownames(traits) = comp_data[[2]]$Short
  new_obj = adjust_data(comp_data[[1]], traits)
  
  mod0 <- estimate_shift_configuration(new_obj$tree, new_obj$Y,
                                       criterion = "AICc",
                                       root.model = "OUrandomRoot")
  
  saveRDS(mod0 , 'output/main_mod0_AICc.rds')
  
  ################################################################################
  ## total_sprmLen_um
  traits <-  data.frame(comp_data[[2]][,c(2)])
  rownames(traits) = comp_data[[2]]$Short
  new_obj = adjust_data(comp_data[[1]], traits)
  
  mod2 <- estimate_shift_configuration(
    new_obj$tree, new_obj$Y, criterion = "AICc", root.model = "OUrandomRoot")
  saveRDS(mod2 , 'output/main_mod2_sprmLen_AICc.rds')
  ################################################################################
  ## Anterior sperm length
  traits <-  data.frame(comp_data[[2]][,c(3)])
  rownames(traits) = comp_data[[2]]$Short
  new_obj = adjust_data(comp_data[[1]], traits)
  
  mod3 <- estimate_shift_configuration(
    new_obj$tree, new_obj$Y, criterion = "AICc", root.model = "OUrandomRoot")
  saveRDS(mod3 , 'output/main_mod3_antSprmLen_AICc.rds')
  ################################################################################
  ## brstlLen
  traits <-  data.frame(comp_data[[2]][,c(4)])
  rownames(traits) = comp_data[[2]]$Short
  new_obj = adjust_data(comp_data[[1]], traits)
  
  mod4 <- estimate_shift_configuration(
    new_obj$tree, new_obj$Y, criterion = "AICc", root.model = "OUrandomRoot")
  saveRDS(mod4 , 'output/main_mod4_brstlLen_AICc.rds')
  ################################################################################
  ## sperm ratio
  traits <-  data.frame(comp_data[[2]][,c(5)])
  rownames(traits) = comp_data[[2]]$Short
  new_obj = adjust_data(comp_data[[1]], traits)
  
  mod5 <- estimate_shift_configuration(
    new_obj$tree, new_obj$Y, criterion = "AICc", root.model = "OUrandomRoot")
  saveRDS(mod5 , 'output/main_mod5_sprmratio_AICc.rds')
} else {
  mod0 <- readRDS('output/main_mod0_AICc.rds')
}

################################################################################

## sperm length only big set

################################################################################

comp_data <- readRDS('output/comp_data_sprmLen_only.rds')
traits <-  data.frame(comp_data[[2]][,2])
rownames(traits) = comp_data[[2]]$Short
new_obj = adjust_data(comp_data[[1]], traits)

# l1ou with sperm length
mod6 <- estimate_shift_configuration(new_obj$tree, new_obj$Y,
                                     criterion = "AICc",
                                     root.model = "OUrandomRoot",
                                     alpha.upper = 500)

saveRDS(mod6 , 'output/sprmLen_only_mod0_AICc.rds')

mod6 <- readRDS('output/sprmLen_only_mod0_AICc.rds')

################################################################################

## 

################################################################################
mod_lst = list(mod0,
               mod2,
               mod3,
               mod4,
               mod5,
               mod6)

var_names = list(
  c('total sperm length', 'sperm ratio', 'bristle length'),
  c('total sperm length'),
  c('anterior sperm length'),
  c('bristle length'),
  c('sperm ratio'),
  c('total sperm length')
)

mod_names <-  c(
   'full model',
  'total sperm length',
  'anterior sperm length',
  'bristle length',
  'sperm ratio',
  'total sperm length - more species'
)

names(mod_lst) <- mod_names

# get the number of models within delta 10 AICc
included_models <- lapply(mod_lst, profile_plot) 

################################################################################
## For loop to print diagnostics and results.
res_par = list()
for (i in 1:length(mod_lst)) {
  mod = mod_lst[[i]]
  tmp = list()
  for (j in 1:included_models[[i]]){
    mod_i <- fit_OU(mod$tree, mod$Y, 
                  mod$profile$configurations[[j]], l1ou.options=mod$l1ou.options)
    tmp[[j]] <- get_result(mod_i, var_names[[i]])
  }
  tmp2 <- as.data.frame(do.call(rbind, tmp))
  tmp2$data <- mod_names[[i]]
  nrow(tmp2)/4
  tmp2$parameter <- rep(x = c("adaptation rate (alpha)", "variance (sigma2)", 
                              "stationary variance (gamma)", "logLik",
                              "AICc", "nShifts"), included_models[[i]])
  
  tmp2$model <- rep(1:included_models[[i]], each = 6)
  # tmp2$AICc <- mod$score
  # tmp2$nShifts <- mod$nShifts
  res_par[[i]] <- tmp2
}

df <- res_par[[2]]

for (i in 2:6) {
  df <- res_par[[i]]
  res_par[[i]] <- df |> 
    pivot_wider(values_from = 1, names_from = parameter) |> 
    select(data, model, nShifts, AICc, logLik, everything())
}

out_lst <- res_par

names(out_lst) <- mod_names

# output xlsx
write_xlsx(
  x = out_lst,
  path = 'output/results_l1ou.xlsx'
)

################################################################################

## Bootstrap

################################################################################

outname = 'output/main_mod0_AICc_'

# Bootstrapping the models
for (i in 1:2) {
  bs_model(mod0, i, BS_iterations, cores, outname)
}


