source('./R/functions.R')
library(phylolm)
library(tidyverse)
library(writexl)

comp_data <- readRDS('data/comp_data_main.rds')

t_clade = c(
  'Mac056', 'Mac057', 'Mac058', 'Mac059', 'Mac060',
  'Mac061', 'Mac064', 'Mac065', 'Mac067', 'Mac068',
  'Mac069', 'Mac070', 'Mac071', 'Mac072', 'Mac074',
  'Mac077', 'Mac117')

tre <- comp_data[[1]]
df <- comp_data[[2]]
df$clade <- ifelse(df$Short %in% t_clade, 'Tanganyika', 'Other')
row.names(df) <- df$Short

responses <- c('total_sprmLen_um', "antSpermLength_um", "sperm_ratio", "brstlLen_um_per_sperm")
res <- lapply(responses, fit_lambda, df, tre)

out <- do.call(rbind, res)

write_xlsx(out, path = 'output/PGLS.xlsx')