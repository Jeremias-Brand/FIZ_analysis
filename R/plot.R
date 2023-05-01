source('./R/functions.R')
library(l1ou)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

# similar pairs
col = c('#9800D3', '#CC9C8D',
'#F9AE54', '#B1CC97',
'#FA9584', '#749696',
'#C4E2E2', '#FD79B9',
'#688052', '#FD7849',
'#806158', '#E00B7A')


cols = col

models = c(
  'mod0_AICc',
  'mod2_sprmLen_AICc',
  'mod3_antSprmLen_AICc',
  'mod4_brstlLen_AICc', 
  'mod5_sprmratio_AICc'
)

for (i in 1:length(models)) {
  mod0 <- readRDS(paste0('output/main_', models[[i]], '.rds'))
  
  pdf(paste0('fig/main_', models[[i]], '_profile.pdf'))
  passed_models <- profile_plot(mod0)
  dev.off()
  
  pdf(paste0('fig/main_', models[[i]], '_diagnostic.pdf'), height = 9)
  plot_diag(mod0 = mod0, passed_models = passed_models)
  dev.off()
  
}


mod0 <- readRDS('output/sprmLen_only_mod0_AICc.rds')

pdf('fig/sprmLen_only_mod0_AICc_profile.pdf')
passed_models <- profile_plot(mod0)
dev.off()

pdf('fig/main_sprmLen_only_mod0_AICc_diagnostic.pdf', height = 9)
plot_diag(mod0 = mod0, passed_models = passed_models)
dev.off()

################################################################################

## main model diagnostics

################################################################################

pdf('fig/sprmLen_only_mod0_AICc_diagnostic.pdf')
mod0 <- readRDS('data/sprmLen_only_mod0_AICc.rds')

passed_models <- profile_plot(mod0)
for (index in 1:passed_models) {
  mod <- fit_OU(mod0$tree, mod0$Y, 
                mod0$profile$configurations[[index]], l1ou.options=mod0$l1ou.options)
  
  mod$nShifts
  if (mod$nShifts <= length(cols)) {
    pal = c(cols[1:mod$nShifts],'#808080')
  } else {
    pal = c(brewer.pal(mod$nShifts, 'Spectral'),'#808080')
  }
  
  plot(mod, plot.bar = FALSE,  asterisk = TRUE, edge.label.ann = FALSE,
       edge.shift.ann = FALSE,
       show.tip.label = FALSE)
}
dev.off()




