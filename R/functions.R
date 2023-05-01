fit_lambda <- function(response, df, tre) {
  model <- phylolm(as.formula(paste0(response , '~ clade')), data = df, phy = tre, model="lambda")
  model_summary <- summary(model)
  
  # Coefficients
  coefficients_df <- model_summary$coefficients %>% as.data.frame() %>% rownames_to_column("Predictor")
  coefficients_df$AIC <- model_summary$aic
  coefficients_df$df <- model_summary$df
  coefficients_df$logLik <- model_summary$logLik
  coefficients_df$lambda <- model_summary$optpar
  coefficients_df$sigma2 <- model_summary$sigma2
  
  coefficients_df$Response <- response
  
  coefficients_df <- coefficients_df |> 
    select(Response, AIC, df, logLik, lambda, sigma2, Predictor:p.value) |> 
    mutate(across(Estimate:t.value, \(x) round(x, 2)),
           p.value = format.pval(p.value, eps = 0.001, digits = 1),
           across(AIC:sigma2, \(x) round(x, 2)))
  
  return(coefficients_df)
}


force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}


bs_model <- function(mod0, index, BS_iterations, cores, outname) {
  mod_i = fit_OU(mod0$tree, mod0$Y, 
                 mod0$profile$configurations[[index]], l1ou.options=mod0$l1ou.options)
  mod_i_bs <- l1ou_bootstrap_support(
    mod_i, nItrs = BS_iterations, multicore = TRUE, 
    nCores = cores, quietly = FALSE)
  saveRDS(mod_i, paste0(outname, index, '.rds'))
  saveRDS(mod_i_bs, paste0(outname, index, '_bs.rds'))
}

get_result <- function(mod, var_names) {
  tmp.mat = rbind(mod$alpha, mod$sigma2, mod$sigma2/(2 * 
                                                       mod$alpha), mod$logLik,
                  mod$score, mod$nShifts)
  rownames(tmp.mat) = c("adaptation rate (alpha)", "variance (sigma2)", 
                        "stationary variance (gamma)", "logLik",
                        "AICc", "nShifts")
  colnames(tmp.mat) = var_names
  return(round(tmp.mat,2))
}

profile_plot <- function(mod, cutoff = 10) {
  mod.profile <- profile(mod)
  df <- data.frame(nshifts = unlist(mod.profile$nShifts),
                   scores = unlist(mod.profile$scores)) %>% 
    mutate(sig = case_when(
      scores < min(scores) + cutoff ~ TRUE,
      TRUE ~ FALSE
    ))
}

get_clades_for_regime <- function(row_regime, t) {
  # check if we end in a tip node
  endnode = t$edge[row_regime,2]
  if (endnode < t$Nnode) {
    return(t$tip.label[endnode])
  }
  subclade <- extract.clade(t, endnode)
  return(subclade$tip.label)
}

# plot precomputed models
plot_bs <- function(modRDS, bsRDS, cutoff = 0.10, cols = NA) {
  mod <- readRDS(modRDS)
  bs <- readRDS(bsRDS)
  
  nEdges <- Nedge(mod$tree)
  e.w <- rep(1,nEdges)
  e.w[mod$shift.configuration] <- 2
  e.l <- round(bs$detection.rate, digits=2)
  # to avoid annotating edges with support at or below 10%
  e.l <- ifelse(e.l>cutoff, e.l, NA)
  
  mod$nShifts
  if (mod$nShifts <= length(cols)) {
    pal = c(cols[1:mod$nShifts],'#808080')
  } else {
    pal = c(brewer.pal(mod$nShifts, 'Spectral'),'#808080')
  }
  plot(mod, 
       palette = pal,
       edge.label=e.l, edge.ann.cex=1.5,
       edge.label.ann=TRUE, cex=1.2, label.offset=0.02,
       edge.label.adj = c(0.5, 1),
       edge.width=e.w * 2, 
       asterisk = FALSE,
       edge.shift.ann = FALSE)
  
}


plot_top_model <- function(modRDS, cutoff = 0.10, cols = NA) {
  mod <- readRDS(modRDS)
  
  nEdges <- Nedge(mod$tree)
  e.w <- rep(1,nEdges)
  e.w[mod$shift.configuration] <- 2
  e.l <- round(bs$detection.rate, digits=2)
  # to avoid annotating edges with support at or below 10%
  e.l <- ifelse(e.l>cutoff, e.l, NA)
  
  mod$nShifts
  if (mod$nShifts <= length(cols)) {
    pal = c(cols[1:mod$nShifts],'#808080')
  } else {
    pal = c(brewer.pal(mod$nShifts, 'Spectral'),'#808080')
  }
  plot(mod, 
       palette = pal,
       edge.label=e.l, edge.ann.cex=1.5,
       edge.label.ann=TRUE, cex=1.2, label.offset=0.02,
       edge.label.adj = c(0.5, 1),
       edge.width=e.w * 2, 
       asterisk = FALSE,
       edge.shift.ann = FALSE)
  
}


profile_plot <- function(mod, cutoff = 10) {
  mod.profile <- profile(mod)
  df <- data.frame(nshifts = unlist(mod.profile$nShifts),
                   scores = unlist(mod.profile$scores)) %>% 
    mutate(sig = case_when(
      scores < min(scores) + cutoff ~ TRUE,
      TRUE ~ FALSE
    ))
  
  cc = c('#9800D3', 
         '#F9AE54', '#B1CC97',
         '#FA9584', '#749696',
         '#C4E2E2', '#FD79B9',
         '#688052', '#FD7849',
         '#E00B7A')
  
  p <- df %>% 
    filter(scores < 1000) %>% 
    ggplot(aes(x = nshifts, y = scores, color = sig)) +
    geom_hline(yintercept = min(df$scores) + cutoff) +
    geom_point(size = 4) +
    scale_color_manual(values = cc) +
    scale_x_continuous(breaks = seq(0,40,5)) +
    xlab('Number of shifts') +
    ylab('AICc') +
    theme_minimal() +
    theme(legend.position = 'none') +
    theme(text = element_text(size = 20))
  
  print(p)
  
  return(nrow(df[df$sig == TRUE,]))
}

plot_diag <- function(mod0, passed_models) {
  for (index in 1:passed_models) {
    mod <- fit_OU(mod0$tree, mod0$Y, 
                  mod0$profile$configurations[[index]], l1ou.options=mod0$l1ou.options)
    
    mod$nShifts
    if (mod$nShifts <= length(cols)) {
      pal = c(cols[1:mod$nShifts],'#808080')
    } else {
      pal = c(brewer.pal(mod$nShifts, 'Spectral'),'#808080')
    }
    
    plot(mod, plot.bar = FALSE,  asterisk = FALSE, edge.label.ann = FALSE,
         show.tip.label = FALSE, palette = pal, 
         edge.width=4, edge.shift.adj = c(0.9, -0.02), edge.ann.cex=1.5,
         margin = FALSE)
  }
  
}
