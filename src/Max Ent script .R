library(tidyverse)
library(maxent.ot)

#set wd
setwd("~/Desktop/git_repo_persian_epenthesis")

#fleischhaker global 
fh_global <- read_csv("data/tableaux/fleischhacker/fh_global.csv")

#fit model 
fh_model <- optimize_weights(fh_global)
fh_model$weights
fh_model$loglik


results <- predict_probabilities(fh_global, fh_model$weights)
results$loglik
results$predictions



#gouskova_simple global 
gs_global <- read_csv("data/tableaux/gouskova_simple/gs_global.csv")

#fit model 
gs_model <- optimize_weights(gs_global)
gs_model$weights
gs_model$loglik


results <- predict_probabilities(gs_global, gs_model$weights)
results$loglik
results$predictions

#gouskova_complex global 
gc_global <- read_csv("data/tableaux/gouskova_complex/gc_global.csv")

#fit model 
gc_model <- optimize_weights(gc_global)
gc_model$weights
gc_model$loglik


results <- predict_probabilities(gc_global, gc_model$weights)
results$loglik
results$predictions

#compare models 
compare_models(fh_model, gs_model, gc_model, method = "bic")

