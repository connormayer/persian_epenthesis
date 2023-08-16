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


#run max ent on each individual participant
#fleischhacker
library(fs)

fh_file_path <- fs::dir_ls("tableaux/fleischhacker")
fh_file_path

fh_file_contents <- list()
for(i in seq_along(fh_file_path)) {
  fh_file_contents[[i]] <- read_csv(
    file = fh_file_path[[i]]
  )
}

fh_file_contents <- set_names(fh_file_contents, c("fh_global.csv", "fh_p1.csv", "fh_p2.csv", "fh_p3.csv", "fh_p4.csv", "fh_p5.csv", "fh_p6.csv", "fh_p7.csv", "fh_p8.csv", "fh_p9.csv", "fh_p10.csv", "fh_p13.csv", "fh_p14.csv", "fh_p15.csv", "fh_p16.csv", "fh_p17.csv", "fh_p18.csv", "fh_p19.csv", "fh_p20.csv", "fh_p21.csv","fh_template.csv"))


lapply(fh_file_contents, optimize_weights)


