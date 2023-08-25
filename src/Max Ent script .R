library(tidyverse)
library(maxent.ot)

#set wd
#setwd("~/Desktop/git_repo_persian_epenthesis")
setwd("E:/git_repos/persian_epenthesis")

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



fit_models <- function(tableaux_folder, output) {
  # Make empty tibble to hold model results
  results_df <- tibble()
  
  # Loop over files in folder
  for (f in list.files(tableaux_folder)) {
    # Concatenate folder path and filename to get full path to file
    full_path <- file.path(tableaux_folder, f)
    # Read file
    tableau_df <- read_csv(full_path)
    # Optimize weights
    model <- optimize_weights(tableau_df)
    # Store results in dataframe
    results_df <- rbind(results_df, c(f, model$weights, model$loglik, model$k, model$n))
  }
  # Add column names
  colnames(results_df) <- c('participant', names(model$weights), 'loglik', 'k', 'n')
  # Write results
  write_csv(results_df, output)
  return(results_df)
}

fh_results_df <- fit_models(
  "data/tableaux/fleischhacker", 
  "data/individual_results_fleischhacker.csv"
)

gouskova_results_df <- fit_models(
  "data/tableaux/gouskova_simple", 
  "data/individual_results_gouskova_simple.csv"
)

gouskova_complex_results_df <- fit_models(
  "data/tableaux/gouskova_complex", 
  "data/individual_results_gouskova_complex.csv"
)

# Add profiency information
experimental_df <- read_csv('data/experimental_results.csv') %>%
  select(participant, PC1) %>%
  unique()

regex <- "fh_p(\\d+)\\.csv"
fh_results_df$participant <- as.double(str_match(fh_results_df$participant, regex)[,2])
joined_df <- inner_join(fh_results_df, experimental_df, by=c("participant"))
view(joined_df)

#table4a %>% 
##pivot_longer(c(`1999`, `2000`), names_to = "year", values_to = "cases")

#pivot joined_df so that it has the columns Participant, Constraint, Weight, k, n, loglik, PC1.
fh_pivot <- joined_df %>%
  pivot_longer(c('DEP-[ə]/S_T', 'C//V', 'C/V', 'L-ANCHOR', 'CONTIGUITY', 'DEP-[ə]/S_N', 'DEP-[ə]/S_L', 'DEP-[ə]/S_W', 'DEP-[ə]/T_R', 'COMPLEX'),
            names_to = "Constraint", values_to = "Weight")
fh_pivot

class(fh_pivot$Weight)

#normalize constraint weights

fh_normalized <- fh_pivot %>% 
  group_by(participant) %>% 
  mutate(Max_weight = max(Weight), 
         weights_normalized = as.numeric(Weight) / as.numeric(Max_weight)) 


#run max ent on each individual participant
#fleischhacker
# library(fs)
# 
# fh_file_path <- fs::dir_ls("tableaux/fleischhacker")
# fh_file_path
# 
# fh_file_contents <- list()
# for(i in seq_along(fh_file_path)) {
#   fh_file_contents[[i]] <- read_csv(
#     file = fh_file_path[[i]]
#   )
# }
# 
# fh_file_contents <- set_names(fh_file_contents, c("fh_global.csv", "fh_p1.csv", "fh_p2.csv", "fh_p3.csv", "fh_p4.csv", "fh_p5.csv", "fh_p6.csv", "fh_p7.csv", "fh_p8.csv", "fh_p9.csv", "fh_p10.csv", "fh_p13.csv", "fh_p14.csv", "fh_p15.csv", "fh_p16.csv", "fh_p17.csv", "fh_p18.csv", "fh_p19.csv", "fh_p20.csv", "fh_p21.csv","fh_template.csv"))
# 
# 
# lapply(fh_file_contents, optimize_weights)

#run max ent on each individual participant
#fleischhacker
#

