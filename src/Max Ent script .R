library(tidyverse)
library(maxent.ot)

#set wd
#setwd("~/Desktop/git_repo_persian_epenthesis")
setwd("E:/git_repos/persian_epenthesis")
#setwd("C:/Users/conno/git_repos/persian_epenthesis")

##################################################
# Fit weights to data pooled across participants #
##################################################

# Fleischhacker global 
fh_global <- read_csv("data/tableaux/fh_global.csv")

# fit model 
fh_model <- optimize_weights(fh_global)
fh_model$weights
fh_model$loglik

# Look at predicted frequencies
fh_results <- predict_probabilities(fh_global, fh_model$weights)
fh_results$predictions

# Gouskova simple global 
gs_global <- read_csv("data/tableaux/gs_global.csv")

#fit model 
gs_model <- optimize_weights(gs_global)
gs_model$weights
gs_model$loglik

gs_results <- predict_probabilities(gs_global, gs_model$weights)
gs_results$predictions

#gouskova_complex global 
gc_global <- read_csv("data/tableaux/gc_global.csv")

#fit model 
gc_model <- optimize_weights(gc_global)
gc_model$weights
gc_model$loglik

gc_results <- predict_probabilities(gc_global, gc_model$weights)
gc_results$predictions

#compare models 
compare_models(fh_model, gs_model, gc_model, method = "bic")

#################################################
# Look at models fit to individual participants #
#################################################

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
  # Convert to numeric
  results_df <- as_tibble(type.convert(results_df))
  # Write results
  write_csv(results_df, output)
  return(results_df)
}

fh_results_df <- fit_models(
  "data/tableaux/fleischhacker", 
  "data/individual_results_fleischhacker.csv"
)

gs_results_df <- fit_models(
  "data/tableaux/gouskova_simple", 
  "data/individual_results_gouskova_simple.csv"
)

gc_results_df <- fit_models(
  "data/tableaux/gouskova_complex", 
  "data/individual_results_gouskova_complex.csv"
)

experimental_df <- read_csv('data/experimental_results.csv') %>%
  select(participant, PC1) %>%
  unique()

process_data <- function(model_df, experimental_df, prefix, constraint_names) {
  regex <- str_glue("{prefix}_p(\\d+)\\.csv")
  new_df <- model_df
  new_df$participant <- as.double(str_match(new_df$participant, regex)[,2])
  joined_df <- inner_join(new_df, experimental_df, by=c("participant"))
  joined_df$participant <- as_factor(joined_df$participant)
  pivot_df <- joined_df %>% 
    pivot_longer(constraint_names, names_to = "Constraint", values_to = "Weight")
  normalized_df <- pivot_df %>% 
    group_by(participant) %>% 
    mutate(weights_normalized = as.numeric(Weight) / max(Weight)) 
  final_df <- normalized_df 
  final_df$participant <- factor(normalized_df$participant)   
  final_df$participant <- fct_reorder(normalized_df$participant, normalized_df$PC1)
  return(final_df)
}

plot_results <- function(df) {
  ggplot(df, aes(x=PC1, y=weights_normalized)) +
    geom_point() +
    geom_smooth(method='lm') +
    facet_wrap(~Constraint) +
    xlab("Relative Farsi Dominance") + 
    ylab("Normalized constraint weight")
}

fh_prefix <- 'fh'
fh_colnames <- colnames(fh_results_df)
fh_names <- fh_colnames[2:(length(fh_colnames)-3)]
fh_final <- process_data(fh_results_df, experimental_df, fh_prefix, fh_names)
plot_results(fh_final)

gs_prefix <- 'gs'
gs_colnames <- colnames(gs_results_df)
gs_names <- gs_colnames[2:(length(gs_colnames)-3)]
gs_final <- process_data(gs_results_df, experimental_df, gs_prefix, gs_names)
plot_results(gs_final)

gc_prefix <- 'gc'
gc_colnames <- colnames(gc_results_df)
gc_names <- gc_colnames[2:(length(gc_colnames)-3)]
gc_final <- process_data(gc_results_df, experimental_df, gc_prefix, gc_names)
plot_results(gc_final)
