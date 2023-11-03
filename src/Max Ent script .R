library(tidyverse)
library(maxent.ot)
install.packages("ggrepel")
library(ggrepel)

#set wd
#setwd("~/Desktop/git_repo_persian_epenthesis")
setwd("E:/git_repos/persian_epenthesis")
#setwd("C:/Users/conno/git_repos/persian_epenthesis")

##################################################
# Fit weights to data pooled across participants #
##################################################

# Fleischhacker global 
fh_global <- read_csv("tableaux/fh_global.csv")

# fit model 
fh_model <- optimize_weights(fh_global)
fh_model$weights
fh_model$loglik

# Look at predicted frequencies
fh_results <- predict_probabilities(fh_global, fh_model$weights)
fh_results$predictions

# Gouskova simple global 
gs_global <- read_csv("tableaux/gs_global.csv")

#fit model 
gs_model <- optimize_weights(gs_global)
gs_model$weights
gs_model$loglik

gs_results <- predict_probabilities(gs_global, gs_model$weights)
gs_results$predictions

#gouskova_complex global 
gc_global <- read_csv("tableaux/gc_global.csv")

#fit model 
gc_model <- optimize_weights(gc_global)
gc_model$weights
gc_model$loglik

gc_results <- predict_probabilities(gc_global, gc_model$weights)
gc_results$predictions

#compare models without *complex split 
compare_models(fh_model, gs_model, gc_model, method = "bic")



#Models with split up *complex 
#Fleischhacker 
fh_star_complex <- read_csv("fh_star_complex.csv")
fh_complex_model <- optimize_weights(fh_star_complex) 
fh_complex_model$weights

#Gouskova Complex
gc_star_complex <- read_csv("gc_star_complex.csv")
gc_complex_model <- optimize_weights(gc_star_complex)
gc_complex_model$loglik
gc_complex_model$weights


#Gouskova simple 
gs_star_complex <- read_csv("gs_star_complex.csv")
gs_complex_model <- optimize_weights(gs_star_complex)
gs_complex_model$weights

#compare models with star complex split
compare_models(fh_model, gs_model, gc_model,fh_complex_model, gc_complex_model, gs_complex_model, method = "bic")

#predict probabilities and save to fh/gs/gc_results
fh_complex_results <- predict_probabilities(fh_star_complex, fh_complex_model$weights)
gs_complex_results <- predict_probabilities(gs_star_complex, gs_complex_model$weights)
gc_complex_results <- predict_probabilities(gc_star_complex, gc_complex_model$weights)


#visualization 

# Load experimental results to get pairs of words and onsets
onsets <- read_csv('experimental_results.csv') %>%
  select(word, onset) %>%
  mutate(Input=word) %>%
  select(-word) %>%
  unique() %>%
  filter(Input != 'spreading') 

# Join the predicted results with the onset column
fh_preds_with_onsets <- inner_join(fh_complex_results$predictions, onsets, by=c('Input')) %>% 
  unique() %>% 
  mutate(ep_type = rep(c("none", "anaptyxis", "prothesis"), times=71))


view(fh_preds_with_onsets)
view(onsets)
# Aggregate error rates by onset type and make a plot showing them
fh_error_by_onset <- fh_preds_with_onsets %>%
  group_by(onset) %>% 
  mutate(total_error=sum(abs(Error))) %>%
  select(onset, total_error) %>%
  unique() %>%
  arrange(-total_error)

ggplot(fh_error_by_onset, aes(y=total_error, x=fct_reorder(onset, total_error))) +
  geom_bar(stat='identity')

# Aggregate predicted vs. observed frequencies for each onset type and make
# a plot showing them
fh_obs_pred_by_onset_ep_type <- fh_preds_with_onsets %>%
  mutate(onset_ep=str_c(onset, '-', str_sub(ep_type, 1, 1))) %>%
  group_by(onset_ep) %>% 
  mutate(pred=mean(Predicted),
         obs=mean(Observed),
         error=mean(Error)) %>%
  select(onset_ep, obs, pred, error) %>%
  unique() 

ggplot(fh_obs_pred_by_onset_ep_type, aes(x=obs, y=pred, label=onset_ep)) +
  geom_point(size=2) +
  geom_label_repel(max.overlaps = Inf) +
  geom_abline(intercept=0, slope=1, size=1) + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25)) +
  xlim(0,1) +
  ylim(0,1)

# Prot predicted vs. observed frequencies for infrequent outcomes
# I just did this to make it a bit easier to read, this is just the left
# half of the above graph
ggplot(fh_obs_pred_by_onset_ep_type %>%
         filter(obs < 0.5), aes(x=obs, y=pred, label=onset_ep)) +
  geom_point(size=2) +
  geom_label_repel(max.overlaps = Inf) +
  geom_abline(intercept=0, slope=1, size=1) + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25)) +
  xlim(0,0.3) +
  ylim(0,0.3)

# Same thing but for the right half of the graph
ggplot(fh_obs_pred_by_onset_ep_type %>%
         filter(obs >= 0.5), aes(x=obs, y=pred, label=onset_ep)) +
  geom_point(size=2) +
  geom_label_repel(max.overlaps = Inf) +
  geom_abline(intercept=0, slope=1, size=1) + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25)) +
  xlim(0.7, 1) +
  ylim(0.7, 1)

# You'll have to do the same thing for the other data sets. I would actually
# see if you can turn the code above into a function where you just pass in the
# data and you save the appropriate plots. We could also go through how to do 
# this next week







# #################################################
# # Look at models fit to individual participants #
# #################################################
# 
# fit_models <- function(tableaux_folder, output) {
#   # Make empty tibble to hold model results
#   results_df <- tibble()
#   
#   # Loop over files in folder
#   for (f in list.files(tableaux_folder)) {
#     # Concatenate folder path and filename to get full path to file
#     full_path <- file.path(tableaux_folder, f)
#     # Read file
#     tableau_df <- read_csv(full_path)
#     # Optimize weights
#     model <- optimize_weights(tableau_df)
#     # Store results in dataframe
#     results_df <- rbind(results_df, c(f, model$weights, model$loglik, model$k, model$n))
#   }
#   # Add column names
#   colnames(results_df) <- c('participant', names(model$weights), 'loglik', 'k', 'n')
#   # Convert to numeric
#   results_df <- as_tibble(type.convert(results_df))
#   # Write results
#   write_csv(results_df, output)
#   return(results_df)
# }
# 
# fh_results_df <- fit_models(
#   "data/tableaux/fleischhacker", 
#   "data/individual_results_fleischhacker.csv"
# )
# 
# gs_results_df <- fit_models(
#   "data/tableaux/gouskova_simple", 
#   "data/individual_results_gouskova_simple.csv"
# )
# 
# gc_results_df <- fit_models(
#   "data/tableaux/gouskova_complex", 
#   "data/individual_results_gouskova_complex.csv"
# )
# 
# experimental_df <- read_csv('data/experimental_results.csv') %>%
#   select(participant, PC1) %>%
#   unique()
# 
# process_data <- function(model_df, experimental_df, prefix, constraint_names) {
#   regex <- str_glue("{prefix}_p(\\d+)\\.csv")
#   new_df <- model_df
#   new_df$participant <- as.double(str_match(new_df$participant, regex)[,2])
#   joined_df <- inner_join(new_df, experimental_df, by=c("participant"))
#   joined_df$participant <- as_factor(joined_df$participant)
#   pivot_df <- joined_df %>% 
#     pivot_longer(constraint_names, names_to = "Constraint", values_to = "Weight")
#   normalized_df <- pivot_df %>% 
#     group_by(participant) %>% 
#     mutate(weights_normalized = as.numeric(Weight) / max(Weight)) 
#   final_df <- normalized_df 
#   final_df$participant <- factor(normalized_df$participant)   
#   final_df$participant <- fct_reorder(normalized_df$participant, normalized_df$PC1)
#   return(final_df)
# }
# 
# plot_results <- function(df) {
#   ggplot(df, aes(x=PC1, y=weights_normalized)) +
#     geom_point() +
#     geom_smooth(method='lm') +
#     facet_wrap(~Constraint) +
#     xlab("Relative Farsi Dominance") + 
#     ylab("Normalized constraint weight")
# }
# 
# fh_prefix <- 'fh'
# fh_colnames <- colnames(fh_results_df)
# fh_names <- fh_colnames[2:(length(fh_colnames)-3)]
# fh_final <- process_data(fh_results_df, experimental_df, fh_prefix, fh_names)
# plot_results(fh_final)
# fh_ll <- sum((fh_final %>% select(participant, loglik) %>% unique())$loglik)
# 
# gs_prefix <- 'gs'
# gs_colnames <- colnames(gs_results_df)
# gs_names <- gs_colnames[2:(length(gs_colnames)-3)]
# gs_final <- process_data(gs_results_df, experimental_df, gs_prefix, gs_names)
# plot_results(gs_final)
# gs_ll <- sum((gs_final %>% select(participant, loglik) %>% unique())$loglik)
# 
# gc_prefix <- 'gc'
# gc_colnames <- colnames(gc_results_df)
# gc_names <- gc_colnames[2:(length(gc_colnames)-3)]
# gc_final <- process_data(gc_results_df, experimental_df, gc_prefix, gc_names)
# plot_results(gc_final)
# gc_ll <- sum((gc_final %>% select(participant, loglik) %>% unique())$loglik)
