library(tidyverse)
library(maxent.ot)
library(ggrepel)

##################################################
# Fit weights to data pooled across participants #
##################################################

# Fleischhacker global 
fh_global <- read_csv("data/tableaux/fleischhacker_global.csv")
fh_model <- optimize_weights(fh_global)

# Gouskova simple global 
gs_global <- read_csv("data/tableaux/gouskova_simple_global.csv")
gs_model <- optimize_weights(gs_global)

# Gouskova_complex global 
gc_global <- read_csv("data/tableaux/gouskova_complex_global.csv")
gc_model <- optimize_weights(gc_global)

# Fleischhacker global split *complex
fh_global_split <- read_csv("data/tableaux/fleischhacker_global_split.csv")
fh_split_model <- optimize_weights(fh_global_split) 

# Gouskova simple global split *complex
gs_global_split <- read_csv("data/tableaux/gouskova_simple_global_split.csv")
gs_split_model <- optimize_weights(gs_global_split)

# Gouskova complex global split *complex 
gc_global_split <- read_csv("data/tableaux/gouskova_complex_global_split.csv")
gc_split_model <- optimize_weights(gc_global_split)

# compare models
compare_models(fh_model, gs_model, gc_model, fh_split_model, gs_split_model, gc_split_model, method = "bic")

# predict probabilities and save to fh/gs/gc_results
fh_results <- predict_probabilities(fh_global, fh_model$weights)
gs_results <- predict_probabilities(gs_global, gs_model$weights)
gc_results <- predict_probabilities(gc_global, gc_model$weights)
fh_split_results <- predict_probabilities(fh_global_split, fh_split_model$weights)
gs_split_results <- predict_probabilities(gs_global_split, gs_split_model$weights)
gc_split_results <- predict_probabilities(gc_global_split, gc_split_model$weights)

#################
# VISUALIZATION #
#################
plot_model <- function(model_preds, onsets, outname) {
  outdir <- file.path('figures', outname)
  dir.create(outdir, showWarnings = FALSE)
  
  # Join the predicted results with the onset column
  preds_with_onsets <- inner_join(model_preds$predictions, onsets, by=c('Input')) %>% 
    mutate(ep_type = str_split_fixed(Output, '-', 2)[,2])
  
  # Error rates by onset types
  preds_with_onsets %>%
    group_by(onset) %>% 
    mutate(total_error=sum(abs(Error))) %>%
    select(onset, total_error) %>%
    unique() %>%
    arrange(-total_error) %>% 
    ggplot(aes(y=total_error, x=fct_reorder(onset, total_error))) +
    geom_bar(stat='identity')
  
  ggsave(file.path(outdir, str_c(outname, '_error_by_onset.png')))
  
  # Aggregate predicted vs. observed frequencies for each onset type and make
  # a plot showing them
  obs_pred_by_onset_ep_type <- preds_with_onsets %>%
    mutate(onset_ep=str_c(onset, '-', str_sub(ep_type, 1, 1))) %>%
    group_by(onset_ep) %>% 
    mutate(pred=mean(Predicted),
           obs=mean(Observed),
           error=mean(Error)) %>%
    select(onset_ep, obs, pred, error) %>%
    unique() 
  
  ggplot(obs_pred_by_onset_ep_type, aes(x=obs, y=pred, label=onset_ep)) +
    geom_point(size=2) +
    geom_label_repel(max.overlaps = Inf) +
    geom_abline(intercept=0, slope=1, size=1) + 
    theme(axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.x = element_text(size=25),
          axis.title.y = element_text(size=25)) +
    xlim(0,1) +
    ylim(0,1)
  ggsave(file.path(outdir, str_c(outname, '_pred_obs_full.png')))
  
  # Plot predicted vs. observed frequencies for infrequent outcomes
  # I just did this to make it a bit easier to read, this is just the left
  # half of the above graph
  ggplot(obs_pred_by_onset_ep_type %>%
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
  ggsave(file.path(outdir, str_c(outname, '_pred_obs_lower.png')))
  
  # Same thing but for the right half of the graph
  ggplot(obs_pred_by_onset_ep_type %>%
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
  ggsave(file.path(outdir, str_c(outname, '_pred_obs_upper.png')))
}

# Load experimental results to get pairs of words and onsets
onsets <- read_csv('data/experiment/experimental_results.csv') %>%
  select(word, onset) %>%
  mutate(Input=word) %>%
  select(-word) %>%
  unique() %>%
  filter(Input != 'spreading') 

plot_model(fh_results, onsets, 'fh_global')
plot_model(gs_results, onsets, 'gs_global')
plot_model(gc_results, onsets, 'gc_global')
plot_model(fh_split_results, onsets, 'fh_split_global')
plot_model(gs_split_results, onsets, 'gs_split_global')
plot_model(gc_split_results, onsets, 'gc_split_global')

# inner join experimental results and *_results, then group by sonority
# delta (or whatever) and get the mean error

#load experimental results 
experimental_results <- read_csv('data/experiment/experimental_results.csv') 

#create errors function
sonority_errors <- function(results) {
  #join df 
  joined_df <- full_join(results$predictions, experimental_results, by = c("Input" = "word")) %>%
    distinct(Input, Output, .keep_all = TRUE) %>%
    select(Input, Output, Predicted, Observed, Error, nap_sonority, sonority_delta)
  
  #mean errors 
  joined_df <- joined_df%>%
    group_by(sonority_delta) %>%
    mutate(delta_error = mean(Error))
  
  joined_df <- joined_df %>%
    group_by(nap_sonority) %>%
    mutate(mean_nap_error = mean(Error)) }

#run on results 
gs_sonority <- sonority_errors(gs_results)
gc_sonority <- sonority_errors(gc_results)
fh_sonority <- sonority_errors(fh_results)
gs_split_results <- sonority_errors(gs_split_results)
gc_split_results <- sonority_errors(gc_split_results)
fh_split_results <- sonority_errors(fh_split_results)
