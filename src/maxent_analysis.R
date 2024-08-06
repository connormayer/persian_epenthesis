library(tidyverse)
library(maxent.ot)
library(ggrepel)

options(dplyr.summarise.inform = FALSE)

##################################################
# Fit weights to data pooled across participants #
##################################################

# Not currently used, but this fits separate models to each individual
# participant. Provides an upper bound for model performance
fit_best_models <- function(full_model, pc_scores, name_template, folder) {
  total_ll <- 0
  predictions <- data.frame()
  # Go through tableaux for each participant
  for (p in pc_scores$participant) {
    filename <- file.path(folder, str_glue("{name_template}_p{p}.csv"))
    p_tableau <- read_csv(filename, show_col_types=FALSE)
    best_p_model <- optimize_weights(p_tableau)
    cur_preds <- predict_probabilities(p_tableau, best_p_model$weights) 
    total_ll <- total_ll + best_p_model$loglik
    predictions <- rbind(predictions, 
                         cur_preds$predictions %>%
                           mutate(participant=p))
  }
  return(list(
    name = str_glue('{name_template}_ind'),
    loglik=total_ll, 
    k = full_model$k * nrow(pc_scores),
    n = full_model$n,
    predictions = predictions,
    w = full_model$weights
  ))
}

fit_individual_models <- function(full_model, pc_scores, name_template, folder, separate_rho=FALSE) {
  # Rho is the amount by which we scale the Persian dominance score. Rho of 0
  # means no effect of dominance.
  rhos_to_test <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
  best_ll <- -Inf
  best_rho <- -1
  best_preds <- NULL
  
  # H = SUM_i (w_i + rho * pc_score) * constraint_violation_i
  
  if (!separate_rho) {
    # Try every possible rho value
    for (rho in rhos_to_test) {
      total_ll <- 0
      predictions <- data.frame()
      # Go through tableaux for each participant
      for (p in pc_scores$participant) {
          filename <- file.path(folder, str_glue("{name_template}_p{p}.csv"))
          p_tableau <- read_csv(filename, show_col_types=FALSE)
          
          # Get Persian dominance for participant
          scaling_factor <- pc_scores %>%
            filter(participant == p)
          
          # Scale weight of *Complex constraint(s) up or down depending on
          # scaling factor
          weights <- full_model$weights
          if ('*Complex' %in% names(weights)) {
            weights['*Complex'] <- max(weights['*Complex'] + rho * scaling_factor$PC1, 0)
          } else {
            weights['*Complex-S'] <- max(weights['*Complex-S'] + rho * scaling_factor$PC1, 0)
            weights['*Complex-T'] <- max(weights['*Complex-T'] + rho * scaling_factor$PC1, 0)
          }
          # Generate predictions and LL under scaled weights
          cur_preds <- predict_probabilities(p_tableau, weights) 
          total_ll <- total_ll + cur_preds$loglik
          predictions <- rbind(predictions, 
                               cur_preds$predictions %>%
                                  mutate(participant=p))
      }
      if (total_ll > best_ll) {
        # If overall LL is better for this value of rho, store it
        best_ll <- total_ll
        best_rho <- rho
        best_preds <- predictions
      }
    }
    return(list(
      name = str_glue('{name_template}_ind'),
      loglik=best_ll, 
      rho=best_rho,
      k = full_model$k + 1,
      n = full_model$n,
      predictions = best_preds,
      w = full_model$weights
    ))
  } else{
    # Try every possible rho value
    for (rho_s in rhos_to_test) {
      for (rho_t in rhos_to_test) {
        print(str_glue("Rho S: {rho_s}, Rho T: {rho_t}"))
        total_ll <- 0
        predictions <- data.frame()
        # Go through tableaux for each participant
        for (p in pc_scores$participant) {
          filename <- file.path(folder, str_glue("{name_template}_p{p}.csv"))
          p_tableau <- read_csv(filename, show_col_types=FALSE)
          
          # Get Persian dominance for participant
          scaling_factor <- pc_scores %>%
            filter(participant == p)
          
          # Scale weight of *Complex constraint(s) up or down depending on
          # scaling factor
          weights <- full_model$weights
          weights['*Complex-S'] <- max(weights['*Complex-S'] + rho_s * scaling_factor$PC1, 0)
          weights['*Complex-T'] <- max(weights['*Complex-T'] + rho_t * scaling_factor$PC1, 0)
          # Generate predictions and LL under scaled weights
          cur_preds <- predict_probabilities(p_tableau, weights) 
          total_ll <- total_ll + cur_preds$loglik
          predictions <- rbind(predictions, 
                               cur_preds$predictions %>%
                                 mutate(participant=p))
        }
        if (total_ll > best_ll) {
          # If overall LL is better for this value of rho, store it
          best_ll <- total_ll
          best_rho_s <- rho_s
          best_rho_t <- rho_t
          best_preds <- predictions
        }
      }
    }
    return(list(
      name = str_glue('{name_template}_ind_rho'),
      loglik=best_ll, 
      rho_s=best_rho_s,
      rho_t=best_rho_t,
      k = full_model$k + 2,
      n = full_model$n,
      predictions = best_preds,
      w = full_model$weights
    ))
  }
}

pc_scores <- read_csv('data/experiment/experimental_revised_results.csv') %>%
  group_by(participant) %>%
  summarize(PC1_acquisition_exposure = mean(PC1_acquisition_exposure),
            PC1_original = mean(PC1_original))

pc_orig <- pc_scores %>% 
  select(-PC1_acquisition_exposure) %>%
  mutate(PC1 = PC1_original)

pc_acq <- pc_scores %>% 
  select(-PC1_original) %>%
  mutate(PC1 = PC1_acquisition_exposure)

# Fleischhacker global 
fh_global <- read_csv("data/tableaux/fleischhacker_global.csv")
#cross_validate(fh_global, 5, 0, c(100, 5), grid_search=TRUE)
fh_model <- optimize_weights(fh_global, mu = 0, sigma = 100, model_name='fh')
# Fit individual models
fh_ind_model_acq <- fit_individual_models(
  fh_model, pc_acq, 'fh', 'data/tableaux/fleischhacker_ind'
)

# Gouskova simple global 
gs_global <- read_csv("data/tableaux/gouskova_simple_global.csv")
#cross_validate(gs_global, 5, 0, c(2000, 1000, 5, 2, 1), grid_search=TRUE)
gs_model <- optimize_weights(gs_global, mu = 0, sigma = 5, model_name='gs')
# Fit individual models
gs_ind_model_acq <- fit_individual_models(
  gs_model, pc_acq, 'gs', 'data/tableaux/gouskova_simple_ind'
)

# Gouskova_complex global 
gc_global <- read_csv("data/tableaux/gouskova_complex_global.csv")
#cross_validate(gc_global, 5, 0, c(1000, 500, 200, 100, 50, 20, 5, 2, 1), grid_search=TRUE)
gc_model <- optimize_weights(gc_global, mu = 0, sigma = 500, model_name='gc', upper_bound = 70)
gc_ind_model_acq <- fit_individual_models(
  gc_model, pc_acq, 'gc', 'data/tableaux/gouskova_complex_ind'
)

# Fleischhacker global split *complex
fh_global_split <- read_csv("data/tableaux/fleischhacker_global_split.csv")
#cross_validate(fh_global_split, 5, 0, c(1000, 500, 200, 100, 50, 20, 5, 2, 1), grid_search=TRUE)
fh_split_model <- optimize_weights(fh_global_split, mu = 0, sigma = 1000, model_name='fh_split') 
fh_split_ind_model_acq <- fit_individual_models(
  fh_split_model, pc_acq, 'fh_split', 'data/tableaux/fleischhacker_split_ind'
)
fh_split_ind_model_rho_acq <- fit_individual_models(
  fh_split_model, pc_acq, 'fh_split', 'data/tableaux/fleischhacker_split_ind',
  separate_rho=TRUE
)



# Gouskova simple global split *complex
gs_global_split <- read_csv("data/tableaux/gouskova_simple_global_split.csv")
#cross_validate(gs_global_split, 5, 0, c(1000, 500, 200, 100, 50, 20, 5, 2, 1), grid_search=TRUE)
gs_split_model <- optimize_weights(gs_global_split, mu = 0, sigma = 5, model_name='gs_split')
gs_split_ind_model_acq <- fit_individual_models(
  gs_split_model, pc_acq, 'gs_split', 'data/tableaux/gouskova_simple_split_ind'
)
gs_split_ind_model_rho_acq <- fit_individual_models(
  gs_split_model, pc_acq, 'gs_split', 'data/tableaux/gouskova_simple_split_ind',
  separate_rho=TRUE
)

# Gouskova complex global split *complex 
gc_global_split <- read_csv("data/tableaux/gouskova_complex_global_split.csv")
#cross_validate(gc_global_split, 5, 0, c(1000, 500, 200, 100, 50, 20, 5, 2, 1), grid_search=TRUE)
gc_split_model <- optimize_weights(gc_global_split, mu = 0, sigma = 200, model_name='gc_split', upper_bound = 70)
gc_split_ind_model_orig <- fit_individual_models(
  gc_split_model, pc_orig, 'gc_split_orig', 'data/tableaux/gouskova_complex_split_ind'
)
gc_split_ind_model_acq <- fit_individual_models(
  gc_split_model, pc_acq, 'gc_split', 'data/tableaux/gouskova_complex_split_ind'
)
gc_split_ind_model_rho_acq <- fit_individual_models(
  gc_split_model, pc_acq, 'gc_split', 'data/tableaux/gouskova_complex_split_ind',
  separate_rho=TRUE
)

# compare models
compare_models(
  fh_model, 
  fh_split_model,
  fh_ind_model_acq,
  fh_split_ind_model_acq,
  fh_split_ind_model_rho_acq,
  
  gs_model,
  gs_split_model,
  gs_ind_model_acq,
  gs_split_ind_model_acq,
  gs_split_ind_model_rho_acq,
  
  gc_model,
  gc_split_model,
  gc_ind_model_acq,
  gc_split_ind_model_acq,
  gc_split_ind_model_rho_acq,
  
  method = "bic"
)

#--------

# predict probabilities and save to fh/gs/gc_results
fh_results <- predict_probabilities(fh_global, fh_model$weights)
gs_results <- predict_probabilities(gs_global, gs_model$weights)
gc_results <- predict_probabilities(gc_global, gc_model$weights)


fh_split_results <- predict_probabilities(fh_global_split, fh_split_model$weights)
gs_split_results <- predict_probabilities(gs_global_split, gs_split_model$weights)
gc_split_results <- predict_probabilities(gc_global_split, gc_split_model$weights)


compare <- fh_split_results
compare$fh_gs_diff <- fh_split_results$predictions$Predicted - gs_split_results$predictions$Predicted
compare$fh_gc_diff <- fh_split_results$predictions$Predicted - gc_split_results$predictions$Predicted
compare$gs_gc_diff <- gs_split_results$predictions$Predicted - gc_split_results$predictions$Predicted

#################
# VISUALIZATION #
#################
plot_model <- function(model_preds, onsets, outname) {
  outdir <- file.path('figures', outname)
  dir.create(outdir, showWarnings = FALSE)
  
  # Join the predicted results with the onset column
  preds_with_onsets <- inner_join(model_preds$predictions, onsets, by=c('Input')) %>% 
    mutate(ep_type = str_split_fixed(Output, '-', 2)[,2]) %>%
    filter(!is.na(Observed))
  
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

plot_model(fh_results, onsets, 'fh')
plot_model(gs_results, onsets, 'gs')
plot_model(gc_results, onsets, 'gc')
plot_model(fh_split_results, onsets, 'fh_split')
plot_model(gs_split_results, onsets, 'gs_split')
plot_model(gc_split_results, onsets, 'gc_split')
plot_model(fh_ind_model, onsets, 'fh')
plot_model(gs_ind_model, onsets, 'gs')
plot_model(gc_ind_model, onsets, 'gc')
plot_model(fh_split_ind_model, onsets, 'fh_split_ind')
plot_model(gs_split_ind_model, onsets, 'gs_split_ind')
plot_model(gc_split_ind_model, onsets, 'gc_split_ind')

# inner join experimental results and *_results, then group by sonority
# delta (or whatever) and get the mean error

#load experimental results 
experimental_results <- read_csv('data/experiment/experimental_results.csv') 

#create errors function
sonority_errors <- function(results) {
  #join df 
  # In order for this to work, experimental_results already needs to be defined
  sonority_scores <- experimental_results %>%
    select(word, nap_sonority, sonority_delta) %>%
    distinct()
  
  joined_df <- inner_join(results$predictions, sonority_scores, by = c("Input" = "word"))
  
  #mean errors 
  joined_df <- joined_df %>%
    group_by(sonority_delta) %>%
    mutate(mean_delta_error = mean(Error^2))
  
  joined_df <- joined_df %>%
    group_by(nap_sonority) %>%
    mutate(mean_nap_error = mean(Error^2)) }

#run on results 
gs_sonority <- sonority_errors(gs_results)
gc_sonority <- sonority_errors(gc_results)
fh_sonority <- sonority_errors(fh_results)
gs_split_results <- sonority_errors(gs_split_results)
gc_split_results <- sonority_errors(gc_split_results)
fh_split_results <- sonority_errors(fh_split_results)

fh_best <- fh_split_ind_model_rho$predictions %>%
  inner_join(onsets, by=c("Input")) %>%
  separate(Output, c('word', 'ep_type'), sep='-') %>%
  mutate(s_initial = str_starts(onset, 's')) %>%
  filter(!is.na(Error))

fh_best_onset <- fh_best %>%
  group_by(onset, ep_type) %>%
  summarize(mean_error = mean(Error))

fh_best_s <- fh_best %>%
  group_by(s_initial, ep_type) %>%
  summarize(mean_error = mean(Error))

fh_best_onset %>%
  ggplot() +
  geom_bar(aes(x=fct_reorder(ep_type, mean_error), y=mean_error), stat='identity') +
  facet_wrap(~ onset)

fh_best_s %>% 
  mutate(s_initial = ifelse(s_initial, 'sC', 'TR')) %>%
  ggplot() +
  geom_bar(aes(x=fct_reorder(ep_type, mean_error), y=mean_error, fill=ep_type), stat='identity') +
  facet_wrap(~ s_initial) +
  xlab("Cluster type") +
  ylab("Mean error") + 
  scale_fill_discrete(guide="none") +
  ggtitle("Perceptual model") +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16, angle = 45, vjust = 0.5, hjust = 1),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  ylim(-0.15, 0.15)
ggsave('figures/fh_best.png', height = 7, width = 7, units='in')

gs_best <- gs_split_ind_model_rho$predictions %>%
  inner_join(onsets, by=c("Input")) %>%
  separate(Output, c('word', 'ep_type'), sep='-') %>%
  mutate(s_initial = str_starts(onset, 's')) %>%
  filter(!is.na(Error))

gs_best_onset <- gs_best %>%
  group_by(onset, ep_type) %>%
  summarize(mean_error = mean(Error))

gs_best_s <- gs_best %>%
  group_by(s_initial, ep_type) %>%
  summarize(mean_error = mean(Error))

gs_best_onset %>%
  ggplot() +
  geom_bar(aes(x=fct_reorder(ep_type, mean_error), y=mean_error), stat='identity') +
  facet_wrap(~ onset)

gs_best_s %>%
  mutate(s_initial = ifelse(s_initial, 'sC', 'TR')) %>%
  ggplot() +
  geom_bar(aes(x=fct_reorder(ep_type, mean_error), y=mean_error, fill=ep_type), stat='identity') +
  facet_wrap(~ s_initial) +
  xlab("Cluster type") +
  ylab("Mean error") + 
  scale_fill_discrete(guide="none") +
  ggtitle("Syl. Cont. Simple Model") +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16, angle = 45, vjust = 0.5, hjust = 1),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) + 
  ylim(-0.15, 0.15)
ggsave('figures/gs_best.png', height = 7, width = 7, units='in')

gc_best <- gc_split_ind_model_rho$predictions %>%
  inner_join(onsets, by=c("Input")) %>%
  separate(Output, c('word', 'ep_type'), sep='-') %>%
  mutate(s_initial = str_starts(onset, 's')) %>%
  filter(!is.na(Error))

gc_best_onset <- gc_best %>%
  group_by(onset, ep_type) %>%
  summarize(mean_error = mean(Error))

gc_best_s <- gc_best %>%
  group_by(s_initial, ep_type) %>%
  summarize(mean_error = mean(Error))

gc_best_onset %>%
  ggplot() +
  geom_bar(aes(x=fct_reorder(ep_type, mean_error), y=mean_error), stat='identity') +
  facet_wrap(~ onset)

gc_best_s %>%
  mutate(s_initial = ifelse(s_initial, 'sC', 'TR')) %>%
  ggplot() +
  geom_bar(aes(x=fct_reorder(ep_type, mean_error), y=mean_error, fill=ep_type), stat='identity') +
  facet_wrap(~ s_initial) + 
  xlab("Cluster type") +
  ylab("Mean error") + 
  scale_fill_discrete(guide="none") +
  ggtitle("Syl. Cont. Hierarchical Model") +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16, angle = 45, vjust = 0.5, hjust = 1),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  ylim(-0.15, 0.15)
ggsave('figures/gc_best.png', height = 7, width = 7, units='in')
