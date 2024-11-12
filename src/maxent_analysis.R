library(tidyverse)
library(maxent.ot)
library(ggrepel)
library(gridExtra)
library(ggpubr)

options(dplyr.summarise.inform = FALSE)

##################################################
# Fit weights to data pooled across participants #
##################################################

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

pc_acq <- pc_scores %>% 
  select(-PC1_original) %>%
  mutate(PC1 = PC1_acquisition_exposure)

# Fleischhacker global 
fh_global <- read_csv("data/tableaux/fleischhacker_global.csv")
fh_model <- optimize_weights(fh_global, mu = 0, sigma = 100, model_name='fh')
# Fit individual models
fh_ind_model_acq <- fit_individual_models(
  fh_model, pc_acq, 'fh', 'data/tableaux/fleischhacker_ind'
)

# Gouskova simple global 
gs_global <- read_csv("data/tableaux/gouskova_simple_global.csv")
gs_model <- optimize_weights(gs_global, mu = 0, sigma = 5, model_name='gs')
# Fit individual models
gs_ind_model_acq <- fit_individual_models(
  gs_model, pc_acq, 'gs', 'data/tableaux/gouskova_simple_ind'
)

# Gouskova_complex global 
gc_global <- read_csv("data/tableaux/gouskova_complex_global.csv")
gc_model <- optimize_weights(gc_global, mu = 0, sigma = 500, model_name='gc', upper_bound = 70)
gc_ind_model_acq <- fit_individual_models(
  gc_model, pc_acq, 'gc', 'data/tableaux/gouskova_complex_ind'
)

# Fleischhacker global split *complex
fh_global_split <- read_csv("data/tableaux/fleischhacker_global_split.csv")
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
gc_split_model <- optimize_weights(gc_global_split, mu = 0, sigma = 200, model_name='gc_split', upper_bound = 70)
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

onsets <- read_csv('data/experiment/experimental_raw_results.csv') %>%
  select(word, onset) %>%
  mutate(Input=word) %>%
  select(-word) %>%
  unique() %>%
  filter(Input != 'spreading') 

fh_best <- fh_split_ind_model_rho_acq$predictions %>%
  inner_join(onsets, by=c("Input")) %>%
  separate(Output, c('word', 'ep_type'), sep='-') %>%
  mutate(s_initial = str_starts(onset, 's')) %>%
  filter(!is.na(Error)) %>% 
  mutate(ep_type = plyr::revalue(ep_type, 
                          c("none" = "none", 
                            "anaptyxis" = 'medial\nepenthesis',
                            "prothesis" = 'pre-\nepenthesis')),
         ep_type = fct_relevel(ep_type, "none", 'medial\nepenthesis', "pre-\nepenthesis"))

fh_best_s <- fh_best %>%
  group_by(s_initial, ep_type) %>%
  summarize(mean_error = mean(Error))

fh_plot <- fh_best_s %>% 
  mutate(s_initial = ifelse(s_initial, 'sC', 'TR')) %>%
  ggplot() +
  geom_bar(aes(x=ep_type, y=mean_error, fill=ep_type), stat='identity') +
  facet_wrap(~ s_initial) +
  xlab("Cluster type") +
  ylab("Mean error") + 
  scale_fill_discrete(guide="none") +
  ggtitle("\nPerceptual\nCost") +
  theme_classic(base_size=20) +
  theme(axis.text=element_text(size=16, angle = 45, vjust = 0.5, hjust = 1),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  ylim(-0.15, 0.15) +
  geom_hline(yintercept=0, size=0.75)
fh_plot
ggsave('figures/fh_best.png', height = 7, width = 7, units='in')

gs_best <- gs_split_ind_model_rho_acq$predictions %>%
  inner_join(onsets, by=c("Input")) %>%
  separate(Output, c('word', 'ep_type'), sep='-') %>%
  mutate(s_initial = str_starts(onset, 's')) %>%
  filter(!is.na(Error)) %>% 
  mutate(ep_type = plyr::revalue(ep_type, 
                                 c("none" = "none", 
                                   "anaptyxis" = 'medial\nepenthesis',
                                   "prothesis" = 'pre-\nepenthesis')),
         ep_type = fct_relevel(ep_type, "none", 'medial\nepenthesis', "pre-\nepenthesis"))

gs_best_s <- gs_best %>%
  group_by(s_initial, ep_type) %>%
  summarize(mean_error = mean(Error))

gs_plot <- gs_best_s %>%
  mutate(s_initial = ifelse(s_initial, 'sC', 'TR')) %>%
  ggplot() +
  geom_bar(aes(x=ep_type, y=mean_error, fill=ep_type), stat='identity') +
  facet_wrap(~ s_initial) +
  xlab("Cluster type") +
  ylab("Mean error") + 
  scale_fill_discrete(guide="none") +
  ggtitle("Syllable\nContact\nSimple") +
  theme_classic(base_size=20) +
  theme(axis.text=element_text(size=16, angle = 45, vjust = 0.5, hjust = 1),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) + 
  ylim(-0.15, 0.15) +
  geom_hline(yintercept=0, size=0.75)
gs_plot
ggsave('figures/gs_best.png', height = 7, width = 7, units='in')

gc_best <- gc_split_ind_model_rho_acq$predictions %>%
  inner_join(onsets, by=c("Input")) %>%
  separate(Output, c('word', 'ep_type'), sep='-') %>%
  mutate(s_initial = str_starts(onset, 's')) %>%
  filter(!is.na(Error)) %>% 
  mutate(ep_type = plyr::revalue(ep_type, 
                                 c("none" = "none", 
                                   "anaptyxis" = 'medial\nepenthesis',
                                   "prothesis" = 'pre-\nepenthesis')),
         ep_type = fct_relevel(ep_type, "none", 'medial\nepenthesis', "pre-\nepenthesis"))

gc_best_s <- gc_best %>%
  group_by(s_initial, ep_type) %>%
  summarize(mean_error = mean(Error))

gc_plot <- gc_best_s %>%
  mutate(s_initial = ifelse(s_initial, 'sC', 'TR')) %>%
  ggplot() +
  geom_bar(aes(x=ep_type, y=mean_error, fill=ep_type), stat='identity') +
  facet_wrap(~ s_initial) + 
  xlab("Cluster type") +
  ylab("Mean error") + 
  scale_fill_discrete(guide="none") +
  ggtitle("Syllable\nContact\nComplex") +
  theme_classic(base_size=20) +
  theme(axis.text=element_text(size=16, angle = 45, vjust = 0.5, hjust = 1),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  ylim(-0.15, 0.15) +
  geom_hline(yintercept=0, size=0.75)
gc_plot
ggsave('figures/gc_best.png', height = 7, width = 7, units='in')


ggarrange(gs_plot, gc_plot, nrow=1)
ggsave('figures/best_errors.png', height=7, width=14, units='in')
