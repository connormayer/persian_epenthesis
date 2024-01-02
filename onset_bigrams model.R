library(tidyverse)
library(brms)

experimental_results <- read_csv("data/experiment/experimental_results.csv")

experimental_results <- experimental_results %>%
  mutate(onset_bigram = case_when(onset == "tr" ~ "T R", 
                                                        onset == "sn" ~ "S N", 
                                                        onset == "dw" ~ "D W",
                                                        onset == "kj" ~ "K Y", 
                                                        onset == "sw" ~ "S W",
                                                        onset == "gl" ~ "G L", 
                                                        onset == "sm" ~ "S M", 
                                                        onset == "st" ~ "S T", 
                                                        onset == "sl" ~ "S L", 
                                                        onset == "kw" ~ "K W", 
                                                        onset == "sk" ~ "S K", 
                                                        onset == "fl" ~ "F L", 
                                                        onset == "dr" ~ "D R", 
                                                        onset == "gr" ~ "G R", 
                                                        onset == "bj" ~ "B Y", 
                                                        onset == "br" ~ "B R", 
                                                        onset == "sp" ~ "S P", 
                                                        onset == "θr" ~ "TH R", 
                                                        onset == "pl" ~ "P L", 
                                                        onset == "tw" ~ "T W", 
                                                        onset == "kl" ~ "K L",
                                                        onset == "kr" ~ "K R", 
                                                        onset == "fr" ~ "F R")) 
  # select(onset_bigram) %>%
  # distinct() %>%
  # write_csv("test_data.csv", col_names = FALSE)

  
output <- read_csv("data/phonotactic_probability/output.csv") %>%
  select(word, bi_prob_smoothed) %>%
  rename("onset_bigram" = word)

experimental_results <- left_join(experimental_results, output, by = "onset_bigram") %>%
  write_csv("experimental_results_2.csv")
#moved to phonotactic probability folder 


experimental_results <-  experimental_results %>%
  mutate(participant = as_factor(participant),
                                ep_type = fct_relevel(ep_type, "none"))


m_bigrams <- brm(
  ep_type ~ sonority_delta + onset + preceding_v + PC1 + context + bi_prob_smoothed + (1|participant) + (1|word),
  data = experimental_results,
  family="categorical",
  prior=c(set_prior("normal(0,3)")),
  chains=4, cores=4)

summary(m_bigrams)

hyp_test <- hypothesis(m_bigrams, 'muprothesis_PC1 < muanaptyxis_PC1')
hyp_test

plot(hyp_test)
