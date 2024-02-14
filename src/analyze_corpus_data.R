library(tidyverse)
library(brms)
library(scales)

vowels <- c("[u]", "[ə]", "[ɔɪ]")

#########################
# LOAD AND PROCESS DATA #
#########################

df <- read_csv("data/corpus/corpus_data.csv") |>
  # Rename the columns with spaces since they're annoying to reference
  rename(
    delta = `sonority delta`,
    last_sound = `last sound`,
    onset_age = `age of English onset`,
    res_len = `length of English residence`,
    add_lang = `additonal language`
  ) |>
  # Add two new columns: preceding_v indicates whether previous sound is a
  # vowel. has_ep indicates whether epenthesis occurs or not.
  # Convert ep_type from a string to a factor and set the reference level to be
  # 'none'
  mutate(preceding_v = last_sound %in% vowels,
         has_ep = ep_type %in% c("anaptyxis", "prothesis"),
         ep_type = fct_relevel(ep_type, "none"),
         s_initial = str_starts(onset, 's'),
         onset = str_replace(onset, 'sc', 'sk'))

##add NAP sonority delta 
nap_sonority <- c(t = 1, 
                  d = 2,
                  p = 1,
                  b = 2,
                  k = 1,
                  g = 2,
                  s = 1,
                  f = 1,
                  θ = 1,
                  n = 3,
                  m = 3,
                  l = 3,
                  r = 3,
                  j = 4,
                  w = 4, 
                  e = 4,
                  i = 4,
                  o = 4,
                  u = 4, 
                  a = 4)

df <- df %>%
  separate(onset, into = c("first", "second"), sep = 1, remove = FALSE) %>%
  mutate(nap_sonority = (4 - nap_sonority[first]) + 
           (nap_sonority[second] - nap_sonority[first]))
view(df)


#################
# VISUALIZATION #
#################
# Let's remove some of these once we know what we want to include
# in the paper.

# Simple bar plot of counts of different epenthesis times
df %>%
  ggplot(aes(x=ep_type)) +
  geom_bar()

# Heatmap plotting counts of epenthesis type by sonority delta
df %>%
  group_by(delta, ep_type) %>%
  summarize(count = n()) %>%
  ggplot(aes(delta, ep_type, fill=count)) +
  geom_tile() +
  geom_text(aes(label=count)) +
  scale_fill_gradient(trans='log', breaks=c(0, 2, 8, 32, 128, 512),
                      low="white", high="darkblue")

foo <-df %>%
  group_by(s_initial, ep_type) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(s_initial) %>%
  mutate(count = count / sum(count),
         s_initial = ifelse(s_initial, 'sC onsets', 'OR onsets'))

foo %>%
  ggplot(aes(x=ep_type, y=count, fill=ep_type)) +
  geom_bar(stat='identity') +
  facet_wrap(~ s_initial) +
  ylab("Proportion of responses") +
  xlab("Epenthesis type") + 
  # ggtitle("Overall rate of epenthesis \n decreases with L2 proficiency") +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  scale_fill_discrete(guide="none")
ggsave('figures/corpus_epenthesis_by_onset_type.png', height = 7, width = 10, units='in')

foo2 <-df %>%
  group_by(s_initial, has_ep) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(s_initial) %>%
  mutate(count = count / sum(count),
         s_initial = ifelse(s_initial, 'sC onsets', 'OR onsets'),
         has_ep = ifelse(has_ep, 'epenthesis', 'no epenthesis'))

foo2 %>%
  ggplot(aes(x=has_ep, y=count, fill=has_ep)) +
  geom_bar(stat='identity') +
  facet_wrap(~ s_initial) +
  ylab("Proportion of responses") +
  xlab("Epenthesis type") + 
  # ggtitle("Overall rate of epenthesis \n decreases with L2 proficiency") +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  scale_fill_discrete(guide="none")
ggsave('figures/corpus_overall_epenthesis_by_onset_type.png', height = 7, width = 10, units='in')

# Line plot plotting counts of epenthesis type by sonority delta
df %>%
  count(delta, ep_type, .drop=FALSE) %>%
  ggplot(aes(delta, n, color=ep_type)) +
  geom_line(linewidth=1.5)  +
  geom_point(size=4) +
  xlab("Sonority delta") + 
  ylab("Count") +
  labs(color = "Ep. type") + 
  theme_classic() + 
  theme(text=element_text(size=20))

# Heatmap plotting epenthesis type by proceding vowel
df %>%
  group_by(preceding_v, ep_type) %>%
  summarize(count = n()) %>%
  ggplot(aes(preceding_v, ep_type, fill=count)) +
  geom_tile() +
  geom_text(aes(label=count)) +
  scale_fill_gradient(trans='log', breaks=c(0, 2, 8, 32, 128, 512), 
                      low="white", high="darkblue")

df %>%
  group_by(preceding_v, ep_type) %>%
  summarize(count = n()) %>%
  filter(ep_type != 'none') %>%
  ggplot(aes(x=ep_type, y=count, fill=preceding_v)) +
  geom_bar(stat='identity', position=position_dodge()) +
  xlab("Epenthesis type") + 
  ylab("Count") +
  labs(fill = "Preceding V?") + 
  theme_classic() + 
  theme(text=element_text(size=20))

# Plot showing relationship between onset age and epenthesis rate
# The first step is to make a new df where ep_rate is the mean of has_ep
# This means it's the proportion of cases where epenthesis occurs
ep_rate_df <- df %>%
  group_by(speaker, onset_age, res_len) %>%
  summarize(ep_rate = mean(has_ep),
            ep_count = sum(has_ep),
            ana_rate = sum(ep_type == 'anaptyxis') / sum(!is.na(ep_type)),
            pro_rate = sum(ep_type == 'prothesis') / sum(!is.na(ep_type)),
            ana_count = sum(ep_type == 'anaptyxis'),
            pro_count = sum(ep_type == 'prothesis'))

ep_rate_df %>%
  ggplot(aes(x=onset_age, y=ep_rate)) +
  geom_point() +
  geom_smooth(method='lm')

ep_rate_df %>%
  ggplot(aes(x=onset_age, y=ep_count)) +
  geom_point(size=3) +
  geom_smooth(method='lm') + 
  xlab("Age of English onset") + 
  ylab("Epenthesis count") +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"))
ggsave('figures/corpus_age_by_ep_count.png', height = 7, width = 10, units='in')


ep_rate_df %>%
  ggplot(aes(x=onset_age, y=ana_rate)) +
  geom_point() +
  geom_smooth(method='lm')

ep_rate_df %>%
  ggplot(aes(x=onset_age, y=ana_count)) +
  geom_point(size=3) +
  geom_smooth(method='lm') +
  xlab("Age of English onset") + 
  ylab("Anaptyxis count") +
  theme_classic() + 
  theme(text=element_text(size=20))

ep_rate_df %>%
  ggplot(aes(x=onset_age, y=pro_rate)) +
  geom_point() +
  geom_smooth(method='lm')

ep_rate_df %>%
  ggplot(aes(x=onset_age, y=pro_count)) +
  geom_point(size=3) +
  geom_smooth(method='lm') +
  xlab("Age of English onset") + 
  ylab("Prothesis count") +
  theme_classic() + 
  theme(text=element_text(size=20))

# Look at correlation between the two
cor(ep_rate_df$onset_age, ep_rate_df$ep_rate)

# Compare residency length against ep_rate
ep_rate_df %>%
  ggplot(aes(x=res_len, y=ep_rate)) +
  geom_point() +
  geom_smooth(method='lm')

# Look at proportion of anaptyxis vs. prothesis as a function of onset age
# This transformation is a bit more complex, but the idea is that we're converting
# counts of anaptyxis and prothesis into proportions of each within each speaker.
# We then only keep prothesis rate (because anaptyxis rate is just 1 - prothesis rate)
temp_df_2 <- df %>%
  filter(has_ep) %>%
  group_by(speaker, onset_age, ep_type) %>%
  # Count number of epenthesis types
  summarize(count = n()) %>%
  # Normalize counts by total count for each speaker
  mutate(freq = count / sum(count)) %>%
  # Keep only proportions for prothesis
  filter(ep_type == 'prothesis')

# Plot prothesis rate against onset age
temp_df_2 %>%
  ggplot(aes(x=onset_age, y=freq)) +
  geom_point(size=3) +
  geom_smooth(method='lm') +
  xlab("Age of English onset") + 
  ylab("Prothesis count / Total epenthesis count") +
  theme_classic() + 
  theme(text=element_text(size=25))
ggsave('figures/age_by_prothesis_prop_corpus.png', height = 7, width = 10, units='in')


bar <- df %>%
  pivot_longer(cols = c("delta", "nap_sonority", "s_initial"), 
               names_to = "sonority_type",
               values_to = 'sonority_val') %>%
  group_by(sonority_type, sonority_val, ep_type) %>%
  summarize(count = n()) %>%
  filter(ep_type != 'none') %>%
  ungroup() %>%
  group_by(sonority_type, sonority_val) %>%
  mutate(count = count / sum(count)) %>%
  ungroup() %>%
  add_row(sonority_type = 'delta',
          sonority_val = 3,
          ep_type = 'prothesis',
          count=0) %>%
  add_row(sonority_type = 'delta',
          sonority_val = -1,
          ep_type = 'anaptyxis',
          count=0) %>%
  add_row(sonority_type = 'delta',
          sonority_val = 1,
          ep_type = 'anaptyxis',
          count=0) %>%
  add_row(sonority_type = 's_initial',
          sonority_val = 0,
          ep_type = 'prothesis',
          count=0) %>%
  add_row(sonority_type = 's_initial',
          sonority_val = 1,
          ep_type = 'anaptyxis',
          count=0) %>%
  # add_row(sonority_type = 'nap_sonority',
  #         sonority_val = 4,
  #         ep_type = 'prothesis',
  #         count=0) %>%
  # add_row(sonority_type = 's_initial',
  #         sonority_val = FALSE,
  #         ep_type = 'prothesis',
  #         count=0) %>%
  mutate(sonority_type = fct_relevel(
    fct_recode(
      sonority_type, 
      `Traditional Sonority Δ` = "delta",
      `NAP Sonority Δ` = "nap_sonority",
      `Is sC cluster?` = 's_initial'), 
    "Traditional Sonority Δ"))

bar %>%
  ggplot(aes(as.integer(sonority_val), count, color=ep_type)) +
  geom_line(lwd=3) +
  geom_point(size=5) + 
  xlab("Onset value") + 
  ylab("Relative proportion") +
  #ggtitle("Traditional sonority Δ predicts epenthesis type, \n but NAP sonority Δ does not") +
  theme_classic(base_size=30) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  labs(color='Epenthesis \ntype') +
  facet_grid(~ sonority_type, scales = 'free_x') +
  scale_x_continuous(breaks = breaks_width(1))
ggsave('figures/sonority_measures_corpus.png', height = 7, width = 18, units='in')

ep_rate_df_2 <- df %>%
  group_by(onset, delta, nap_sonority, s_initial) %>%
  summarize(non_ep_rate = 1 - mean(has_ep),
            ana_rate = sum(ep_type == 'anaptyxis') / sum(!is.na(ep_type)),
            pro_rate = sum(ep_type == 'prothesis') / sum(!is.na(ep_type))) %>%
  pivot_longer(cols = c("non_ep_rate", "ana_rate", "pro_rate"), 
               names_to = "ep_group",
               values_to = 'rate') %>%
  mutate(ep_group = fct_relevel(
    fct_recode(
      ep_group, 
      `No epenthesis` = "non_ep_rate",
      `Anaptyxis` = "ana_rate",
      `Prothesis` = 'pro_rate'), 
    "No epenthesis"))


ep_rate_df_2 %>% filter(delta == 2) %>%
  ggplot(aes(x=onset, y=rate, fill=ep_group)) +
  geom_bar(stat='identity') +
  xlab("Onset") + 
  ylab("Relative proportion") +
  #ggtitle("Traditional sonority Δ predicts epenthesis type, \n but NAP sonority Δ does not") +
  theme_classic(base_size=30) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  labs(color='Epenthesis \ntype')
ggsave('figures/traditional_delta_2_corpus.png', height = 7, width = 10, units='in')

########################
# STATISTICAL ANALYSIS #
########################

# Fit Bayesian multinomial logistic regression model
m_base <- brm(ep_type ~ delta + preceding_v + onset_age + onset + (1|speaker) + (1|word),
         data=df, family="categorical", prior=c(set_prior("normal(0,3)")),
         chains=4, cores=4,
         save_pars = save_pars(all = TRUE))
summary(m_base)
m_base_loo <- loo(m_base, k_threshold=0.7)

#model with NAP sonority 
m_base_2 <- brm(ep_type ~ nap_sonority + preceding_v + onset_age + onset + (1|speaker) + (1|word),
              data=df, family="categorical", prior=c(set_prior("normal(0,3)")),
              chains=4, cores=4,
              save_pars = save_pars(all = TRUE))
summary(m_base_2)
m_base_2_loo <- loo(m_base_2, k_threshold=0.7)

#model with binary sonority
m_base_3 <- brm(ep_type ~ s_initial + preceding_v + onset_age + onset + (1|speaker) + (1|word),
                data=df, family="categorical", prior=c(set_prior("normal(0,3)")),
                chains=4, cores=4,
                save_pars = save_pars(all = TRUE))
summary(m_base_3)
m_base_3_loo <- loo(m_base_3, k_threshold=0.7)

loo_list_corp <- list(m_base_loo, m_base_2_loo, m_base_3_loo)
loo_ws_corp <- loo_model_weights(loo_list_corp)
loo_model_weights(loo_list_corp, method='pseudobma')
loo_model_weights(loo_list_corp, method='pseudobma', BB=FALSE)

# Plot posteriors
plot(m_base)

# Test whether rate of improvement at anaptyxis is slower than at prothesis
hyp_test <- hypothesis(m_base, 'muprothesis_onset_age < muanaptyxis_onset_age')
hyp_test
plot(hyp_test)

hyp_test_2 <- hypothesis(m_base_2, 'muprothesis_onset_age < muanaptyxis_onset_age')
hyp_test_2
plot(hyp_test_2)

hyp_test_3 <- hypothesis(m_base_3, 'muprothesis_onset_age < muanaptyxis_onset_age')
hyp_test_3
plot(hyp_test_3)
