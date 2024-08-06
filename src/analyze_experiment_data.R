library(tidyverse)
library(lme4)
library(car)



# Define preceding vowwls 
vowels <- c("e", "i", "u", "o")

# Standard sonority scale
sonority <- c(t = 0,
              d = 0,
              p = 0,
              b = 0,
              k = 0,
              g = 0,
              s = 1,
              f = 1,
              θ = 1,
              n = 2,
              m = 2,
              l = 3,
              r = 3,
              j = 4,
              w = 4)

# NAP sonority 
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

# add sonority deltas
#df <- read_csv("data/experiment/experiment_responses.csv") %>%
df <- read_csv("data/experiment/revised_results.csv") %>%
  separate(onset, into = c("first", "second"), sep = 1, remove = FALSE) %>%
  # C2 - C1
  mutate(sonority_delta = sonority[second] - sonority[first]) %>%
  # (V - C1) + (C2 - C1)
  mutate(nap_sonority = (4 - nap_sonority[first]) + 
      (nap_sonority[second] - nap_sonority[first]))

# load leap q csv 
leap_q_df <- read_csv("data/experiment/leap_q.csv")

# join leap q csv to df 
epenthesis_df <- df %>%
  inner_join(leap_q_df, by = "participant") %>%
# Add two new columns: preceding_v indicates whether previous sound is a
# vowel. has_ep indicates whether epenthesis occurs or not
  mutate(preceding_v = last_sound %in% vowels,
         has_ep = ep_type %in% c("anaptyxis", "prothesis"),
         s_initial = str_starts(onset, 's')) %>%
  # Rename a few columns
  rename(
    dominant_language = `dominant language`,
    current_farsi_exposure = `current farsi exposure`,
    current_english_exposure = `current english exposure`,
    age_of_farsi_acquisition = `age of farsi acquisition`,
    age_of_farsi_fluency = `age of farsi fluency`,
    length_of_farsi_residence = `length of farsi residence`,
    farsi_speaking_proficiency = `self reported farsi speaking proficiency`,
    farsi_understanding_proficiency = `self reported farsi understadning proficiency`,
    farsi_reading_proficiency = `self reported farsi reading`,
    age_of_english_acquisition = `age of english acquisition`,
    age_of_english_fluency = `age of english fluency`, 
    length_of_english_residence = `length of english residence`, 
    english_speaking_proficiency = `self reported english speaking proficency`,
    english_understanding_proficiency = `self reported english understanding proficiency`,
    english_reading_proficiency = `self reported english reading proficiency`
  ) %>%
  # Keep only valid ep types, remove NAs
  filter(ep_type %in% c("anaptyxis", "prothesis", "none")) %>%
  # Keep only columns we need
  select(participant,
         word,
         onset,
         ep_type,
         has_ep,
         preceding_v,
         context,
         nap_sonority,
         sonority_delta,
         s_initial,
         gender,
         dominant_language,
         current_farsi_exposure,
         current_english_exposure,
         length_of_farsi_residence,
         farsi_speaking_proficiency,
         farsi_understanding_proficiency,
         farsi_reading_proficiency,
         age_of_english_acquisition,
         age_of_english_fluency,
         length_of_english_residence,
         english_speaking_proficiency,
         english_understanding_proficiency,
         english_reading_proficiency
  ) %>%
  mutate(participant = as_factor(participant),
         ep_type = fct_relevel(ep_type, "none")) %>%
  filter(word != 'spreading')

# Get subset of columns we'll do PCA on
pca_input <- epenthesis_df %>%
  select(participant,
         current_farsi_exposure,
         current_english_exposure,
         length_of_farsi_residence,
         farsi_speaking_proficiency,
         farsi_understanding_proficiency,
         farsi_reading_proficiency,
         age_of_english_acquisition,
         age_of_english_fluency,
         length_of_english_residence,
         english_speaking_proficiency,
         english_understanding_proficiency,
         english_reading_proficiency
  ) %>%
  # Scale columns
  mutate(across(where(is.numeric), scale))

#create grouped leap_q inputs
pca_input_acquisition_exposure <- pca_input %>%
  select(participant, 
         current_farsi_exposure,
         age_of_english_acquisition,
         age_of_english_fluency,
         current_english_exposure)

pca_input_immersion <- pca_input %>%
  select(participant, 
         length_of_farsi_residence,
         length_of_english_residence) 

pca_input_self_reported_proficiency <- pca_input %>%
  select(participant, 
         farsi_speaking_proficiency,
         farsi_understanding_proficiency,
         farsi_reading_proficiency,
         english_speaking_proficiency,
         english_understanding_proficiency,
         english_reading_proficiency)
  

# Do the PCA
pca <- pca_input %>% 
  select(-participant) %>%
  prcomp(scale. = TRUE, center = TRUE)


#run PCA on different leap_q inputs
pca_acquisiton_exposure <-pca_input_acquisition_exposure %>% 
  select(-participant) %>%
  prcomp() 

pca_acquisiton_exposure_2 <-pca_input_acquisition_exposure %>% 
  select(-participant) %>%
  prcomp(scale. = TRUE) 

pca_immersion <- pca_input_immersion %>%
  select(-participant) %>%
  prcomp()

pca_self_reported_proficiency <- pca_input_self_reported_proficiency %>%
  select(-participant) %>%
  prcomp()

# Convert ungrouped LEAP-Q scores into corresponding PC scores
predicted <- data.frame(predict(pca, pca_input)) %>%
  mutate(participant=pca_input$participant) %>% 
  select(participant, PC1, PC2) %>% 
  unique()


#convert grouped leap_q scores into corresponding PC scores 
predicted_acquisition_exposure <- data.frame(predict(pca_acquisiton_exposure, pca_input_acquisition_exposure)) %>%
  mutate(participant=pca_input$participant) %>% 
  select(participant, PC1) %>% 
  rename(PC1_acquisition_exposure = PC1) %>%
  unique()


predicted_immersion <- data.frame(predict(pca_immersion, pca_input_immersion)) %>%
  mutate(participant=pca_input$participant) %>% 
  select(participant, PC1) %>% 
  unique() %>%
  rename(PC1_immersion = PC1) 


predicted_self_reported_proficiency <- data.frame(predict(pca_self_reported_proficiency, 
                                                          pca_input_self_reported_proficiency)) %>%
  mutate(participant=pca_input$participant) %>% 
  select(participant, PC1) %>% 
  unique() %>%
  rename(PC1_self_report = PC1)

#create df with all the PC scores
PC_df <- inner_join(predicted_acquisition_exposure, predicted_immersion)
PC_df <- inner_join(PC_df, predicted_self_reported_proficiency)
PC_df <- inner_join(PC_df, predicted)

PC_df <- PC_df %>%
  rename(PC1_original = PC1, 
         PC2_original = PC2)


# Add all the PC scores to original dataframe
epenthesis_df <- inner_join(epenthesis_df, PC_df)


# Simple bar plot of counts of different epenthesis types
epenthesis_df %>%
  ggplot(aes(x = ep_type)) +
  geom_bar()

# Barplot of counts of anaptyxis vs prothesis 
epenthesis_df %>%
  filter(ep_type %in% c("anaptyxis", "prothesis")) %>%
  ggplot(aes(x = ep_type)) +
  geom_bar()

# Heatmap plotting counts of epenthesis type by sonority delta
epenthesis_df %>%
  group_by(sonority_delta, ep_type) %>%
  summarize(count = n()) %>%
  ggplot(aes(sonority_delta, ep_type, fill=count)) +
  geom_tile() +
  geom_text(aes(label=count)) +
  scale_fill_gradient(trans='log', breaks=c(0, 2, 8, 32, 128, 512),
                      low="white", high="darkblue")

epenthesis_df %>%
  pivot_longer(cols = c("sonority_delta", "nap_sonority"), 
               names_to = "sonority_type",
               values_to = 'sonority_val') %>%
  group_by(sonority_type, sonority_val, ep_type) %>%
  summarize(count = n()) %>%
  filter(ep_type != 'none') 

epenthesis_df %>%
  pivot_longer(cols = c("sonority_delta", "nap_sonority", "s_initial"), 
               names_to = "sonority_type",
               values_to = 'sonority_val') %>%
  group_by(sonority_type, sonority_val, ep_type) %>%
  summarize(count = n()) %>%
  filter(ep_type != 'none') %>%
  ungroup() %>%
  group_by(sonority_type, sonority_val) %>%
  mutate(count = count / sum(count)) %>%
  ungroup() %>%
  add_row(sonority_type = 'sonority_delta',
          sonority_val = 4,
          ep_type = 'prothesis',
          count=0) %>%
  add_row(sonority_type = 'sonority_delta',
          sonority_val = 1,
          ep_type = 'anaptyxis',
          count=0) %>%
  add_row(sonority_type = 'sonority_delta',
          sonority_val = -1,
          ep_type = 'anaptyxis',
          count=0) %>%
  add_row(sonority_type = 'nap_sonority',
          sonority_val = 4,
          ep_type = 'prothesis',
          count=0) %>%
  add_row(sonority_type = 's_initial',
          sonority_val = FALSE,
          ep_type = 'prothesis',
          count=0) %>%
  mutate(sonority_type = fct_relevel(
                          fct_recode(
                            sonority_type, 
                            `Traditional Sonority Δ` = "sonority_delta",
                            `NAP Sonority Δ` = "nap_sonority",
                            `sC onset?` = 's_initial'), 
                          "Traditional Sonority Δ")) %>%
  ggplot(aes(as.integer(sonority_val), count, color=ep_type)) +
  geom_line(lwd=3) +
  geom_point(size=5) + 
  xlab("Onset value") + 
  ylab("Relative proportion") +
  #ggtitle("Traditional sonority Δ predicts epenthesis type, \n but NAP sonority Δ does not") +
  theme_classic(base_size=30) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  labs(color='Epenthesis \ntype') +
  facet_grid(~ sonority_type, scales = 'free_x') +
  scale_x_continuous(breaks = breaks_width(1))
ggsave('figures/sonority_measures_experiment.png', height = 7, width = 18, units='in')

ep_rate_df_2 <- epenthesis_df %>%
  group_by(onset, sonority_delta, nap_sonority, s_initial) %>%
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


ep_rate_df_2 %>% filter(sonority_delta == 2) %>%
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
ggsave('figures/traditional_delta_2_experiment.png', height = 7, width = 10, units='in')

blah <- ep_rate_df_2 %>% filter(s_initial)
blah$onset <- fct_relevel(blah$onset, "sk", "sp", "st", "sn", "sm", "sl", "sw")

blah %>%
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
ggsave('figures/s_onsets_experiment.png', height = 7, width = 10, units='in')


# Heatmap plotting epenthesis type by proceding vowel
epenthesis_df %>%
  group_by(preceding_v, ep_type) %>%
  summarize(count = n()) %>%
  ggplot(aes(preceding_v, ep_type, fill=count)) +
  geom_tile() +
  geom_text(aes(label=count)) +
  scale_fill_gradient(trans='log', breaks=c(0, 2, 8, 32, 128, 512), low="white", high="darkblue")

# Plot showing relationship between LEAP-Q PCs and epenthesis rate
# The first step is to make a new df where ep_rate is the mean of has_ep
# This means it's the proportion of cases where epenthesis occurs
ep_rate_df <- epenthesis_df %>%
  group_by(participant, PC1_original, PC2_original) %>%
  summarize(ep_rate = mean(has_ep),
            ep_count = sum(has_ep),
            ana_rate = sum(ep_type == 'anaptyxis') / sum(!is.na(ep_type)),
            pro_rate = sum(ep_type == 'prothesis') / sum(!is.na(ep_type)),
            ana_count = sum(ep_type == 'anaptyxis'),
            pro_count = sum(ep_type == 'prothesis'))

# Look at correlation between the two
cor(ep_rate_df$PC1_original, ep_rate_df$ep_rate)
cor(ep_rate_df$PC2, ep_rate_df$ep_rate)

# Look at proportion of anaptyxis vs. prothesis as a function of onset age
# This transformation is a bit more complex, but the idea is that we're converting
# counts of anaptyxis and prothesis into proportions of each within each speaker.
# We then only keep prothesis rate (because anaptyxis rate is just 1 - prothesis rate)
temp_df_2 <- epenthesis_df %>%
  filter(has_ep) %>% 
  group_by(participant, PC1_original, PC2_original, ep_type) %>%
  # Count number of epenthesis types
  summarize(count = n()) %>%
  # Normalize counts by total count for each speaker
  mutate(freq = count / sum(count)) %>%
  # Keep only proportions for prothesis
  filter(ep_type == 'prothesis')

# Plot prothesis proportion against PC1
temp_df_2 %>%
  ggplot(aes(x=-PC1_original, y=freq)) +
  geom_point(size=8) +
  geom_smooth(method='lm') +
  xlab("Relative English Dominance") + 
  ylab(expression(frac("Prothesis count", "Total epenthesis count"))) +
  ggtitle("Relative prothesis rate \nincreases with L2 proficiency") +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"))
ggsave('figures/relative_prothesis.png', height = 7, width = 7, units='in')

# I think these plots are actually easier to interpret without the regression
# lines

# Plot epenthesis rate against PC1
ep_rate_df %>%
  ggplot(aes(x=-PC1_original, y=ep_rate)) +
  geom_point(size = 6) +
  geom_smooth(method='lm') +
  ylab("Proportion of epenthesis") +
  xlab("Relative English Dominance") + 
  ggtitle("Overall rate of epenthesis \n decreases with L2 proficiency") +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"))
ggsave('figures/overall_epenthesis.png', height = 7, width = 7, units='in')

# Plot ana rate against PC1
ep_rate_df %>%
  ggplot(aes(x=PC1, y=ana_rate)) +
  geom_point() +
  #geom_smooth(method='lm') +
  ylab("Proportion of total complex onsets with anaptyxis")

# Plot pro rate against PC1
ep_rate_df %>%
  ggplot(aes(x=PC1, y=pro_rate)) +
  geom_point() +
  #geom_smooth(method='lm') +
  ylab("Proportion of total complex onsets with prothesis")

pca_df <- epenthesis_df %>%
  group_by(participant) %>%
  mutate(RED = -mean(PC1)) %>%
  select(participant, RED) %>%
  unique()

pca_df %>%
  ggplot() +
  geom_histogram(aes(x=RED), binwidth=1) +
  ylab("Number of participants") +
  xlab("Relative English Dominance") + 
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"))
ggsave('figures/relative_english_dominance.png', height = 7, width = 10, units='in')

foo <- epenthesis_df %>%
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
ggsave('figures/exp_epenthesis_by_onset_type.png', height = 7, width = 10, units='in')

foo2 <- epenthesis_df %>%
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
ggsave('figures/exp_overall_epenthesis_by_onset_type.png', height = 7, width = 10, units='in')



# Create experimental results df
write_csv(epenthesis_df, 'data/experiment/experimental_revised_results.csv')

# # Fit a model with nap_sonority
# m_nap <- brm(
#   ep_type ~ nap_sonority + onset + preceding_v + PC1_original + context + (1|participant) + (1|word),
#   data = epenthesis_df,
#   family="categorical",
#   prior=c(set_prior("normal(0,3)")),
#   chains=4, cores=4,
#   save_pars = save_pars(all = TRUE))
# summary(m_nap)
# nap_loo <- loo(m_nap, k_threshold = 0.7)
# 
# 
# #fit a model with sonority_delta
# m_sonority_d <- brm(
#   ep_type ~ sonority_delta + onset + preceding_v + PC1_original + context + (1|participant) + (1|word),
#   data = epenthesis_df,
#   family="categorical",
#   prior=c(set_prior("normal(0,3)")),
#   chains=4, cores=4,
#   save_pars = save_pars(all = TRUE))
# summary(m_sonority_d)
# sonority_d_loo <- loo(m_sonority_d, k_threshold=0.7)


epenthesis_df <- read_csv('data/experiment/experimental_revised_results.csv')

simple_model <- glmer(
  has_ep ~ preceding_v + scale(PC1_original) * s_initial + context + 
    (1|participant) + (1|onset/word), 
  data = epenthesis_df, family = 'binomial',
  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

summary(simple_model)

model_no_im <- glmer(
  has_ep ~ preceding_v + s_initial + scale(PC1_acquisition_exposure) + 
    scale(PC1_self_report) + scale(PC1_acquisition_exposure):s_initial + 
    scale(PC1_self_report):s_initial + context + (1|participant) + (1|onset/word), 
  data = epenthesis_df, family = 'binomial',
  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

summary(model_no_im)

model_no_self <- glmer(
  has_ep ~ preceding_v + s_initial + scale(PC1_acquisition_exposure) + 
    scale(PC1_immersion) + scale(PC1_acquisition_exposure):s_initial + 
    scale(PC1_immersion):s_initial + context + (1|participant) + (1|onset/word), 
  data = epenthesis_df, family = 'binomial',
  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

summary(model_no_self)

model_no_acq <- glmer(
  has_ep ~ preceding_v + s_initial + scale(PC1_immersion) + 
    scale(PC1_self_report) + scale(PC1_immersion):s_initial + 
    scale(PC1_self_report):s_initial + context + (1|participant) + (1|onset/word), 
  data = epenthesis_df, family = 'binomial',
  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

summary(model_no_acq)

model_self <- glmer(
  has_ep ~ preceding_v + scale(PC1_self_report) * s_initial + context + 
    (1|participant) + (1|onset/word), 
  data = epenthesis_df, family = 'binomial',
  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

summary(model_self)

model_im <- glmer(
  has_ep ~ preceding_v + scale(PC1_immersion) * s_initial + context + 
    (1|participant) + (1|onset/word), 
  data = epenthesis_df, family = 'binomial',
  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

summary(model_im)

model_acq<- glmer(
  has_ep ~ preceding_v + scale(PC1_acquisition_exposure) * s_initial + context +
    (1|participant) + (1|onset/word), 
  data = epenthesis_df, family = 'binomial',
  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

summary(model_acq)

full_model <- glmer(
  has_ep ~ preceding_v + scale(PC1_acquisition_exposure) + scale(PC1_immersion) + 
    scale(PC1_self_report) + s_initial + context + scale(PC1_self_report):s_initial + 
    scale(PC1_immersion):s_initial + scale(PC1_acquisition_exposure):s_initial + 
    (1|participant) + (1|onset/word), 
  data = epenthesis_df, family = 'binomial',
  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

summary(full_model)

anova(model_acq, model_im, model_self, model_no_acq, model_no_self, model_no_im, simple_model, full_model)


# Model preferred by BIC
vif(model_acq)
# Model preferred by AIC, but predictors are highly correlated...
vif(model_no_im)

