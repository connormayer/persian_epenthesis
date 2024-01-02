require(tidyverse)
require(brms)

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
df <- read_csv("data/experiment/experiment_responses.csv") %>%
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
         has_ep = ep_type %in% c("anaptyxis", "prothesis")) %>%
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

# Do the PCA
pca <- pca_input %>% 
  select(-participant) %>%
  prcomp()

# Convert LEAP-Q scores into corresponding PC scores
predicted <- data.frame(predict(pca, pca_input)) %>%
  mutate(participant=pca_input$participant) %>% 
  select(participant, PC1, PC2) %>% 
  unique()

# Add these scores to original dataframe
epenthesis_df <- inner_join(epenthesis_df, predicted, by='participant')

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
  pivot_longer(cols = c("sonority_delta", "nap_sonority"), 
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
  mutate(sonority_type = fct_relevel(
                          fct_recode(
                            sonority_type, 
                            `Traditional Sonority Δ` = "sonority_delta",
                            `NAP Sonority Δ` = "nap_sonority"), 
                          "Traditional Sonority Δ")) %>%
  ggplot(aes(sonority_val, count, color=ep_type)) +
  geom_line(lwd=3) +
  geom_point(size=5) + 
  xlab("Onset sonority Δ") + 
  ylab("Relative proportion") +
  ggtitle("Traditional sonority Δ predicts epenthesis type, \n but NAP sonority Δ does not") +
  theme_classic(base_size=30) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  labs(color='Epenthesis \ntype') +
  facet_grid(~ sonority_type, scales = 'free_x')

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
  group_by(participant, PC1, PC2) %>%
  summarize(ep_rate = mean(has_ep),
            ep_count = sum(has_ep),
            ana_rate = sum(ep_type == 'anaptyxis') / sum(!is.na(ep_type)),
            pro_rate = sum(ep_type == 'prothesis') / sum(!is.na(ep_type)),
            ana_count = sum(ep_type == 'anaptyxis'),
            pro_count = sum(ep_type == 'prothesis'))

# Look at correlation between the two
cor(ep_rate_df$PC1, ep_rate_df$ep_rate)
cor(ep_rate_df$PC2, ep_rate_df$ep_rate)

# Look at proportion of anaptyxis vs. prothesis as a function of onset age
# This transformation is a bit more complex, but the idea is that we're converting
# counts of anaptyxis and prothesis into proportions of each within each speaker.
# We then only keep prothesis rate (because anaptyxis rate is just 1 - prothesis rate)
temp_df_2 <- epenthesis_df %>%
  filter(has_ep) %>% 
  group_by(participant, PC1, PC2, ep_type) %>%
  # Count number of epenthesis types
  summarize(count = n()) %>%
  # Normalize counts by total count for each speaker
  mutate(freq = count / sum(count)) %>%
  # Keep only proportions for prothesis
  filter(ep_type == 'prothesis')

# Plot prothesis proportion against PC1
temp_df_2 %>%
  ggplot(aes(x=-PC1, y=freq)) +
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
  ggplot(aes(x=-PC1, y=ep_rate)) +
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

# Create experimental results df
write_csv(epenthesis_df, 'data/experiment/experimental_results.csv')

# Fit a model with nap_sonority
m_nap <- brm(
  ep_type ~ nap_sonority + onset + preceding_v + PC1 + context + (1|participant) + (1|word),
  data = epenthesis_df,
  family="categorical",
  prior=c(set_prior("normal(0,3)")),
  chains=4, cores=4)
summary(m_nap)
m_nap <- add_criterion(m_nap, criterion = c("loo", "waic"))


#fit a model with sonority_delta
m_sonority_d <- brm(
  ep_type ~ sonority_delta + onset + preceding_v + PC1 + context + (1|participant) + (1|word),
  data = epenthesis_df,
  family="categorical",
  prior=c(set_prior("normal(0,3)")),
  chains=4, cores=4)
summary(m_sonority_d)
m_sonority_d <- add_criterion(m_sonority_d, criterion = c("loo", "waic"))

# This makes a series of plots of the probability distributions
# the model has calculated for each parameter
plot(m_sonority_d)

#hypothesis testing
hyp_test <- hypothesis(m_nap, 'muprothesis_PC1 < muanaptyxis_PC1')
hyp_test

hyp_test_2 <- hypothesis(m_sonority_d, 'muprothesis_PC1 < muanaptyxis_PC1')
hyp_test_2

# This plots the probability distribution over the difference between these two
# coefficients. You can see the peak at about -0.03
plot(hyp_test_2)
