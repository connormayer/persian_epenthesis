library(tidyverse)
library(brms)

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
         ep_type = fct_relevel(ep_type, "none"))

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
  theme_classic() + 
  theme(text=element_text(size=20))

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





########################
# STATISTICAL ANALYSIS #
########################

# Fit Bayesian multinomial logistic regression model
m_base <- brm(ep_type ~ delta + preceding_v + onset_age + onset + (1|speaker) + (1|word),
         data=df, family="categorical", prior=c(set_prior("normal(0,3)")),
         chains=4, cores=4)
summary(m_base)

#model with NAP sonority 
m_base_2 <- brm(ep_type ~ nap_sonority + preceding_v + onset_age + onset + (1|speaker) + (1|word),
              data=df, family="categorical", prior=c(set_prior("normal(0,3)")),
              chains=4, cores=4)
summary(m_base_2)

# Plot posteriors
plot(m_base)

# Test whether rate of improvement at anaptyxis is slower than at prothesis
hyp_test <- hypothesis(m_base, 'muprothesis_onset_age < muanaptyxis_onset_age')
hyp_test
plot(hyp_test)
