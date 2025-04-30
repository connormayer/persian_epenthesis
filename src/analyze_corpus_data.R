library(tidyverse)
#library(brms)
library(scales)
library(lme4)

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

df <- df %>%
  separate(onset, into = c("first", "second"), sep = 1, remove = FALSE) %>%
  mutate(ep_type = plyr::revalue(ep_type, 
                           c("none" = "none", 
                             "anaptyxis" = 'medial\nepenthesis',
                             "prothesis" = 'pre-\nepenthesis')))

#################
# VISUALIZATION #
#################

# Simple bar plot of counts of different epenthesis tokens
df %>%
  ggplot(aes(x=ep_type, fill=ep_type)) +
  geom_bar() + 
  ylab("Number of tokens") +
  xlab("Epenthesis outcome") + 
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  scale_fill_discrete(guide="none")
ggsave('figures/corpus_epenthesis_counts.png', height=7, width=10, units='in')

# Broken down by onset
df %>%
  group_by(onset, s_initial) %>%
  summarize(ep_rate = mean(has_ep),
            ep_err = sd(has_ep)/sqrt(length(has_ep))) %>%
  ggplot() +
  geom_bar(aes(x=fct_reorder(onset, ep_rate), y=ep_rate, fill=s_initial), stat='identity') +
  geom_errorbar(aes(x=fct_reorder(onset, ep_rate), ymin=ep_rate-ep_err, ymax=ep_rate+ep_err), 
                width=0.5, alpha=0.9, size=0.5, position=position_dodge(width=0.9)) +
  xlab("Onset identity") +
  ylab("Mean epenthesis rate") +
  scale_fill_discrete(name = "sC onset?") +
  theme_classic(base_size=30) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"))
ggsave('figures/corpus_onset_ep_rates.png', height = 6, width = 12, units='in')


temp_df <-df %>%
  group_by(s_initial, ep_type) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(s_initial) %>%
  mutate(count = count / sum(count),
         s_initial = ifelse(s_initial, 'sC onsets', 'TR onsets'))

temp_df %>%
  ggplot(aes(x=ep_type, y=count, fill=ep_type)) +
  geom_bar(stat='identity') +
  facet_wrap(~ s_initial) +
  ylab("Proportion of tokens") +
  xlab("Epenthesis outcome") + 
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  scale_fill_discrete(guide="none") +
  ylim(0, 1)
ggsave('figures/corpus_epenthesis_by_onset_type.png', height = 7, width = 15, units='in')

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
  ggplot(aes(x=onset_age, y=ep_count)) +
  geom_point(size=4) +
  geom_smooth(method='lm') + 
  xlab("Age of English onset") + 
  ylab("Epenthesis count") +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"))
ggsave('figures/corpus_age_by_ep_count.png', height = 7, width = 10, units='in')

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
  filter(ep_type == 'pre-\nepenthesis')

# Plot prothesis rate against onset age
temp_df_2 %>%
  ggplot(aes(x=onset_age, y=freq)) +
  geom_point(size=4) +
  geom_smooth(method='lm') +
  xlab("Age of English onset") + 
  ylab(expression(paste("", frac("Pre-epenthesis count", "Total epenthesis count")))) +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  ylim(0, 1.1)
ggsave('figures/corpus_age_by_prothesis_prop.png', height = 7, width = 10, units='in')

########################
# STATISTICAL ANALYSIS #
########################

#simple logistic regression
simple_model <- glmer(
  has_ep ~ preceding_v + scale(onset_age) * s_initial + (1|speaker) + (1|onset), 
  data = df, family = 'binomial',
  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
)

summary(simple_model)
