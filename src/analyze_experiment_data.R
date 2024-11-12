library(tidyverse)
library(lme4)
library(car)

################
# PROCESS DATA #
################

# Define preceding vowwls 
vowels <- c("e", "i", "u", "o")

df <- read_csv("data/experiment/experimental_raw_results.csv") %>%
  separate(onset, into = c("first", "second"), sep = 1, remove = FALSE)

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
         ep_type = fct_recode(ep_type, 
                              'none' = 'none', 
                              'medial\nepenthesis' = 'anaptyxis',
                              'pre-\nepenthesis' = 'prothesis'),
         ep_type = fct_relevel(ep_type, "none")) %>%
  # Remove one trisegmental onset
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
  prcomp()

#run PCA on different leap_q inputs
pca_acquisiton_exposure <-pca_input_acquisition_exposure %>% 
  select(-participant) %>%
  prcomp() 

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
# Create experimental results df
write_csv(epenthesis_df, 'data/experiment/experimental_revised_results.csv')

#########
# PLOTS #
#########

# Plot general epenthesis rates
epenthesis_df %>%
  ggplot(aes(x=ep_type, fill=ep_type)) +
  geom_bar() + 
  ylab("Number of tokens") +
  xlab("Epenthesis outcome") + 
  # ggtitle("Overall rate of epenthesis \n decreases with L2 proficiency") +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  scale_fill_discrete(guide="none")
ggsave('figures/experiment_epenthesis_counts.png', height=7, width=10, units='in')

# Plot epenthesis rates by onset type
epenthesis_df %>%
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
ggsave('figures/experiment_ep_rates_by_onset.png', height = 6, width = 12, units='in')

# Plot showing relationship between LEAP-Q PCs and epenthesis rate
# The first step is to make a new df where ep_rate is the mean of has_ep
# This means it's the proportion of cases where epenthesis occurs
ep_rate_df <- epenthesis_df %>%
  group_by(participant, PC1_acquisition_exposure) %>%
  summarize(ep_rate = mean(has_ep),
            ep_count = sum(has_ep),
            ana_rate = sum(ep_type == 'anaptyxis') / sum(!is.na(ep_type)),
            pro_rate = sum(ep_type == 'prothesis') / sum(!is.na(ep_type)),
            ana_count = sum(ep_type == 'anaptyxis'),
            pro_count = sum(ep_type == 'prothesis'))

# Look at proportion of anaptyxis vs. prothesis as a function of onset age
# This transformation is a bit more complex, but the idea is that we're converting
# counts of anaptyxis and prothesis into proportions of each within each speaker.
# We then only keep prothesis rate (because anaptyxis rate is just 1 - prothesis rate)
temp_df <- epenthesis_df %>%
  filter(has_ep) %>% 
  group_by(participant, PC1_acquisition_exposure, ep_type) %>%
  # Count number of epenthesis types
  summarize(count = n()) %>%
  # Normalize counts by total count for each speaker
  mutate(freq = count / sum(count)) %>%
  # Keep only proportions for prothesis
  filter(ep_type == 'pre-\nepenthesis')


# Plot epenthesis rate against RED
ep_rate_df %>%
  ggplot(aes(x=-PC1_acquisition_exposure, y=ep_rate)) +
  geom_point(size = 6) +
  geom_smooth(method='lm') +
  ylab("Epenthesis rate") +
  xlab("Relative English Dominance (acquisition)") + 
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"))
ggsave('figures/experiment_overall_epenthesis_RED.png', height = 7, width = 7, units='in')

# Plot pre-epenthesis rate by RED
temp_df %>%
  ggplot(aes(x=-PC1_acquisition_exposure, y=freq)) +
  geom_point(size=4) +
  geom_smooth(method='lm') +
  xlab("Relative English Dominance (acquisition)") + 
  ylab("Pre-epenthesis rate") +
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  ylim(0, 1.1)
ggsave('figures/experiment_prothesis_by_red.png', height = 7, width = 10, units='in')

temp_df_2 <- epenthesis_df %>%
  group_by(s_initial, ep_type) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(s_initial) %>%
  mutate(count = count / sum(count),
         s_initial = ifelse(s_initial, 'sC onsets', 'TR onsets'))

temp_df_2 %>%
  ggplot(aes(x=ep_type, y=count, fill=ep_type)) +
  geom_bar(stat='identity') +
  facet_wrap(~ s_initial) +
  ylab("Proportion of tokens") +
  xlab("Epenthesis outcome") + 
  theme_classic(base_size=22) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold")) +
  scale_fill_discrete(guide="none") +
  ylim(0,1)
ggsave('figures/experiment_epenthesis_by_onset_type.png', height = 7, width = 10, units='in')

#########
# STATS #
#########
epenthesis_df <- read_csv('data/experiment/experimental_revised_results.csv')

# Compare different 
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

model_acq <- glmer(
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

# Best model
summary(model_acq)