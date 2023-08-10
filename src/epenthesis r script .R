require(nnet)
require(tidyverse)
require(brms)

#setwd("C:/Users/conno/Dropbox/ling/research/noah_project")
setwd("E:/Dropbox/ling/research/noah_project")
df <- read_csv("Persian_epenthesis.csv")

#define preceding vowwls 
vowels <- c("e", "i", "u", "o")

#create dictionary for sonority values 
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


#load leap q csv 
leap_q_df <- read_csv("leap_q.csv")

#join leap q csv to df 
epenthesis_df <- df %>%
  left_join(leap_q_df, by = "participant")

#Add two new columns: preceding_v indicates whether previous sound is a
# vowel. has_ep indicates whether epenthesis occurs or not
epenthesis_df <- epenthesis_df %>%
  mutate(preceding_v = ifelse(last_sound %in% vowels, TRUE, FALSE),
       has_ep = ifelse(ep_type %in% c("anaptyxis", "prothesis"), TRUE, FALSE))


#add sonority delta 
epenthesis_df <- epenthesis_df %>%
  separate(onset, into = c("first", "second"), sep = 1, remove = FALSE) %>%
  mutate(sonority_delta = sonority[second] - sonority[first])

#rename columns with spaces
epenthesis_df <- epenthesis_df %>%
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
  )

# Filter down to columns we're interested in and remove na values for ep_type
epenthesis_df <- epenthesis_df %>%
  filter(ep_type == "anaptyxis" | ep_type == "prothesis" |ep_type == "none") %>%
  select(participant,
         word,
         onset,
         ep_type,
         has_ep,
         preceding_v,
         context,
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
    )

# Convert ep_type and participant to factors and set the reference level to be
# 'none' for ep_type
epenthesis_df$participant <- as_factor(epenthesis_df$participant)
epenthesis_df$ep_type <- relevel(as.factor(epenthesis_df$ep_type), ref = "none")

# Get subset of columns we'll do PCA on
leap_q_subset <- epenthesis_df %>%
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
  )

# Scale columns (good practice before doing a PCA)
scaled_df <- leap_q_subset %>%
  mutate(across(where(is.numeric), scale))

# Do the PCA
pca_input <- leap_q_subset %>% select(-participant)
pca <- prcomp(pca_input)

# Convert LEAP-Q scores into corresponding PC scores
predicted <- data.frame(predict(pca, pca_input))
predicted$participant <- leap_q_subset$participant
predicted <- predicted %>% select(participant, PC1, PC2) %>% unique()

# Add these scores to original dataframe
epenthesis_df <- inner_join(epenthesis_df, predicted, by='participant')

# Simple bar plot of counts of different epenthesis times
epenthesis_df %>%
  filter(ep_type == "anaptyxis" | ep_type == "prothesis" |ep_type == "none") %>%
  ggplot(aes(x = ep_type)) +
  geom_bar()

#barplot of counts of anaptyxis vs prothesis 
epenthesis_df %>%
  filter(ep_type == "anaptyxis" | ep_type == "prothesis") %>%
  ggplot(aes(x = ep_type)) +
  geom_bar()

# Heatmap plotting counts of epenthesis type by sonority delta
epenthesis_df %>%
  group_by(sonority_delta, ep_type) %>%
  filter(ep_type == "anaptyxis" | ep_type == "prothesis" |ep_type == "none") %>%
  summarize(count = n()) %>%
  ggplot(aes(sonority_delta, ep_type, fill=count)) +
  geom_tile() +
  geom_text(aes(label=count)) +
  scale_fill_gradient(trans='log', breaks=c(0, 2, 8, 32, 128, 512),
                      low="white", high="darkblue")

# Heatmap plotting epenthesis type by proceding vowel
epenthesis_df %>%
  filter(ep_type == "anaptyxis" | ep_type == "prothesis" |ep_type == "none") %>%
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
# This transformation is a bit more complex, but the idea is that we're convert
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
  ggplot(aes(x=PC1, y=freq)) +
  geom_point() +
  geom_smooth(method='lm') +
  xlab("Relative Farsi Dominance") + 
  ylab("Prothesis count / Total epenthesis count") +
  ggtitle("Relative rate of prothesis increases with L2 proficiency") +
  theme_classic(base_size=15)

# I think these plots are actually easier to interpret without the regression
# lines

# Plot epenthesis rate against PC1
ep_rate_df %>%
  ggplot(aes(x=PC1, y=ep_rate)) +
  geom_point() +
  #geom_smooth(method='lm') +
  ylab("Proportion of complex onsets epenthesized")

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

# Fit a model with none as the reference level
m_base <- brm(
  ep_type ~ sonority_delta + onset + preceding_v + PC1 + context + (1|participant) + (1|word), 
  data = epenthesis_df, 
  family="categorical", 
  prior=c(set_prior("normal(0,3)")),
  chains=4, cores=4)
summary(m_base)
# This makes a series of plots of the probability distributions
# the model has calculated for each parameter
plot(m_base)

# I can explain how this works in more detail later, but in terms of interpreting
# the output, Estimate is the estimated difference between these two parameters,
# CI.Lower and CI.Upper are the 95% credible intervals, and post.prob is the
# probability under this model that muanaptyxis_PC1 > muprothesis_PC1. For the
# data we have now this comes out as 'weakly significant'. The CI is [-0.08, 0.02]
# and the posterior probability is just over 86%.
hyp_test <- hypothesis(m_base, 'muprothesis_PC1 < muanaptyxis_PC1')
hyp_test

# This plots the probability distribution over the difference between these two
# coefficients. You can see the peak at about -0.03
plot(hyp_test)

epenthesis_df$word <- as_factor(epenthesis_df$word)

epenthesis_df %>%
  group_by(word, ep_type) %>%
  count(.drop=FALSE)

