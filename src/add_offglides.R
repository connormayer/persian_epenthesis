library(tidyverse)

cmu <- read_csv('data/english_training_data.txt', col_names = "word")
training_data <- read_csv('data/onset_type_frequencies.txt', col_names = 'onset')

cmu <- cmu %>%
  mutate(onset = str_extract(word, '^\\w Y')) %>%
  filter(!is.na(onset)) %>%
  select(onset)

final_data <- training_data %>%
  rbind(cmu) %>%
  write_csv('data/training_data.csv', col_names = FALSE)
