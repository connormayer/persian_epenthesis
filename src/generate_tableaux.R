library(tidyverse)

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

# Generate tableaux
experiment_df <- read_csv('data/experiment/experimental_results.csv') %>%
  separate(onset, into = c("first", "second"), sep = 1, remove = FALSE)

# Participant-specific epenthesis counts
counts_df <- experiment_df %>%
  group_by(word, ep_type, participant) %>%
  summarize(count=n()) %>%
  ungroup() %>%
  complete(word, ep_type, participant, fill=list(count = 0)) %>% 
  arrange(-participant) %>%
  inner_join(experiment_df %>% 
               select(word, sonority_delta, first, second) %>%
               distinct(),
             by=c('word'))


#reconfigure df 
counts_df <- counts_df %>%
  mutate(syllable_contact = ifelse(
    ep_type == 'prothesis',
    sonority[second] - sonority[first],
    ifelse(ep_type == 'anaptyxis',
           sonority[second] - 5, NA)))

# Global counts
global_counts <- counts_df %>%
  group_by(word, ep_type, sonority_delta, first, second, syllable_contact) %>%
  summarize(count=sum(count)) %>%
  arrange(word, fct_relevel(ep_type, "none"))

create_individual_tableaux <- function(data, constraints, out_folder, name_template) {
  dir.create(out_folder, showWarnings = FALSE)
  for (p in unique(data$participant)) {
    participant_counts <- data %>%
      filter(participant == p)
    filename <- file.path(out_folder, str_glue("{name_template}_p{p}.csv"))
    create_tableaux(participant_counts, constraints, filename)
  }
}


# Function that builds relevant tableaux given a list of constraints
create_tableaux <- function(data, constraints, output_file) {
  headers <- c('Input', 'Output', 'Frequency', constraints)
  
  write.table(
    matrix(headers, nrow=1),
    file=output_file,
    row.names=FALSE,
    col.names = FALSE,
    sep=','
  )
  
  for (i in 1:nrow(data)) {
    row <- data[i,]
    candidate <- c(row$word, str_c(row$word, row$ep_type, sep='-'), row$count)
    
    for (constraint in constraints) {
      if (constraint == '*Complex') {
        violation <- ifelse(row$ep_type == 'none', 1, '')
      } else if (constraint == 'Dep') {
        violation <- ifelse(row$ep_type == 'none', '', 1)
      } else if (constraint == 'Contiguity') {
        violation <- ifelse(row$ep_type == 'anaptyxis', 1, '')
      } else if (constraint == 'SyllableContact') {
        violation <- ifelse(
          row$ep_type == 'prothesis' & row$sonority_delta > 0, 1, ''
        )
      } else if (constraint == 'SyllableContact_4') {
        violation <- ifelse( ! is.na(row$syllable_contact) &
        row$syllable_contact >= 4, 1, ''
        )
      } else if (constraint == 'SyllableContact_3') {
        violation <- ifelse(! is.na(row$syllable_contact) &
          row$syllable_contact >= 3, 1, ''
        )
      } else if (constraint == 'SyllableContact_2') {
          violation <- ifelse(! is.na(row$syllable_contact) &
            row$syllable_contact >= 2, 1, ''
          )
      } else if (constraint == 'SyllableContact_1') {
        violation <- ifelse(! is.na(row$syllable_contact) &
           row$syllable_contact >= 1, 1, ''
        )
      } else if (constraint == 'SyllableContact_-1') {
        violation <- ifelse(! is.na(row$syllable_contact) &
           row$syllable_contact >= -1,
          1, ''
        )
      } else if (constraint == 'SyllableContact_-2') {
        violation <- ifelse(! is.na(row$syllable_contact) &
          row$syllable_contact >= -2, 
          1, ''
        )
      } else if (constraint == 'SyllableContact_-3') {
        violation <- ifelse(! is.na(row$syllable_contact) &
         row$syllable_contact >= -3, 
          1, ''
        )
      } else if (constraint == 'SyllableContact_-4') {
        violation <- ifelse(! is.na(row$syllable_contact) &
          row$syllable_contact >= -4, 
          1, ''
        )
      } else if (constraint == 'SyllableContact_-5') {
        violation <- ifelse(! is.na(row$syllable_contact) &
          row$syllable_contact >= -5, 
          1, ''
        )
     
      } else if (constraint == 'C/V') {
        violation <- ifelse(
          row$ep_type == 'anaptyxis', '', 1
        )
      } else if (constraint == 'L-Anchor') {
        violation <- ifelse(
          row$ep_type == 'prothesis', 1, ''
        )
      } else if (constraint == 'Dep-[ə]/S_T') {
        violation <- ifelse(
          row$first == 's' & row$ep_type == 'anaptyxis' & row$syllable_contact == -5, 
          1, ''
        )
      } else if (constraint == 'Dep-[ə]/S_N') {
        violation <- ifelse(
          row$first == 's' & row$ep_type == 'anaptyxis' & row$syllable_contact <= -3, 1, ''
        )
      } else if (constraint == 'Dep-[ə]/S_L') {
        violation <- ifelse(
          row$ep_type == 'anaptyxis' & 
            row$first == 's' & row$syllable_contact <= -2, 
          1, ''
        )
      } else if (constraint == 'Dep-[ə]/S_W') {
        violation <- ifelse(
          row$ep_type == 'anaptyxis' & row$first == 's' & 
            row$syllable_contact <= -1,
          1, ''
        )
      } else if (constraint == 'Dep-[ə]/O_R') {
        violation <- ifelse(
          row$ep_type == 'anaptyxis',
          1, ''
        )
      } else if (constraint == '*Complex-S') {
        violation <- ifelse(row$ep_type == 'none' & row$first == 's', 1, '')
      } else if (constraint == '*Complex-T') {
        violation <- ifelse(row$ep_type == 'none' & row$first != 's', 1, '')
      } else {
        print(str_c("Unknown constraint ", constraint))
        stop()
      }
      candidate <- c(candidate, violation)
    }
    
    write.table(
      matrix(candidate, nrow=1),
      file=output_file,
      row.names=FALSE,
      col.names = FALSE,
      append=TRUE,
      sep=','
    )
  }
}

# Create Gouskova simple tableau with global counts
create_tableaux(
  global_counts,
  c('SyllableContact', '*Complex', 'Dep', 'Contiguity'),
  'data/tableaux/gouskova_simple_global.csv')

# And with participant-specific counts
create_individual_tableaux(
  counts_df,
  c('SyllableContact', '*Complex', 'Dep', 'Contiguity'),
  'data/tableaux/gouskova_simple_ind', 'gs')

# Create Gouskova complex tableau with global counts
create_tableaux(
  global_counts, 
  c('SyllableContact_4', 'SyllableContact_3', 'SyllableContact_2', 
    'SyllableContact_1', 'SyllableContact_-1', 'SyllableContact_-2',
    'SyllableContact_-3', 'SyllableContact_-4', 'SyllableContact_-5',
    '*Complex', 'Dep', 'Contiguity'),
  'data/tableaux/gouskova_complex_global.csv')

# And with participant-specific counts
create_individual_tableaux(
  counts_df,
  c('SyllableContact_4', 'SyllableContact_3', 'SyllableContact_2', 
    'SyllableContact_1', 'SyllableContact_-1', 'SyllableContact_-2',
    'SyllableContact_-3', 'SyllableContact_-4', 'SyllableContact_-5',
    '*Complex', 'Dep', 'Contiguity'),
  'data/tableaux/gouskova_complex_ind', 'gc')

# Create Fleischhacker tableau with global counts
create_tableaux(
  global_counts, 
  c('C/V', 'L-Anchor', 'Contiguity', 'Dep-[ə]/S_T', 
    'Dep-[ə]/S_N', 'Dep-[ə]/S_L', 'Dep-[ə]/S_W', 'Dep-[ə]/O_R', '*Complex'),
  'data/tableaux/fleischhacker_global.csv')

# And with participant-specific counts
create_individual_tableaux(
  counts_df,
  c('C/V', 'L-Anchor', 'Contiguity', 'Dep-[ə]/S_T',
    'Dep-[ə]/S_N', 'Dep-[ə]/S_L', 'Dep-[ə]/S_W', 'Dep-[ə]/O_R', '*Complex'),
  'data/tableaux/fleischhacker_ind', 'fh')

# Create Gouskova simple split complex tableau with global counts
create_tableaux(
  global_counts,
  c('SyllableContact', '*Complex-S', '*Complex-T', 'Dep', 'Contiguity'),
  'data/tableaux/gouskova_simple_global_split.csv')

# And with participant-specific counts
create_individual_tableaux(
  counts_df,
  c('SyllableContact', '*Complex-S', '*Complex-T', 'Dep', 'Contiguity'),
  'data/tableaux/gouskova_simple_split_ind', 'gs_split')

# Create Gouskova complex split complex tableau with global counts
create_tableaux(
  global_counts, 
  c('SyllableContact_4', 'SyllableContact_3', 'SyllableContact_2', 
    'SyllableContact_1', 'SyllableContact_-1', 'SyllableContact_-2',
    'SyllableContact_-3', 'SyllableContact_-4', 'SyllableContact_-5',
    '*Complex-S', '*Complex-T', 'Dep', 'Contiguity'),
  'data/tableaux/gouskova_complex_global_split.csv')

# And with participant-specific counts
create_individual_tableaux(
  counts_df,
  c('SyllableContact_4', 'SyllableContact_3', 'SyllableContact_2', 
    'SyllableContact_1', 'SyllableContact_-1', 'SyllableContact_-2',
    'SyllableContact_-3', 'SyllableContact_-4', 'SyllableContact_-5',
    '*Complex-S', '*Complex-T', 'Dep', 'Contiguity'),
  'data/tableaux/gouskova_complex_split_ind', 'gc_split')

# Create Fleischhacker split complex tableau with global counts
create_tableaux(
  global_counts, 
  c('C/V', 'L-Anchor', 'Contiguity', 'Dep-[ə]/S_T', 
    'Dep-[ə]/S_N', 'Dep-[ə]/S_L', 'Dep-[ə]/S_W', 'Dep-[ə]/O_R',
    '*Complex-S', '*Complex-T'),
  'data/tableaux/fleischhacker_global_split.csv')

# And with participant-specific counts
create_individual_tableaux(
  counts_df,
  c('C/V', 'L-Anchor', 'Contiguity', 'Dep-[ə]/S_T', 
    'Dep-[ə]/S_N', 'Dep-[ə]/S_L', 'Dep-[ə]/S_W', 'Dep-[ə]/O_R',
    '*Complex-S', '*Complex-T'),
  'data/tableaux/fleischhacker_split_ind', 'fh_split')
