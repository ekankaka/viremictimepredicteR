
rm(list=ls())

library(coda)
library(rjags)
library(Biostrings)

my_scripts = list.files(pattern = "\\.R$")

lapply(my_scripts, source)

#mcmc_samples = readRDS("../data/samples_linearRegressionOnly.rds")

# read fasta
dnaset_gp41 = readDNAStringSet("../data/example_gp41_outgrowth.fasta")
dnaset_RT = readDNAStringSet("../data/example_RT_outgrowth.fasta")

# remove reference sequences
noRefs_gp41 = remove_reference_sequences(dnaset = dnaset_gp41)
noRefs_RT = remove_reference_sequences(dnaset = dnaset_RT)

# Trim leading and trailing gappy codons, where all sequences in alignment have gappy codons
trimmed_gp41 = trim_terminal_gappy_codons(dnaset = noRefs_gp41)
trimmed_RT = trim_terminal_gappy_codons(dnaset = noRefs_RT)

# Ensure trimmed alignment width is at least the required minimum (returns TRUE or FALSE)
pass1_gp41 = check_alignment_width(dnaset = trimmed_gp41, min_alignment_width = 9)
pass1_RT = check_alignment_width(dnaset = trimmed_RT, min_alignment_width = 9)

# Remove sequences with non-nucleotide characters in trimmed alignment, beyond a specified threshold
notgappy_gp41 = filter_sequences_by_non_nucleotides(dnaset = trimmed_gp41, non_nucleotide_threshold = 0.25)
notgappy_RT = filter_sequences_by_non_nucleotides(dnaset = trimmed_RT, non_nucleotide_threshold = 0.25)

# Ensure remaining sequence count is at least the required minimum (returns TRUE or FALSE)
pass2_gp41 = count_eligible_sequences(dnaset = notgappy_gp41, min_eligible_count = 2)
pass2_RT = count_eligible_sequences(dnaset = notgappy_RT, min_eligible_count = 2)

# calculate distance
dist_gp41 <- calculate_distance(dnaset = notgappy_gp41)
dist_RT <- calculate_distance(dnaset = notgappy_RT)
dist_gp41_and_RT_Mean <- (dist_gp41 + dist_RT) / 2

# predict viremic time
viremic_time_gp41 <- predict_viremic_time(distances = dist_gp41, sequence_type = "outgrowth", hiv_region = "gp41")
viremic_time_RT <- predict_viremic_time(distances = dist_RT, sequence_type = "outgrowth", hiv_region = "RT")
viremic_time_gp41_and_RT_Mean <- predict_viremic_time(distances = dist_gp41_and_RT_Mean, sequence_type = "outgrowth", hiv_region = "gp41_and_RT_Mean")

