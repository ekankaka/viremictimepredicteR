
# viremictimepredicteR

<!-- badges: start -->
<!-- badges: end -->

The goal of viremictimepredicteR is to predict viremic time using unique sequence diversity for a specified region of aligned HIV proviral or outgrowth sequences. Predictions are based on markov chain monte carlo samples from Bayesian models fitted by the authors, for these regions and diversity metrics.

HIV regions: The currently recommended regions to use include:
- `gp41` (predict using diversity for the gp41 region of env)
- `RT` (predict using diversity for the reverse transcriptase region of Pol)
- `gp41_and_RT_Mean` (predict using mean diversity for the gp41 and RT regions)
- `Matrix` (Matrix region of gag, especially p17)

Diversity metrics: The currently recommended metrics include:
- mean pairwise distance from the "raw" model (`rawMPD`)
- mean pairwise distance from the "TN93" model (`tn93MPD`)
- nucleotide diversity from the "raw" model (`rawPI`)
- nucleotide diversity from the "TN93" model (`tn93PI`)

Weights: The currently recommended options include:
- `None`: Model with no weights (all data points carry the same weight = 1)

## Installation

You can install the development version of viremictimepredicteR from [GitHub](https://github.com/) with: 

``` r
devtools::install_github("ekankaka/viremictimepredicteR")
```

## Example

This is a basic example of how to predict viremic time using diversity in gp41 and RT combined:

``` r
## load package
library(viremictimepredicteR)
library(Biostrings) # for readDNAStringSet
library(pegas) # for nuc.div and dist.dna in calculate_distances function
ibrary(rjags) # for as.mcmc to load Bayesian posterior samples

## input fasta file. For best results, ensure:
## - No dual (or multiple) infections
## - No hypermutation
## - 2 or more unique sequences
## - Codon-aligned sequences, cut to a specific hiv region 
## (e.g., using Gene cutter from Los Alamos) 
gp41_path <- system.file("extdata", "example_gp41_outgrowth.fasta", package = "viremictimepredicteR")
dnaset_gp41 <- readDNAStringSet(gp41_path)

RT_path <- system.file("extdata", "example_RT_outgrowth.fasta", package = "viremictimepredicteR")
dnaset_RT <- readDNAStringSet(RT_path)

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
notgappy_gp41 = filter_non_nucleotides(dnaset = trimmed_gp41, threshold = 0.25)
notgappy_RT = filter_non_nucleotides(dnaset = trimmed_RT, threshold = 0.25)

# Ensure remaining sequence count is at least the required minimum (returns TRUE or FALSE)
pass2_gp41 = count_eligible_sequences(dnaset = notgappy_gp41, min_eligible_count = 2)
pass2_RT = count_eligible_sequences(dnaset = notgappy_RT, min_eligible_count = 2)

# calculate distance
dist_gp41 <- calculate_distance(dnaset = notgappy_gp41)
dist_RT <- calculate_distance(dnaset = notgappy_RT)
dist_Mean <- (dist_gp41 + dist_RT) / 2

# predict viremic time using Bayesian simple linear regression, with no weights
viremic_time <- predict_viremic_time(distances = dist_Mean, sequence_type = "outgrowth",  
hiv_region = "gp41_and_RT_Mean")

# Results (predictions for the four recommended diversity metrics, with credible intervals)
View(viremic_time)

```

