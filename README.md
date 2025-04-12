
# viremictimepredicteR

<!-- badges: start -->
<!-- badges: end -->

The goal of viremictimepredicteR is to calculate diversity for a specified region of aligned HIV proviral sequences or outgrowth virus sequences, and use this diversity to predict time since infection. Predictions are based on markov chain monte carlo samples from Bayesian models fitted by the authors, for these regions and distance metrics. Predictions include a median estimate and a 95% credible interval, based on intercepts and slopes from the posterior distribution of a simple linear regression model, with or without weights.

The currently supported hiv regions (selected due to better performance) include:
- Matrix (Matrix region of gag, especially p17)
- RT (reverse transcriptase region of Pol)
- gp41 (gp41 region of env)

Currently supported measures of proviral diversity include:
- mean pairwise distance from the "raw" model (rawMPD)
- mean pairwise distance from the "TN93" model (tn93MPD)
- nucleotide diversity from the "raw" model (rawPI)
- nucleotide diversity from the "TN93" model (tn93PI)
- weighted fraction of polymorphic sites (WFPS).
- WFPS at third codon positions only (WFPScodons). Available but not recommended for viremic time prediction.

Note: *For WFPS and WFPScodons, sequencing errorthresholds are required to calculate these metrics. 
 For example:
- errorthreshold = 0 (for prominent outgrowth viruses )
- errorthreshold = 0.01 (for usual NGS e.g proviral DNA sequences that are not outgrowth)

## Installation

You can install the development version of viremictimepredicteR from [GitHub](https://github.com/) with: 

``` r
devtools::install_github("ekankaka/viremictimepredicteR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
## load package
library(viremictimepredicteR)

## specify the input fasta file
## preferably, sequences in this input fasta file should:
## - NOT contain dual (or multiple) infections
## - be filtered for APOBEC-induced G to A hypermutation
## - contain unique sequences only (no identical sequences)
## - contain a minimum of two eligible sequences as per the above criteria.
## - be codon-alined and cut to a specific hiv region (One great tool is Gene cutter from los alamos) 

# read fasta
dnaset = readDNAStringSet("data/example_gp41_outgrowth.fasta")

# remove reference sequences
noRefs = remove_reference_sequences(dnaset = dnaset)

# Trim leading and trailing gappy codons, where all sequences in alignment have gappy codons
trimmed = trim_terminal_gappy_codons(dnaset = noRefs)

# Ensure trimmed alignment width is at least the required minimum (returns TRUE or FALSE)
pass1 = check_alignment_width(dnaset = trimmed)

# Remove sequences with non-nucleotide characters in trimmed alignment, beyond a specified threshold
notgappy = filter_sequences_by_non_nucleotides(dnaset = trimmed, non_nucleotide_threshold = 0.25)

# Ensure remaining sequence count is at least the required minimum (returns TRUE or FALSE)
pass2 = count_eligible_sequences(dnaset = notgappy, min_eligible_count = 2)

# calculate distance
dist <- calculate_distance(dnaset = notgappy, sequence_type = "outgrowth")

# predict viremic time
viremic_time <- predict_viremic_time(distances = dist, hiv_region = "gp41")

```

