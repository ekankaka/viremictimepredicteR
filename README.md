
# viremictimepredicteR

<!-- badges: start -->
<!-- badges: end -->

The goal of viremictimepredicteR is to calculate diversity for a specified region of aligned HIV proviral sequences, and use this diversity to predict time since infection. Predictions are based on markov chain monte carlo samples from Bayesian models fitted by the authors, for these regions and distance metrics. Predictions include a mean estimate and a 95% credible interval, based on intercepts and slopes from the posterior distribution of a simple linear regression model.

The currently supported hiv_regions (selected due to better performance) include:
- Envgp120V5
- Polp51RT
- Gagp17Matrix
- Envgp41

Currently supported measures of proviral diversity include:
- mean pairwise distance from the "raw" model (rawMPD)
- mean pairwise distance from the "TN93" model (tn93MPD)
- nucleotide diversity from the "raw" model (rawPI)
- nucleotide diversity from the "TN93" model (tn93PI)
- weighted fraction of polymorphic sites (WFPS)
- WFPS at third codon positions only (WFPS_codons)

Note: For WFPS and WFPS_codons, variants and errorthreshold must be specified. For example:
- variants = c("A","C","G","T","-"), 
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
library(viremictimepredicteR)
## specify the input fasta file
input_fasta = "data/example_GAG_P17.fasta"

## load fasta file and remove reference sequences if any (sequences with pattern "Ref", or "HXB2" in the headers)
fasta_filtered = remove_reference_sequences(fasta_file_path = input_fasta)

# check eligibility
eligibile = check_eligibility(dna_bin = fasta_filtered, min_sequence_count = 2, min_seq_width = 9)

# calculate distance
dist_rawMPD = calculate_distance(dna_bin = fasta_filtered, distance_metric = "rawMPD", min_sequence_count = 2, min_seq_width = 9)

dist_WFPS = calculate_distance(dna_bin  =  fasta_filtered, distance_metric  =  "WFPS", min_sequence_count = 2, min_seq_width = 9, errorthreshold = 0, variants = c("A","C","G","T","-"))

# predict time since infection (and credible intervals)
TSI_rawMPD = predict_TSI(hiv_region = "Gagp17Matrix", distance_metric = "rawMPD", x_new = 0.003)
TSI_WFPS = predict_TSI(hiv_region = "Gagp17Matrix", distance_metric = "WFPS", x_new = 0.003)
```

