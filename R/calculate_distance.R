
library(ape)
library(phangorn)
library(pegas)
library(Biostrings)
library(strex) # for str_elem
library(rjags)


#' Function to calculate pairwise distances from DNA sequences
#'
#' This function takes a DNAbin object,
#' checks eligibility using the "check_eligibility" function,
#' (if eligible), calculates the distance using the specified evolutionary model.
#' Currently supported distance metrics include:
#' - mean pairwise distance from the "raw" model (rawMPD)
#' - mean pairwise distance from the "TN93" model (tn93MPD)
#' - nucleotide diversity from the "raw" model (rawPI)
#' - nucleotide diversity from the "TN93" model (tn93PI)
#' - weighted fraction of polymorphic sites (WFPS)
#' - WFPS at third codon positions only (WFPS_codons)
#'
#' @param dna_bin. DNAbin object
#' @param distance_metric. specified distance metric
#' @param min_sequence_count. Minimum sequence count, used by "check_eligibility function".
#' @param min_seq_width. Minimum sequence alignment width, used by "check_eligibility function".
#' @param error_threshold. Sequencing error threshold. used by "diversity_at_each_position" function
#' and must be specified for WFPS or WFPS_codons
#' @param variants. Expected characters in DNAbin object. used by "diversity_at_each_position" function
#' and must be specified for WFPS or WFPS_codons
#'
#' @return distance as a numeric value.
#' @examples
#' dist_rawMPD <- calculate_distance(dna_bin, "rawMPD", 2, 9)
#'
#' @export
calculate_distance <- function(dna_bin, distance_metric,
                               min_sequence_count, min_seq_width,
                               errorthreshold = 0,
                               variants = c("A","C","G","T","-")){
  # check eligibility
  eligible = check_eligibility(dna_bin, min_sequence_count, min_seq_width)

  if (eligibile){

    if (distance_metric == "rawMPD"){dist <- mean(dist.dna(dna_bin, model = "raw"))}

    if (distance_metric == "tn93MPD"){dist <- mean(dist.dna(dna_bin, model = "TN93"))}

    if (distance_metric == "rawPI"){dist <- nuc.div(dna_bin, model = "raw")}

    if (distance_metric == "tn93PI"){dist <- nuc.div(dna_bin, model = "TN93")}

    if (distance_metric == "WFPS"){
      div_list <- diversity_at_each_position(dna_bin, errorthreshold,variants)
      div <- unlist(div_list)
      dist <- mean(div)
    }

    if (distance_metric == "WFPS_codons"){
      div_list <- diversity_at_each_position(dna_bin, errorthreshold,variants)
      div <- unlist(div_list)
      # trim to first and last non-gappy codons
      indices = locate_overall_non_gap_codons(dna_bin)
      first = indices[1]
      last = indices[2]
      # diversity at third codon positions
      if (!all(is.na(div))){
        div3 <- div[seq(3,length(div),3)]
        div3_trimmed = div3[first:last]
        dist = mean(div3_trimmed)
      } else {print("Diversity at all codon positions is NA")}

    }

    print(dist)
    return(dist)
  }
}
