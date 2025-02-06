
library(ape)
library(phangorn)
library(pegas)
library(Biostrings)
library(strex) # for str_elem
library(rjags)


#' Function to check eligibility
#'
#' This function counts the number of sequences in the input fasta file (DNAbin object)
#' checks if they are not less than the specified minimum sequence count
#' also checks if the alignment width is not less than the specified minimum sequence length
#'
#' @param dna_bin. DNAbin object
#' @param min_sequence_count. Minimum sequence count allowed
#' @param min_sequence_width. Minimum sequence width allowed
#'
#' @return Returns TRUE if the criteria are passed, otherwise FALSE.
#' @examples
#' eligible <- check_eligibility(dna_bin, 2, 9)
#'
#'@export
check_eligibility <- function(dna_bin, min_sequence_count, min_seq_width){
  # Calculate the number of sequences
  numseq <- nrow(dna_bin)

  # Calculate the sequence alignment length
  seqlen <- ncol(dna_bin)

  if (numseq >= min_sequence_count & seqlen >= min_seq_width){
    print("Passed eligibility")
    return(TRUE)
  } else if (numseq >= min_sequence_count & !(seqlen >= min_seq_width)){
    print(paste("Sequence length is shorter than", min_seq_width))
    return(FALSE)
  } else if ( !(numseq >= min_sequence_count) & seqlen >= min_seq_width){
    print(paste("Sequence count is less than", min_sequence_count))
    return(FALSE)
  }
}
