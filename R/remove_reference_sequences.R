
library(ape)
library(phangorn)
library(pegas)
library(Biostrings)
library(strex) # for str_elem
library(rjags)

#' Function to remove reference sequences, if any
#'
#' This function checks for sequences with the pattern "Ref" or "HXB2"
#' in the sequence header, and removes these sequences
#'
#' @param fasta_file. Character. Fasta file path
#'
#' @return Returns filtered sequences as DNAbin object
#' @examples
#' dna_bin_filtered <- remove_reference_sequences("data/example_GAG_P17.fasta")
#'
#' @export
remove_reference_sequences <- function(fasta_file_path) {
  # Read in the fasta file
  sequences <- read.dna(fasta_file_path, format = "fasta")

  # Filter out sequences with "Ref" or "HXB2" in the headers
  headers <- rownames(sequences)
  filtered_indices <- !grepl("Ref|HXB2", headers, ignore.case = TRUE)
  filtered_sequences <- sequences[filtered_indices, ]

  return(filtered_sequences)
}
