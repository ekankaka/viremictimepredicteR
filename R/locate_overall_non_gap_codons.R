library(ape)
library(phangorn)
library(pegas)
library(Biostrings)
library(strex) # for str_elem
library(rjags)

#' Function to locate the portion of an alignment between the first and last non-gap codons,
#' in any of the sequences in the alignment.
#'
#' @param dna_bin. DNAbin object
#'
#' @return index1,index2.
#' @examples
#' indices = locate_overall_non_gap_codons(dna_bin)
#' first = indices[1]
#' last = indices[2]
#'
#' @export
locate_overall_non_gap_codons <- function(dna_bin) {
  # Convert DNAbin object to character matrix
  dna_char <- as.character(dna_bin)

  # Get sequence length (assuming all sequences are the same length)
  seq_length <- ncol(dna_char)

  # Ensure the sequence length is divisible by 3 for codon structure
  if (seq_length %% 3 != 0) {
    stop("Sequence length is not a multiple of 3. Check your alignment.")
  }

  # Break sequences into codons (3 columns at a time)
  num_codons <- seq_length / 3
  codon_matrix <- matrix(NA, nrow = nrow(dna_char), ncol = num_codons)

  for (i in 1:num_codons) {
    codon_start <- (i - 1) * 3 + 1
    codon_end <- codon_start + 2
    codon_matrix[, i] <- apply(dna_char[, codon_start:codon_end], 1, paste, collapse = "")
  }

  # Identify codon positions with no gaps ('-' or 'N') across all sequences
  non_gap_codons <- apply(codon_matrix, 2, function(codon_column) {
    all(!grepl("[-N]", codon_column))  # TRUE if no gaps in all sequences at this codon
  })

  # Find first and last non-gap codon positions
  if (any(non_gap_codons)) {
    first_non_gap <- which(non_gap_codons)[1]
    last_non_gap <- tail(which(non_gap_codons), 1)
  } else {
    first_non_gap <- NA
    last_non_gap <- NA
  }

  return(c(first_non_gap_codon = first_non_gap, last_non_gap_codon = last_non_gap))
}
