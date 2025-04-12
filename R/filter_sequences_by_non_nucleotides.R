
#' This Function removes sequences containing non-nucleotide characters 
#' beyond a set threshold
#'
#' @param dnaset. DNAStringSet object
#'
#'@param min_alignment_width. Minimum required width for dnaset alignment
#'
#' @param valid_nucleotides. Acceptable nucleotide characters (upper or lower case)
#'
#' @param non_nucleotide_threshold. Max acceptable proportion of non-nucleotide characters in alignment
#'
#' @return Returns DNAStringSet object, filtered
#'
#' @examples
#' notgappy = filter_sequences_by_non_nucleotides(dnaset = trimmed, threshold = 0.25)
#'
#' @export
filter_sequences_by_non_nucleotides <- function(dnaset,
                                                valid_nucleotides = c("A", "C", "G", "T"),
                                                min_alignment_width = 9,
                                                non_nucleotide_threshold = 0.25) {
  # assert that DNAstringset exists
  if (!inherits(dnaset, "DNAStringSet")) stop("Input must be a DNAStringSet")
  
  # assert that alignment width is at least the required minimum
  pass1 = check_alignment_width(dnaset = dnaset, min_alignment_width)
  if (pass1 == FALSE) stop("Alignment width less than the required minimum")
  
  # Function to calculate percent of non-nucleotide characters in a sequence
  percent_non_nucleotides <- function(seq, valid_chars) {
    chars <- strsplit(as.character(seq), "")[[1]]
    total_chars <- length(chars)
    if (total_chars == 0) return(1)  # Handle empty sequences as 100% non-nucleotides
    non_nuc_count <- sum(!toupper(chars) %in% valid_chars)
    return(non_nuc_count / total_chars)
  }
  
  # Apply filtering based on threshold
  filtered <- dnaset[
    sapply(dnaset, function(seq) percent_non_nucleotides(seq, valid_nucleotides) < non_nucleotide_threshold)
  ]
  
  return(filtered)


}
