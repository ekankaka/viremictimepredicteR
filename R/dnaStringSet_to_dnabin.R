
library(Biostrings)

#' Function to convert DNAStringSet objects to DNAbin 
#'
#' @param dnaset. DNAStringset object
#'
#' @return Returns DNAbin version of DNAStringSet
#'
#' @examples
#' dnabin = dnaStringSet_to_dnabin(notgappy)
#'
#' @export
dnaStringSet_to_dnabin <- function(dnaset) {
  # Convert DNAStringSet to a character vector
  seq_names <- names(dnaset)
  char_seqs <- as.character(dnaset)
  
  # Convert to list of character vectors (one per sequence)
  seq_list <- lapply(char_seqs, function(seq) strsplit(seq, "")[[1]])
  names(seq_list) <- seq_names
  
  # Convert to DNAbin
  dnabin_obj <- as.DNAbin(seq_list)
  
  return(dnabin_obj)
}

