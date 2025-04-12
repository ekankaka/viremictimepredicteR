
#' Function to trim left and right terminal codons containing non-nucleotide characters, 
#' if any.
#'
#' @param dnaset. DNAStringset object
#'
#'  @param gap_chars. Vector of non-nucleotide characters 
#'  to trim at terminal ends of alignment (e.g., "-". "n", "N")
#'  
#' @param gap_threshold. Codons that contain this number of gap_characters 
#' or more will be trimmed from the ends of the alignment
#'
#' @return Returns DNAStringSet object
#'
#' @examples
#' trimmed = trim_terminal_gappy_codons(dnaset)
#'
#' @export
trim_terminal_gappy_codons <- function(dnaset, gap_chars = c("-", "n", "N"), 
                                       gap_threshold = 2) {
  if (!inherits(dnaset, "DNAStringSet")) stop("Input must be a DNAStringSet")
  
  # Convert to character matrix
  seqs <- as.character(dnaset)
  mat <- do.call(rbind, strsplit(seqs, split = ""))
  
  aln_len <- ncol(mat)
  n_seq <- nrow(mat)
  
  # Function to determine if a codon is gappy
  is_gappy_codon <- function(codon) {
    sum(codon %in% gap_chars) >= gap_threshold
  }
  
  # Initialize trimming indices
  start <- 1
  end <- aln_len
  
  # Trim from start
  while (start + 2 <= end) {
    codon <- mat[, start:(start + 2), drop = FALSE]
    if (all(apply(codon, 1, is_gappy_codon))) {
      start <- start + 3
    } else {
      break
    }
  }
  
  # Trim from end
  while (end - 2 >= start) {
    codon <- mat[, (end - 2):end, drop = FALSE]
    if (all(apply(codon, 1, is_gappy_codon))) {
      end <- end - 3
    } else {
      break
    }
  }
  
  # Subset the matrix and reconstruct DNAStringSet
  trimmed_mat <- mat[, start:end, drop = FALSE]
  trimmed_seqs <- apply(trimmed_mat, 1, paste0, collapse = "")
  trimmed_set <- DNAStringSet(trimmed_seqs)
  names(trimmed_set) <- names(dnaset)
  
  return(trimmed_set)
}
