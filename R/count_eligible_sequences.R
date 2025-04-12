
#' This Function ensures the remaining sequence count is at least the required minimum 
#'
#' @param dnaset. DNAStringSet object
#'
#' @param min_eligible_count. Minimum number of eligible sequences to accept
#'
#' @return Returns TRUE or FALSE
#'
#' @examples
#' pass1 = good_alignment_width(dnaset = trimmed, min_alignment_width = 9)
#'
#' @export
count_eligible_sequences <- function(dnaset, min_eligible_count = 2){
  # assert input is DNAstringset object
  if (!inherits(dnaset, "DNAStringSet")) stop("Input must be a DNAStringSet")
  
  if (length(dnaset) >= min_eligible_count ){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

