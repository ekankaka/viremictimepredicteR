
#' Function checks if the alignment width is at least the required minimum
#'
#' @param dnaset. DNAStringSet object
#'
#' @param min_alignment_width. Minimum alignment width to accept
#'
#' @return Returns TRUE or FALSE
#'
#' @examples
#' pass1 = good_alignment_width(dnaset = trimmed, min_alignment_width = 9)
#'
#' @export
check_alignment_width <- function(dnaset, min_alignment_width = 9){
  
  if (!inherits(dnaset, "DNAStringSet")) stop("Input must be a DNAStringSet")
  
  if (!(all(width(dnaset) == width(dnaset)[1]))) stop("Sequences are not aligned")
  
  if (unique(width(dnaset)) >= min_alignment_width ){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
