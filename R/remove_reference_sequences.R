
#' Function to remove reference sequences 
#' with specified reference identifiers in the sequence headers, if any
#' 
#' @param dnaset. DNAStringSet object
#'
#' @param ref_identifiers. Vector of pattern(s) in headers that identify reference sequences
#'
#' @return Returns DNAStringSet object
#'
#' @examples
#' noRefs = remove_reference_sequences(dnaset)
#'
#' @export
remove_reference_sequences <- function(dnaset, ref_identifiers = c("Ref", "HXB2")){
  
  # assert that DNAstringset exists
  if (!inherits(dnaset, "DNAStringSet")) stop("Input must be a DNAStringSet")
  
  # paste together the patterns to search for in sequence headers
  patterns = paste(ref_identifiers, collapse = "|")
  
  # keep only sequences that do not have these patterns in their headers
  noRefs = dnaset[!grepl(patterns, names(dnaset))]
  
  return(noRefs)
  
}
