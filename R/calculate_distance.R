
#' Function to calculate pairwise distances from DNA sequences
#'
#' @param dnaset. DNAStringSet object
#' 
#' @param sequence_type. outgrowth or provirus (case-insensitive)
#' Used to determine error thresholds when calculating diversity metrics WFPS and WFPScodons
#' 
#' @param min_eligible_count. Minimum number of eligible sequences to accept
#' 
#' @return Returns a list of input parameters and result of different diversity metrics.
#' Currently supported diversity metrics include:
#' - mean pairwise distance from the "raw" model (rawMPD)
#' - mean pairwise distance from the "TN93" model (tn93MPD)
#' - nucleotide diversity from the "raw" model (rawPI)
#' - nucleotide diversity from the "TN93" model (tn93PI)
#' - weighted fraction of polymorphic sites (WFPS)
#' - WFPS at third codon positions only (WFPScodons)
#' @examples
#' dist <- calculate_distance(dnaset = mydnaset, min_eligible_count = 2)
#'
#' @export
calculate_distance <- function(dnaset, min_eligible_count = 2){
  
  if (!inherits(dnaset, "DNAStringSet")) stop("Input must be a DNAStringSet")
  
  pass2 = count_eligible_sequences(dnaset, min_eligible_count)
  if (pass2 == FALSE) stop("Remaining sequence count is less than the required minimum.")
  
  # convert to dnabin, on DNAbin
  dnabin = dnaStringSet_to_dnabin(dnaset)
  
  # Calculate raw mean pairwise distance, on DNAbin
  raw_distances <- dist.dna(dnabin, model = "raw")
  rawMPD <- mean(raw_distances)
  
  # Calculate TN93 mean pairwise distance, on DNAbin
  tn93_distances <- dist.dna(dnabin, model = "TN93")
  tn93MPD <- mean(tn93_distances)
  
  # Calculate nucleotide diversity (pi) from raw p-distances, on DNAbin
  rawPI <- nuc.div(dnabin, model = "raw")
  
  # Calculate nucleotide diversity (pi) from TN93 distances, on DNAbin
  tn93PI <- nuc.div(dnabin, model = "TN93")
  
  #calculate weighted fraction of polymorphic sites (WFPS or APD), on DNAStringSet
  d = diversity_at_each_position(dnaset = dnaset,
                                 errorthreshold = threshold,
                                 valid_nucleotides = c("A", "C", "G", "T","-") )
  d = unlist(d)
  WFPS = mean(d)
  
  # calculate WFPS at 3rd codon positions (WFPScodons or APD3)
  d3 <- d[seq(3,length(d),3)]
  WFPScodons = mean(d3)
  
  # return result
  result = c(rawMPD = rawMPD, tn93MPD = tn93MPD, rawPI = rawPI, tn93PI = tn93PI,
             WFPS = WFPS, WFPScodons = WFPScodons)
  
  return(result)

}
