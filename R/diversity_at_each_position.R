
library(ape)
library(phangorn)
library(pegas)
library(Biostrings)
library(strex) # for str_elem
library(rjags)


#' Function to calculate diversity at each position in an alignment
#'
#' This function checks the diversity contribution of each variant at each position.
#' If frequency of minor variants at a position is less than the sequencing error threshold,
#' that position contributes a diversity score of zero
#'
#' @param dnaset. DNAStringSet
#' 
#' @param errorthreshold. Sequencing error threshold. 
#' Usually 0.01 for NGS, 0 for prominent outgrowth viruses
#' 
#' @param valid_nucleotide. nucleotide characters to consider when calculating diversity.
#' internal gaps ("-") in the alignment are allowed.
#'
#' @return Diversity score at each position of the alignment
#' 
#' @examples
#' div_list <- diversity_at_each_position(dna_bin, 0, c("A","C","G","T","-"))
#'
#' @export
diversity_at_each_position <- function(dnaset, 
                                       errorthreshold = 0, 
                                       valid_nucleotides = c("A", "C", "G", "T", "-")){
  if (!inherits(dnaset, "DNAStringSet")) stop("Input must be a DNAStringSet")
  
  # width of the alignment
  width = unique(width(dnaset))
  
  # diversity at each position of the alignment
  d_accross = sapply(1:width, function(j){
    
    # valid_nucleotides at position j
    posj = str_elem(as.character(dnaset),j)
    
    # frequency of each variant at position j
    freqj = sapply(valid_nucleotides, function(x){sum(posj == x)/length(posj)})
    
    # frequency of major variant at position j
    freqj_max = freqj[freqj==max(freqj)]
    
    # check if sum of the frequency of minor valid_nucleotides is above the cutoff xc
    indicator_variable = ifelse( (1-freqj_max) > errorthreshold, 1, 0)
    
    # diversity contribution of each variant at position j
    d = sum(sapply(freqj, function(x){ x * (1-x)}))
    
    # result at position j
    resultj = indicator_variable*d
    
  })
}
