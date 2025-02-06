
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
#' @param fasta. DNAbin object
#' @param errorthreshold. Sequencing error threshold. Usually 0.01 for NGS, 0 for prominent outgrowth viruses
#' @param variants. expected nucleotide characters.
#'
#' @return Diversity score at each position of the alignment
#' @examples
#' div_list <- diversity_at_each_position(dna_bin, 0, c("A","C","G","T","-"))
#'
#' @export
diversity_at_each_position <- function(fasta, errorthreshold, variants){
  # width of the alignment
  width = dim(fasta)[2]

  # diversity at each position of the alignment
  d_accross = sapply(1:width, function(j){

    # variants at position j
    posj = str_elem(as.character(fasta),j)

    # frequency of each variant at position j
    freqj = sapply(variants, function(x){sum(posj == x)/length(posj)})

    # frequency of major variant at position j
    freqj_max = freqj[freqj==max(freqj)]

    # check if sum of the frequency of minor variants is above the cutoff xc
    indicator_variable = ifelse( (1-freqj_max) > errorthreshold, 1, 0)

    # diversity contribution of each variant at position j
    d = sum(sapply(freqj, function(x){ x * (1-x)}))

    # result at position j
    resultj = indicator_variable*d

  })
}
