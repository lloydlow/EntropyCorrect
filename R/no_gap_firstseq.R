#' no_gap_firstseq function
#'
#' This function is one of the functions that makes EntropyMSA works.
#' @param MSAfile The input is a multiple sequence alignment in FASTA format.
#' @param no.gap.firstseq This allows the user to set gap removal based on a
#'   reference sequence that must be placed as the first sequence in the
#'   multiple alignment file.
#' @keywords entropy
#' @export
#' @examples
#' #see EntropyMSA function
no_gap_firstseq <- function(MSAfile,no.gap.firstseq){
  
  #Implement my own read fasta
  MSA <- fasta2list(MSAfile)
  
  align_firstseqx.axis <- c()
  
  firstseq <- MSA[[1]]
  
  if (no.gap.firstseq == TRUE){ 
    for (element_aa in 1:length(firstseq)){
      if (firstseq[element_aa] != "-") {
        align_firstseqx.axis <- c(align_firstseqx.axis,element_aa)
      }
    }
  }
  
  if (no.gap.firstseq == FALSE){
    for (element_aa in 1:length(firstseq)){
      align_firstseqx.axis <- c(align_firstseqx.axis,element_aa)
    }
  }
  
  fastaid <- names(MSA)
  MSA2 <- c()
  MSA2 <- as.list(MSA2)
  
  #Looping through each fasta read and subset for gaps removal
  for (l in 1:length(fastaid)){
    MSA2[[l]] <- MSA[[l]][align_firstseqx.axis]
  }
  names(MSA2) <- fastaid
  MSA2
}