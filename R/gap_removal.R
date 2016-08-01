#' gap_removal function
#'
#' This function is one of the functions that makes EntropyMSA works.
#' @param MSA2 The input is a multiple sequence alignment in FASTA format.
#' @param gap This is for gap removal in an alignment column. The default 
#'   accepts 100 percentage of gaps. Users can change the numeric value provided to
#'   this parameter from 0 to 1 to change the gaps that can be tolerated per
#'   alignment column. 0 for 0 percentage gaps, 1 for 100 percentage gaps.
#' @keywords entropy
#' @export
#' @examples
#' #see EntropyMSA function
gap_removal <- function(MSA2,gap){
  #character column for each column
  column <- c()
  #holder for the subset align that meets the gap requirements
  align_x.axis <- c()
  for (i in 1:length(MSA2[[1]])){
    for (j in 1:length(MSA2)){
      column <- c(column,MSA2[[j]][[i]])
    }
    #filter out gaps and record which col in the align they are
    howManyGaps.selc <- column == "-"
    column.gapsvect <- column[howManyGaps.selc]    
    #eqn: number of gaps / total characters per align column
    if (length(column.gapsvect)/length(MSA2) >= gap) {
      column <- c()
    } else{align_x.axis <- c(align_x.axis,i)
           column <- c()}
  }
  
  fastaid <- names(MSA2)
  MSA3 <- c()
  MSA3 <- as.list(MSA3)
  
  #Looping through each fasta read and subset for gaps removal
  for (l in 1:length(fastaid)){
    MSA3[[l]] <- MSA2[[l]][align_x.axis]
  }
  names(MSA3) <- fastaid
  MSA3
}