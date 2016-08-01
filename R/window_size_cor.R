#' window_size_cor function
#'
#' This function is one of the functions that makes EntropyMSA works.
#' @param MSA3 The input is a multiple sequence alignment in FASTA format.
#' @param n By default, n = 1, which means only a single character per alignment
#'   column per sequence. For example, users may change to n = 9 for nonamer 
#'   scoring instead of monomer.
#' @param randomise This parameter allows users to try other corrected entropy
#'   values. By default, set.seed is used in entropy correction and hence, same
#'   values are produced, which made it easier for the empiral study to determine
#'   the effects of correction in various sample and population sizes.
#' @keywords entropy
#' @export
#' @examples
#' #see EntropyMSA function
window_size_cor <- function(MSA3,n,randomise){
  #character column for each column
  column <- c()
  #numeric vector to hold each column H bit
  #H <- c()
  Hcor <- c();slope <- c();r2 <- c() 
  for (o in seq_len(length(MSA3[[1]])-n+1)){
    for (q in 1:length(MSA3)){
      column <- c(column,paste(MSA3[[q]][seq(o,o+n-1)],sep="",collapse=""))
    }
    #H <- c(H,entropy_vect(column))
    H_list_result <- resampler(column,randomise)
    Hcor <- c(Hcor,H_list_result$intercept)
    slope <- c(slope,H_list_result$slope)
    r2 <- c(r2,H_list_result$R2)
    column <- c()
  }
  #H
  Hcor
  #finalList <- list(H=Hcor,slope=slope,R2=r2)
}