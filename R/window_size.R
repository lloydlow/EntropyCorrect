#' window_size function
#'
#' This function is one of the functions that makes EntropyMSA works.
#' @param MSA3 The input is a multiple sequence alignment in FASTA format.
#' @param n By default, n = 1, which means only a single character per alignment
#'   column per sequence. For example, users may change to n = 9 for nonamer 
#'   scoring instead of monomer.
#' @keywords entropy
#' @export
#' @examples
#' #see EntropyMSA function
window_size <- function(MSA3,n){
  #character column for each column
  column <- c()
  #numeric vector to hold each column H bit
  H <- c()
  
  for (o in seq_len(length(MSA3[[1]])-n+1)){
    for (q in 1:length(MSA3)){
      column <- c(column,paste(MSA3[[q]][seq(o,o+n-1)],sep="",collapse=""))
    }
    
    #entrophy calc here
    column.uniq <- unique(column)
    #holder for Pi in the Shannon formula
    p <- c()
    #holder for each term Pi*log2Pi to be summed
    #in the Shannon formula
    h <- c()
    
    for (k in 1:length(column.uniq)){
      select <- column == column.uniq[k]
      column.selc <- column[select]
      p[k] <- length(column.selc)/length(column)
      h[k] <- p[k]*log2(p[k])
    }
    H <- c(H,-sum(h))
    column <- c()
  }
  H
  #finalList <- list(H=H,MSA=MSA3)
}