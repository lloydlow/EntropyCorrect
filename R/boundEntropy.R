#' boundEntropy function
#'
#' This function is one of the functions that makes EntropyMSA works.
#' @param H The is the calculated H, entropy value, to be bounded by the maximum entropy
#'  and minimum entropy possible.
#' @param n By default, n = 1, which means only a single character per alignment
#'   column per sequence. For example, users may change to n = 9 for nonamer 
#'   scoring instead of monomer.
#' @param totalCharacters This depends on whether the standard 20 amino acids plus
#'   a gap scoring is used or not. Users may set this to 20 if they wish to ignore
#'   gap character.
#' @keywords entropy
#' @export
#' @examples
#' #see EntropyMSA function
boundEntropy <- function(H,n,totalCharacters = 21){
  #assume 21 characters as the default
  max <- log2(totalCharacters^n)
  for (i in length(H)){
    if (H[i] > max) H[i] <- max
    if (H[i] < 0) H[i] <- 0
  }
  H
}