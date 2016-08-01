#' entropy_vect function
#'
#' This function is one of the functions that makes EntropyMSA works.
#' @param vect The input for entropy calculation in vector form.
#' @keywords entropy
#' @export
#' @examples
#' #see EntropyMSA function
entropy_vect <- function(vect){
  H <- c();p <- c();h <- c()
  vect.uniq <- unique(vect)
  for (k in 1:length(vect.uniq)){
    select <- vect == vect.uniq[k]
    vect.selc <- vect[select]
    p[k] <- length(vect.selc)/length(vect)
    h[k] <- p[k]*log2(p[k])
  }
  H <- c(H,-sum(h))
}