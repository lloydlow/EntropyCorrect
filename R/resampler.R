#' resampler function
#'
#' This function is one of the functions that makes EntropyMSA works.
#' @param vect The input for entropy calculation in vector form.
#' @param randomise This parameter allows users to try other corrected entropy
#'   values. By default, set.seed is used in entropy correction and hence, same
#'   values are produced, which made it easier for the empiral study to determine
#'   the effects of correction in various sample and population sizes.
#' @keywords entropy
#' @export
#' @examples
#' #see EntropyMSA function
resampler <- function(vect,randomise){
  x.axis <- c(); y.axis <- c()
  for (j in 1:length(vect)){
    n <- j + 9 #start min resample size is 10
    x <- 1/n
    if (n <= length(vect)){
      x.axis <- c(x.axis,x)
      if (randomise == FALSE) {
        set.seed(j) #this is needed for empirical study
      }
      holder <- sample(vect,n,replace=FALSE)
      y.axis <- c(y.axis,entropy_vect(holder))
    }
  }
  fit <- lm(y.axis ~ x.axis)
  finalList <- list(x=x.axis,y=y.axis,intercept=fit[[1]][[1]],
                    slope=fit[[1]][[2]],R2=summary(fit)$r.squared,fit=fit)
}