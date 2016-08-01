#' writealn function
#'
#' This function is one of the functions that makes EntropyMSA works.
#' @param MSA3 The input is a multiple sequence alignment in FASTA format.
#' @param filename This is the filename of the output alignment.
#' @keywords entropy
#' @export
#' @examples
#' #see EntropyMSA function
writealn <- function(MSA3,filename){
  if(file.exists(filename)) unlink(filename)
  fastaid <- names(MSA3)
  for (i in 1:length(fastaid)){
    eachFastalist <- MSA3[i]
    eachName <- names(eachFastalist)
    eachName <- c(">",eachName)
    eachName <- paste(eachName,sep="",collapse = "")
    Tfile <- file(filename, "a")
    cat(eachName,"\n", file = Tfile)
    seq <- unname(unlist(eachFastalist))
    seq <- paste(seq,sep="",collapse = "")
    cat(seq,"\n", file = Tfile)
    close(Tfile)
  }
}