#' fasta2list function
#'
#' This function is one of the functions that makes EntropyMSA works.
#' @param MSAfile The input is a multiple sequence alignment in FASTA format.
#' @keywords entropy
#' @export
#' @examples
#' #see EntropyMSA function
fasta2list <- function(MSAfile){
  if (missing(MSAfile)) stop("File is missing.")
  lines <- readLines(MSAfile)
  fastaheader <- grep(">", lines)
  if (length(fastaheader) == 0) stop("File is not in fasta format.")
  id <- sub("^>(\\S+).*$", "\\1", lines[fastaheader])
  totalseq <- length(id)
  start <- fastaheader + 1;   end <- fastaheader - 1
  end <- c(end[-1], length(lines))
  seq <- sapply(seq_len(totalseq), function(i) {
    paste(lines[start[i]:end[i]], collapse = "")
  })
  seq <- gsub("\\s", "", seq)
  seq <- toupper(seq)
  seq <- strsplit(seq, split = "")
  names(seq) <- id
  return(seq)
}