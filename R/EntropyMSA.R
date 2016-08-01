#' An entropy correction function
#'
#' This function estimates entropy in an infite size population from a sample.
#' @param MSA The input is a multiple sequence alignment in FASTA format.
#' @param gap This is for gap removal in an alignment column. The default 
#'   accepts 100 percentage of gaps. Users can change the numeric value provided to
#'   this parameter from 0 to 1 to change the gaps that can be tolerated per
#'   alignment column. 0 for 0 percentage gaps, 1 for 100 percentage gaps.
#' @param n By default, n = 1, which means only a single character per alignment
#'   column per sequence. For example, users may change to n = 9 for nonamer 
#'   scoring instead of monomer. 
#' @param no.gap.firstseq This allows the user to set gap removal based on a
#'   reference sequence that must be placed as the first sequence in the
#'   multiple alignment file.
#' @param entropy.correct By default, there is no entropy correction. If users
#'   wish to perform a correction, they should set entropy.correct = TRUE.
#' @param totalCharacters This depends on whether the standard 20 amino acids plus
#'   a gap scoring is used or not. Users may set this to 20 if they wish to ignore
#'   gap character.
#' @param write.align After filtering for gaps either with no.gap.firstseq and 
#'   gap parameters, the actual multiple alignment used differ. Users may choose
#'   to output the final alignment to find out which alignment columns correspond
#'   to the entropy values.
#' @param randomise This parameter allows users to try other corrected entropy
#'   values. By default, set.seed is used in entropy correction and hence, same
#'   values are produced, which made it easier for the empiral study to determine
#'   the effects of correction in various sample and population sizes.
#' @keywords entropy
#' @export
#' @examples
#' #for monomer and a comparison of corrected versus uncorrected
#' file1 <- system.file("extdata","exampleMSA1.fasta",package="EntropyCorrect")
#' H.uncor <- EntropyMSA(file1, no.gap.firstseq = TRUE)
#' H.cor1 <- EntropyMSA(file1, no.gap.firstseq = TRUE, entropy.correct = TRUE)
#' H.cor2 <- EntropyMSA(file1, no.gap.firstseq = TRUE, entropy.correct = TRUE, randomise = TRUE)
#' layout(matrix(c(1,2),2,1,byrow = TRUE), widths = 7, 
#' heights = c(8, 2), respect = FALSE)
#' par(mar = c(0, 4, 3, 2))
#' barplot(H.uncor,border="white", ylab = "Entropy",space=0,yaxt="n",ylim = c(0,4))
#' axis(side=2, at=c(0,0.5, 1,1.5,2,2.5,3.0,3.5,4))
#' axis(side=3, at=(seq(0,length(H.uncor),by=5)-0.5),labels=seq(0,length(H.uncor),by=5))
#' lines(1:length(H.cor1),H.cor1,col="red")
#' #lines(1:length(H.cor1),H.cor1,col="black",lty=5) #alternative presentation
#' #lines(1:length(H.cor2),H.cor2,col="blue") #alternative line
#' box()
#' mtext("Alignment position", side=3, line=2)
#' par(mar = c(2, 4, 0, 2))
#' plot(c(0, length(H.cor1)), c(0, 100), type= "n", xlab = "", 
#' ylab = "",yaxt="n",xaxt="n",bty="n")
#' rect(22, 30, 101, 80,density=10)
#' rect(106, 30, 230, 80,density=20, angle=135)
#' lines(c(0,22),c(55,55),lwd=2)
#' lines(c(101,106),c(55,55),lwd=2)
#' lines(c(230,243),c(55,55),lwd=2)
#' text(65, 10, "N terminal")
#' text(165, 10, "C terminal")
#' #for nonamer and a comparison of corrected versus uncorrected
#' H9 <- EntropyMSA(file1, n = 9)
#' H9.cor <- EntropyMSA(file1, n = 9, entropy.correct = TRUE)
#' #for nonamer, a simple example to understand entropy calculation
#' file2 <- system.file("extdata","test.fasta",package="EntropyCorrect")
#' H9.new <- EntropyMSA(file2,n=9)
EntropyMSA <- function(MSA,gap = 1, n = 1, 
                       no.gap.firstseq = FALSE,
                       entropy.correct = FALSE,
                       totalCharacters = 21,
                       write.align = FALSE,
                       randomise = FALSE) {
  MSA2 <- no_gap_firstseq(MSA,no.gap.firstseq)
  MSA3 <- gap_removal(MSA2,gap)
  if (write.align == TRUE) writealn(MSA3,"alignment.faa")
  if (entropy.correct == FALSE) H <- window_size(MSA3,n)
  else H <- window_size_cor(MSA3,n,randomise)
  H <- boundEntropy(H,n,totalCharacters)
}