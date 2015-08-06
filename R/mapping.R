#' Map: maps gDNA and mRNA to contigs
#'
#' Convert converts FASTQ input to FASTA format using perl Bio::SeqIO.
#' Runs \code{BLAT -fastMap} in fast mode and outputs the matches.
#' Matches are redundant ie. Reads are matched to more than one contig.
#'
#' @param blatOutputFile the output file to store results of blat
#'
#' @return a data.frame containing the readID the match and the blat information
#'
#' @export
mapper <- function(blatOutputFile){
    bof <- setNames(data.table::fread(blatOutputFile), c("query", "subject", "qstart", "qend", "sstart", "send"))
    bof$width = bof$send - bof$sstart
    bof$newStart = ifelse(bof$width > 0, bof$sstart, bof$send)
    bof$newEnd   = ifelse(bof$width < 0, bof$send, bof$sstart)
    bof
}
