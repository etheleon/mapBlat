#' Map: maps gDNA and mRNA to contigs
#'
#' Convert converts FASTQ input to FASTA format using perl Bio::SeqIO.
#' Runs BLAT in fast mode and outputs the matches
#'
#' @param fq vector of FASTQ files; read1, read2 or single read
#' @param fa FASTA file to store
#' @param blatOutputFile the output file to store results of blat
#' @param fasta file containing contigs
#'
#' @return a data.frame containing the readID the match and the blat information
#'
#' @export
map <- function(fq, fa, blatOutputFile, contigs){
    #Check if fq and fa have reasonable file extension suffixes
    isFastq <- sum(grepl("(fq|fastq)$", fq, TRUE)) == length(fq)
    isFasta <- sum(grepl("(fa|fasta)$", fa, TRUE)) == length(fa)

    if(! (isFastq & isFasta))
        stop("Fasta and Fastq files do not have the right suffixes")

    if(file.exists(fa)){
        message(sprintf("An existing version of %s exists.
                        It will be removed and overwritten", fa))
        file.remove(fa)
    }

    message("Running Bioperl")
    sprintf(
        "perl -E 'use aliased 'Bio::SeqIO';
        my $i=SeqIO->new(-file=>qq|%s|);
        my $o=SeqIO->new(-file=>qq|>>%s|);
        while(my $n=$i->next_seq){$o->write_seq($n)}'",
    fq, fa) %>%
    sapply(function(x) system(x))

    message("Running blat")
    sprintf(
    "blat -fastMap -out=blast8 %s %s %s",
    contigs, fa, blatOutputFile) %>%
    system

    message("Reading in blatOutput")

    blatOutputFile                                         %>%
    read.table(sep="\t", h=F, comment.char="")             %>%
    setNames(c("query", "subject", "identity",
               "alignLen", "mismatch", "gapOpenings",
               "qstart", "qend", "sstart", "send",
               "evalue", "bitscore"))                      %>%

    mutate(
        newStart = ifelse(send - sstart > 0, sstart, send),
        newEnd   = ifelse(send - sstart > 0, send,sstart)
    )                                                      %>%

    select(query:subject, qstart:send, bitscore:newEnd)
}
