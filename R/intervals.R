#' Builds intervals based on a stringSet Object
#'
#' \code{buildIntervals} returns GRanges object with the ranges of contigs spanning the GLOBAL MSA
#'
#' Takes in a DNAStringSet object from \code{Biostrings} generated from reading in the fasta MSA of contigs generated from pAss pipeline
#'
#' @param stringSet string Object from Biostrings
#' @param minIntervalSize the minimal size of the contig
#'
#' @return A GRanges object from \code{GenomicRanges} showing intervals
#'
#' @export
buildIntervals <- function(stringSet, minIntervalSize){
    contigDF = stringSet %>%
    lapply(function(contig) {
        data.frame(
            start      = toString(contig) %>% regexpr("^-*", .) %>% attributes %$% match.length,
            end        = length(contig) - toString(contig) %>% regexpr("-*$", .) %>% attributes %$% match.length,
            trueLength = gsub('-', "", toString(contig)) %>% nchar)
    }) %>%
    do.call(rbind,.)

    sprintf("%s contigs processed", contigDF %>% nrow)

    contigDF                            %>%
    filter(trueLength>=minIntervalSize) %>%
    mutate(name = rownames(.))          %$%
    GRanges(seqnames   = name,
            ranges     = IRanges(start, end),
            trueLength = trueLength)
}

#' Determines nature of contigs wrt to Max Diversity Region; is this spanning, non-spanning, entirely captured inside
#'
#' \code{spanningOrNot} Figures out the way the contigs map to the max diversity regions
#'
#' This function is part of a series 
#' @param gr A GRanges Object
#' @param startLOC the start location for the region of interest
#' @param endLOC the end location for the region of interest
#'
#' @return Returns a list of GRanges Objects which are 1. spanning, 2. non-spanning or entirely inside
#'
#' @export
spanningOrNot <- function(gr, startLOC, endLOC){

    #spanning
    spanning       =  gr[which(start(gr) <= startLOC  & end(gr) >= endLOC)]

    #non-spanning
    nsAll          =  gr[which(!1:length(gr) %in% which(start(gr) <= startLOC  & end(gr) >= endLOC))]    #inclusive of things outside entirely

    #Fragments
    nsInside       =  gr[which(start(gr) > startLOC  & end(gr) > endLOC)]

    #Truncated
    nsTruncated    = nsAll[which(seqnames(nsAll)@values %in% seqnames(nsAll)@values[!seqnames(nsAll)@values %in% seqnames(nsInside)@values])]
    nsTruncated    =  nsTruncated[-which(end(nsTruncated) < startLOC | start(nsTruncated) > endLOC)]    #remove things which 

    mcols(nsTruncated)$spillover = ifelse(end(nsTruncated)<endLOC, "left", "right")

    list(
        spanning = spanning,
        nsAll    = nsAll,
        inside   = nsInside,
        trunc    = nsTruncated
        )
}

#' Extends MAX diversity region to capture mRNAs mapping to ends of reads when executing 100% percent identity read matches
#'
#' extendRange to allow mapping using Blat. Ie. 100% map to mRNAs which touch the ends of the contigs
#'
#' @param spanningGR
#' @param truncGR
#' @param insideGR
#' @param mRNAsize
#' @param msaLength
#'
#' @return list 
#'
#' @export
extendRange <- function(spanningGR,truncGR,insideGR, mRNAsize, msaLength){    #the mRNA size must be after taking away the 1st read
    #spanning
    start(spanningGR) = sapply(start(spanningGR), function(x) ifelse (x - mRNAsize  <=  0,0,x - mRNAsize) )
    end(spanningGR)   = sapply(end(spanningGR), function(x)   ifelse (x + mRNAsize  >=  msaLength,msaLength,x+ mRNAsize))
#Truncated
    #which end
    startEND = do.call(rbind,mapply(function(spillover,startLOC,endLOC){
        if(spillover == 'left'){
            data.frame(newstart    = ifelse(startLOC-mRNAsize < 0, 0,startLOC-mRNAsize), newend=endLOC)
        }else{
            data.frame(newstart    =startLOC, newend=ifelse(endLOC+mRNAsize> msaLength, msaLength,endLOC+mRNAsize))
        }
    },
       spillover = mcols(truncGR)$spillover,
       startLOC  = start(truncGR),
       endLOC    = end(truncGR),
       SIMPLIFY=F))
    start(truncGR) = startEND$newstart
    end(truncGR) = startEND$newend
mcols(truncGR)@listData = list(trueLength = mcols(truncGR)@listData$trueLength)
#Inside (not required)
#    start(insideGR) = sapply(start(insideGR), function(x) ifelse (x - mRNAsize  <=  0,0,x - mRNAsize) )
#    end(insideGR)   = sapply(end(insideGR), function(x)   ifelse (x + mRNAsize  >=  msaLength,msaLength,x+ mRNAsize))
#Increase the span
####    ~
mcols(spanningGR)$type = 'spanning'
mcols(truncGR)$type = 'trunc'
mcols(insideGR)$type = 'inside'
inTheWindow = c(spanningGR,truncGR,insideGR)    #should include the identity as well.
}
