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
    lapply(function(contig){
        data.frame(
            start      = toString(contig) %>% regexpr("^-*", .) %>% attributes %$% match.length,
            end        = length(contig) - toString(contig) %>% regexpr("-*$", .) %>% attributes %$% match.length,
            trueLength = gsub('-', "", toString(contig)) %>% nchar) #not really the true length, somehow the pAss pipeline either muscle or something else throws away the head and ends of the contigs good
    }) %>%
    do.call(rbind,.)

    cat(    sprintf("%s contigs processed", contigDF %>% nrow)  )

    contigDF$name = rownames(contigDF)
    contigDF                              %>%
    filter(trueLength  >=minIntervalSize) %$%
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

#' Extends MAX diversity region to capture mRNAs mapping to ends of reads when executing 100 percent identity read matches

#'
#' \code{extendRange} to allow mapping using Blat. Ie. Full map to mRNAs which touch the ends of the contigs
#'
#' Function builds the range
#' @param spanningGR the spanning Genomicranges obj
#' @param truncGR the truncated Genomicranges obj
#' @param insideGR the inside Genomicranges obj
#' @param mRNAsize the size of the mRNA
#' @param msaLength the MSA length
#'
#' @return list
#'
#' @export
extendRange <- function
(   spanningGR,     
    truncGR,        
    insideGR,       
    mRNAsize,       
    msaLength       
){    
inTheWindow = c()
    #SPANNING
    if(length(spanningGR) != 0){
        start(spanningGR) = sapply(start(spanningGR), function(x) ifelse (x - mRNAsize  <=  0,0,x - mRNAsize) )
        end(spanningGR)   = sapply(end(spanningGR), function(x)   ifelse (x + mRNAsize  >=  msaLength,msaLength,x+ mRNAsize))
        mcols(spanningGR)$type = 'spanning'
        inTheWindow = c(spanningGR)
    }

    #TRUNCATED
    if(length(truncGR) != 0){
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
        mcols(truncGR)$type = 'trunc'
        inTheWindow = c(inTheWindow, truncGR)
    }
        #Inside (not required)
        #    start(insideGR) = sapply(start(insideGR), function(x) ifelse (x - mRNAsize  <=  0,0,x - mRNAsize) )
        #    end(insideGR)   = sapply(end(insideGR), function(x)   ifelse (x + mRNAsize  >=  msaLength,msaLength,x+ mRNAsize))
        #Increase the span
####    ~
    if(length(insideGR) != 0){
        mcols(insideGR)$type = 'inside'
        inTheWindow = c(inTheWindow, insideGR)
    }
inTheWindow     #should include the identity as well.
}


#' extractLoc: extracts MSA information from the header of the FASTA files
#'
#' takes info in header of FASTA output from pAss pipeline and generates and location of the MAX Diversity region
#'
#' @param fileName the name of the file to inspect
#'
#' @return A list with the start and end location
#'
#' @export
extractLoc <- function(fileName){
    headerUP <- sprintf("head -n1 %s", fileName) %>% system(int = TRUE)
    mdUP     <- lapply(strsplit(strsplit(headerUP, " ## ")[[1]][2], " "), function(x) strsplit(x,":"))
    list(
         startLoc = as.integer(mdUP[[1]][[2]][[2]]),
         endLoc   = as.integer(mdUP[[1]][[3]][[2]])
     )
}
