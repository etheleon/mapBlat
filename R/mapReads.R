#' mapReads
#'
#' mapReads maps reads to a region of interest. 
#' A cut off of 200 basepairs for inclusion of the contigs to be mapped, to facilitate for at least mapping of a 100 BP paired end read (need to verify).
#' the regions of interest can be extracted from the header of the fasta files generated from the pAss pipeline.
#'
#' @param kolist vector of KOs to iterate through
#' @param MSA path to MSA directory
#' @param roi region of interest in the form of a named list with names startLoc and endLoc.
#' @param mapDetails Is a filePath to the mapping file in (tabular blast format) planned feature doesnt work now. Use perl script output
#'
#' @export
mapReads <- function(ko, msa, roi, mapDetails){
    upMSA           <-  readDNAStringSet(sprintf("%s/%s.msa",msaDIR, ko))
    #Cleanse NAME
    upMSA@ranges@NAMES %<>%  stringr::str_extract("contig\\d+")
    #genomicRanges 
    upMSA_int       <-  buildIntervals(stringSet = upMSA, minIntervalSize = 200)
    contigsSpan     <-  spanningOrNot(upMSA_int, roi$startLoc, roi$endLoc)
    #runPerlMapper
    blatDF          <-  mapper(blatOutputFile = sprintf("~/reDiamond/out/blatter/trimmed_blatoutput_%s", ko))

    blatRanges <- GRanges(seqnames=blatDF$subject, ranges=IRanges(blatDF$newStart, blatDF$newEnd), read=blatDF$query)

#Overlap
    cDNARanges <- mcols(blatRanges)@listData[[1]] %>% as.character %>% grepl(pattern="cDNA", x=.) %>% blatRanges[.]
    gDNARanges <- mcols(blatRanges)@listData[[1]] %>% as.character %>% grepl(pattern="gDNA", x=.) %>% blatRanges[.]

    ROI <- extendRange(contigsSpan$spanning, contigsSpan$trunc, contigsSpan$inside,99,width(upMSA)[1])

    countDF_cDNA <-     make.map.table(contigs = ROI, reads = cDNARanges, datype="cDNA")
    countDF_cDNA %<>%   mutate(ko = ko)

    countDF_gDNA <-     make.map.table(contigs = ROI, reads = gDNARanges, datype="gDNA")
    countDF_gDNA %<>%   mutate(ko = ko)


    combinedDF <- merge(
          countDF_cDNA %>% select(-type, -trueLength, -contigID),
          countDF_gDNA %>% select(-type, -trueLength, -contigID),
          by = c("contigName","ko"),
          suffixes = c("_cDNA", "_gDNA"),
          all=T
    )
    combinedDF[is.na(combinedDF)] = 0
    lengthDF = data.frame(contigName = seqnames(upMSA_int)@values, trueLength = mcols(upMSA_int)$trueLength)
    merge(combinedDF, lengthDF, by="contigName", all.x=T)
#}) %>% do.call(rbind,.)
}
