#' make.map.table takes the output from map function (mapping reads onto contigs) and 
#'
#' takes the table which generates the overlap output back to data.frame
#'
#' @param overlaps output from the findOverlaps function from IRanges
#'
#' @return a data.frame object with the names of contigs and the number of reads mapped to them
#'
#' @export
make.map.table <- function(contigs, reads){
    overlapDF = findOverlaps(reads, contigs)
    countDF            = overlapDF@subjectHits %>% table %>% data.frame %>% setNames(c("contigID", "Freq"))
    contigNames = as.character(seqnames(contigs)@values)[countDF$contigID]
    countDF$contigName = contigNames
    mcolDF = regionOfInterest %>% mcols %>% as.data.frame
    cbind(
    countDF,
    mcolDF[countDF$contigID, c("type", "trueLength")]
    )
}
